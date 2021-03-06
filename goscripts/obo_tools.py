#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@author: Pieter Moris
'''

import copy
import os
import re


class goTerm:
    """
    GO term object.

    Stores the ID, name and domain of the GO term and contains dictionaries for child and parent nodes.

    Attributes
    ----------
    id : str
        The identifier of the GO term.
    alt_id : str
        Optional tag for an alternative id.
    name : str
        The GO term name.
    namespace : str
        The domain of the GO term (Cellular Component, Molecular Function or Biological Process).
    parents : set of str
        The parent terms of the GO term, as indicated by the `is_a` relationship.
    children : set of str
        The child terms of the GO term, derived from other GO terms after a complete OBO file is processed initially.


    # https://stackoverflow.com/questions/1336791/dictionary-vs-object-which-is-more-efficient-and-why
    # https://stackoverflow.com/questions/3489071/in-python-when-to-use-a-dictionary-list-or-set
    # When you want to store some values which you'll be iterating over, Python's list constructs are slightly faster.
    # However, if you'll be storing (unique) values in order to check for their existence, then sets are significantly faster.
    '''
    """

    goCount = 0

    __slots__ = ('id', 'name', 'alt_id', 'namespace', 'children', 'parents',
                 'recursive_children', 'recursive_parents', 'depth')

    def __init__(self, GOid):
        self.id = GOid
        self.alt_id = []
        self.name = None
        self.namespace = None
        self.children = set()
        self.parents = set()
        self.recursive_children = set()
        self.recursive_parents = set()
        self.depth = None

        goTerm.goCount += 1


def importOBO(path, ignore_part_of):
    """
    Imports an OBO file and generates a dictionary containing an goTerm object for each GO term.

    Every entry that is referred to via either "is_a" or "relationship: part_of" is considered
    a parent of the referring entry.

    Parameters
    ----------
    path : str
        The path to the file.
    ignore_part_of : bool
        A boolean indicating whether or not the "part_of" relationship should be ignored.

    Returns
    -------
    dict of goTerm objects
        Keys are of the format `GO-0000001` and map to goTerm objects..

    Possible improvements:
        Check for `is_obsolete` and `replaced_by`, although the replacement term should be in OBO file as an entry.
    """

    GOdict = {}

    path = os.path.abspath(path)
    with open(path, 'r') as oboFile:

        # find re pattern to match '[Entry]'
        entryPattern = re.compile('^\[.+\]')
        validEntry = False

        for line in oboFile:

            # Only parse entries preceded by [Entry], not [Typedef]
            if entryPattern.search(line):
                if 'Term' in line:
                    validEntry = True
                else:
                    validEntry = False

            # if [Entry] was encountered previously, parse annotation
            elif validEntry:
                if line.startswith('id'):
                    # Store ID for lookup of other attributes in next lines
                    GOid = line.split(': ')[1].rstrip()

                    if GOid not in GOdict:
                        # check if ID is already stored as a key in dictionary and if not,
                        # create a new GOid object as the value for this key
                        GOdict[GOid] = goTerm(GOid)

                # Store all the other attributes for the current term as an object attribute
                elif line.startswith('name:'):
                    GOdict[GOid].name = line.split(': ')[1].rstrip()
                elif line.startswith('namespace:'):
                    GOdict[GOid].namespace = line.split(': ')[1].rstrip()
                elif line.startswith('alt_id:'):
                    GOdict[GOid].alt_id.append(line.split(': ')[1].rstrip())
                elif line.startswith('is_a:'):
                    GOdict[GOid].parents.add(line.split()[1].rstrip())
                elif line.startswith('relationship: part_of'):
                    if not ignore_part_of:
                        GOdict[GOid].parents.add(line.split()[2].rstrip())

    print('Retrieved', len(GOdict), 'GO terms from', path, '\n')

    # create secondary IDs
    print('Adding secondary GO identifiers...\n')
    secondary_id_dict = {}
    for term_id, term in GOdict.items():
        try:
            # only do this for GO terms which have an alternative ID attribute
            # NOTE: not strictly necessary since alt_id is defined as an empty
            # list upon creating new GOterm objects
            for alt_id in term.alt_id:
                # check if alternative id is already present in GO term dictionary
                if alt_id not in GOdict:
                    # if not, add it as its own instance
                    secondary_id_dict[alt_id] = copy.deepcopy(term)
                    secondary_id_dict[alt_id].id = alt_id
                    secondary_id_dict[alt_id].alt_id.append(term_id)
                    secondary_id_dict[alt_id].alt_id.remove(alt_id)
                else:
                    print('Alternative ID term was already present.\n')
        except AttributeError:
            continue

    # Merge alternative id's with main dictionary
    GOdict = {**GOdict, **secondary_id_dict}

    return GOdict


def filterOnNamespace(GOdict, namespace):
    """
    Reduces the dictionary of goTerm objects to those belonging to a specific namespace.

    Parameters
    ----------
    GOdict : dict
        A dictionary of GO objects generated by importOBO().
        Keys are of the format `GO-0000001` and map to goTerm objects.
    namespace : str
        The namespace to restrict the GO dictionary and enrichment test to. E.g. biological_process, cellular_component
        or molecular_function.

    Returns
    -------
    dict
        A filtered dictionary of GO objects all belonging to the namespace.
    """

    filteredGOdict = {GOid: GOobj for GOid, GOobj in GOdict.items() if GOobj.namespace == namespace}

    if not filteredGOdict:
        print('Namespace', namespace, 'was not found in the obo file. Using all annotations in obo file instead.\n')
        return GOdict

    print('Found', len(filteredGOdict), 'GO terms belonging to', namespace + '.\n')

    return filteredGOdict


def set_namespace_root(namespace):
    """
    Stores the GO ID for the root of the selected namespace.

    Parameters
    ----------
    namespace : str
        A string containing the desired namespace. E.g. biological_process, cellular_component
        or molecular_function.

    Returns
    -------
    list
        The list of GO ID's of the root terms of the selected namespace.
    """
    if namespace == 'biological_process':
        namespace_list = ['GO:0008150']
    elif namespace == 'cellular_component':
        namespace_list = ['GO:0005575']
    elif namespace == 'molecular_function':
        namespace_list = ['GO:0003674']
    else:
        namespace_list = ['GO:0008150', 'GO:0005575', 'GO:0003674']

    return namespace_list


def buildGOtree(GOdict, root_nodes):
    """
    Generates the entire GO tree's parent structure by walking through the hierarchy of each GO entry.

    Performs four main functions:
    A) assign all higher order ancestors (recursive parents) for each GO object
    B) assign immediate children
    C) assign all recursively found children
    D) assign depth to each node.

    Parameters
    ----------
    GOdict : dict
        A dictionary of GO objects generated by importOBO().
        Keys are of the format `GO-0000001` and map to OBO objects.
    root_nodes : list
        A list of nodes that lie at the root of the gene ontologies. These form the starting point for the
        recursive function that walks through the directed acyclic graph and assigns depths to each GO term/node.

    Returns
    -------
    None
        The input GO dictionary will be updated in place.
        Each term object's parent attributes now trace back over the full tree hierarchy.
    """

    # NOTE: The process CANNOT be sped up by working upwards from just the GO terms
    # that are present in the subset association dictionary, because by doing this
    # not all of their child terms will be found. Yet these child terms should also be
    # counted while testing their parent terms.
    # subsetGOids = {GOid for gene, GOids in gafSubset.items() for GOid in GOids}
    # for GOid in subsetGOids:
    #     parentSet = set()
    #     propagateParents(GOid, GOid, GOdict, parentSet)
    #     GOdict[GOid].parents.update(parentSet)

    # Process each GO term in the GO dictionary to A) recursively find parents and...
    for GOid, GOobj in GOdict.items():
        # Define new set to store higher order parents (recursed)
        parentSet = set()
        # Call helper function to propagate through parents recursively
        propagateParents(GOid, GOid, GOdict, parentSet)
        # Update GO term's parents attribute to include all higher order ancestors
        GOobj.recursive_parents = parentSet

        # ... B) set children attribute of GO IDs
        for parent in GOobj.parents:
            try:
                GOdict[parent].children.add(GOid)
            except KeyError:
                print(f'WARNING: The .obo file provided {parent} as a parent of {GOid}, but {parent} itself '
                      f'was not found in the GO dictionary. Check if it exists in the .obo file.\nThis might '
                      f'be the result of namespace filtering and cross-namespace relations in the .obo file.\n')

    # C) After all parents have been found, for each GO ID, add it as a child for all of its parents/ancestors
    completeChildHierarchy(GOdict)

    # D) Assign depth to each node.
    # start from each of the three root nodes
    for i in root_nodes:
        if i in GOdict:
            assign_depth(GOdict.get(i), GOdict)

    return None


def propagateParents(currentTerm, baseGOid, GOdict, parentSet):
    """
    Propagates through the parent hierarchy of a provided GO term to create a set of all higher order parents.

    Each term's recursive_parents attribute will be filled with all recursively found parent terms.

    Parameters
    ----------
    currentTerm : str
        The GO id that is being visited.
    baseGOid : str
        The original GO term id for which the search for its parents was started.
    GOdict : dict
        A dictionary of GO objects generated by importOBO().
        Keys are of the format `GO-0000001` and map to goTerm objects.
    parentSet : set
        An, initially, empty set that gets passed through the recursion.
        It tracks the entire recursive group of parent terms of the original base GO id (i.e. the starting point
        of the function call).

    Returns
    -------
    None
        Updates the parentSet set inplace so that it contains all the (recursive) parents for the baseGOid.
    """

    # If current term has no further parents the recursion will end and move back up the stack,
    # since there are no more parents left to iterate over (because looping through None does nothing).
    parents = GOdict.get(currentTerm).parents

    # For each parent of the current term under consideration
    for parent in parents:
        # Check if parent is present in GO dictionary
        # This is required because namespace filtering might lead to parent terms
        # that are no longer present as GOterm objects themselves.
        if parent in GOdict:
            # Add current term's parents to growing set
            parentSet.add(parent)
            # and recurse function for each parent
            propagateParents(parent, baseGOid, GOdict, parentSet)

        else:
            # Print a warning that a parent term was reported for the original base term,
            # yet the term is absent from the gene ontology file
            print('WARNING!\n' + parent, 'was defined as a parent for',
                  baseGOid, ', but was not found in the OBO file.\n')

    return None


def completeChildHierarchy(GOdict):
    """
    Generates the entire GO tree's child structure.

    By iterating over the parents of each GO object, each term's recursive_children attribute will be filled with
    a set of all recursively found child terms.

    NOTE: completeParentsHierarchy() must be run prior to this function!

    Parameters
    ----------
    GOdict : dict
        A dictionary of GO objects generated by importOBO().
        Keys are of the format `GO-0000001` and map to OBO objects.
        The recursive_parents attribute of the GO objects must be complete.

    Returns
    -------
    None
        Updates the provided GO dictionary inplace so that for each GO term object the "children" attribute points
        to the immediate children of the GO term and the "recursive_children" attribute traces back over the full
        GO hierarchy.
    """

    # For every GO term
    for GOid, GOobj in GOdict.items():
        # add this term as a child for all of its parents
        [GOdict[parent].recursive_children.add(GOid) for parent in GOobj.recursive_parents]

    return None


def assign_depth(node, GOdict, depth=0):
    """Recursive function to assign depth to node and all of its children.

    Starting from the root node of a directed acyclic graph (i.e. the gene ontology hierarchy), provided that all
    nodes have their child nodes as an attribute, the depth of each node will be set.
    Note: the minimal depth is always assigned. The depth of the root node is set to zero.

    Parameters
    ----------
    node : goTerm object
        The GO node that is the starting point for the recursive assignment, normally the root nodes of the
        three GO ontologies.
    GOdict : dict
        A dictionary of GO objects generated by importOBO().
        Keys are of the format `GO-0000001` and map to OBO objects.
    depth : int
        The current depth of the recursion, assigned to each node as its depth in the hierarchy.

    Returns
    -------
    None
        Modifies the goTerm objects' depth attribute inplace.


    """
    # assign depth value if node has no depth assigned yet
    if node.depth is None:
        node.depth = depth
    # if node has depth already, replace it only if the new depth is lower
    else:
        node.depth = min(node.depth, depth)
    # recurse through child nodes
    for i in node.children:
        assign_depth(GOdict.get(i), GOdict, depth+1)
