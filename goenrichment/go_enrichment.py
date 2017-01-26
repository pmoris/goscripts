#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@author: Pieter Moris
'''


import argparse
import os
import re

import pandas as pd
from scipy.stats import hypergeom

# import local modules


def importGAF(path, geneSet):
    """
    Imports a GAF file (gene association format) and generates a dictionary mapping the gene uniprot AC to the GO ID.
    Only imports genes which are present in the provided (background) gene set.

    Information on the GAF 2.1 format can be found at 
        http://geneontology.org/page/go-annotation-file-gaf-format-21 

    Parameters
    ----------
    path : str
        The path to the file.
    geneSet : set
        A set containing the uniprot AC's of all the genes under consideration (background).

    Returns
    -------
    dict of str mapping to set
        A dictionary mapping gene uniprot AC's to GO ID's.


    Possible improvements:
        Check for `is_obsolete` and `replaced_by`, although the replacement term should be in OBO file as an entry.

        Check for inclusion in provided geneset afterwards using:
            gafDict = {key: value for key, value in gafDict.items() if key in geneSet }

        To do: double dictionary? already filter based on gene list? what to do with NOT?
    """

    gafPath = os.path.abspath(path)

    if not geneSet:
        print('Empty gene set was provided during creation of GAF dictionary.')
        exit()

    with open(gafPath, 'r') as gafFile:
        gafDict = {}
        for line in gafFile:
            if line.startswith('UniProtKB'):
                splitLine = line.split('\t')                # split column-wise
                uniprotAC = splitLine[1]
                goTerm = splitLine[4]
                goQualifier = splitLine[3]
                if 'NOT' not in goQualifier:                # ignore annotations with "NOT"
                    if uniprotAC in geneSet:                # only keep genes present in background gene set
                        if uniprotAC not in gafDict:        # Create new key if AC does not already appear in dictionary
                            # make dictionary containing uniprot AC as key and
                            gafDict[uniprotAC] = {goTerm}
                        else:                               # set of GO's as value
                            gafDict[uniprotAC].add(goTerm)

    print('Retrieved', len(gafDict),
          'annotated (background filtered) uniprot AC\'s from', gafPath + '\n')

    if len(gafDict) != len(geneSet):
        print('WARNING!\nNot every uniprot AC that was provided in the background set was found in the GAF file:')
        print([AC for AC in geneSet if AC not in gafDict])
        print('WARNING!\n')

    return gafDict


class goTerm:
    """
    GO term object.

    Stores the ID, name and domain of the GO term and contains dictionaries for child and parent nodes.

    Attributes
    ----------
    id : str
        The identifier of the GO term.
    altid : str
        Optional tag for an alternative id.
    name : str
        The GO term name.
    namespace : str
        The domain of the GO term (Cellular Component, Molecular Function or Biological Process).
    parents : set of str
        The parent terms of the GO term, as indicated by the `is_a` relationship.
    childs : set of str
        The child terms of the GO term, derived from other GO terms after a complete OBO file is processed initially.

    not necessary... Methods
    -------
    returnID
        Returns the ID of the GO term.
    gamma(n=1.0)
        Change the photo's gamma exposure.  
    # https://stackoverflow.com/questions/1336791/dictionary-vs-object-which-is-more-efficient-and-why

    '''
    https://stackoverflow.com/questions/3489071/in-python-when-to-use-a-dictionary-list-or-set
    When you want to store some values which you'll be iterating over, 
    Python's list constructs are slightly faster. 
    However, if you'll be storing (unique) values in order to check for their existence, 
    then sets are significantly faster.
    '''
    """

    goCount = 0

    __slots__ = ('id', 'name', 'altid', 'namespace', 'childs', 'parents')

    def __init__(self, GOid):
        self.id = GOid
        self.altid = []
        self.name = ''
        self.namespace = ''
        self.childs = set()
        self.parents = set()

        goTerm.goCount += 1

    def returnID(self):
        return self.id


def importOBO(path):
    """
    Imports an OBO file and generates a dictionary containing an OBO object for each GO term.

    Parameters
    ----------
    path : str
        The path to the file.

    Returns
    -------
    dict of OBO objects
        Keys are of the format `GO-0000001` and map to OBO objects.

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
                # and hierarchy from subsequent lines
                if line.startswith('id'):
                    # Store ID for lookup of other attributes in next lines
                    GOid = line.split(': ')[1].rstrip()

                    if not GOid in GOdict:               # check if id is already stored as a key in dictionary
                        # if not, create new GOid object as the value for this
                        # key
                        GOdict[GOid] = goTerm(GOid)

                elif line.startswith('name'):
                    GOdict[GOid].name = line.split(': ')[1].rstrip()
                elif line.startswith('namespace'):
                    GOdict[GOid].namespace = line.split(': ')[1].rstrip()
                elif line.startswith('alt_id'):
                    GOdict[GOid].altid.append(line.split(': ')[1].rstrip())
                elif line.startswith('is_a'):
                    GOdict[GOid].parents.add(line.split()[1].rstrip())
    return GOdict


def importBackground(path):
    """
    Imports the background set of genes (uniprot AC).

    Parameters
    ----------
    path : str
        The path to the file.

    Returns
    -------
    set of str
        A set of background uniprot AC's.

    Notes: Gene lists should not contain a header. One gene per line.
    Possible improvement: check for file structure and allow headers, comma separated lists, etc.
    """
    with open(path, 'r') as inGenes:
        backgroundSet = set([line.rstrip() for line in inGenes][1:])
    print('Retrieved', len(backgroundSet),
          'background uniprot AC\'s from', path)
    return backgroundSet


def importSubset(path):
    """
    Imports the gene subset of interest (uniprot AC).

    Parameters
    ----------
    path : str
        The path to the file.

    Returns
    -------
    set of str
        A subset of uniprot ACs of interest.

    Notes: Gene lists should not contain a header. One gene per line.
    """
    with open(path, 'r') as inGenes:
        geneSubset = set([line.rstrip() for line in inGenes][1:])
    print('Retrieved', len(geneSubset), 'subset uniprot AC\'s from', path)
    return geneSubset

def isValidSubset(subset, background):
    if subset.issubset(background):
        return True
    else:
        print('WARNING! Subset contained genes not present in background list. Terminating program now.')
        exit()


def createSubsetGafDict(subset, gafDict):
    """
    Generates a dictionary mapping the subset's gene uniprot AC's to the GO ID's,
    based on the provided gene subset and the gaf dictionary.

    Parameters
    ----------
    subset : set of str
        A subset of uniprot ACs of interest.

    Returns
    -------
    dict of str mapping to set
        A dictionary mapping the subset's gene uniprot AC's to GO ID's.
    """

    GafSubset = {}

    for gene, GOids in gafDict.items():
        if gene in subset:
            GafSubset[gene] = GOids

    if len(subset) != len(GafSubset):
        print('WARNING!\nNot every uniprot AC that was provided in the gene subset was found in the GAF file:')
        print([AC for AC in subset if AC not in GafSubset])
        print('WARNING!\n')

    return GafSubset


def buildGOtree(GOdict):
    """
    Generates the entire GO tree's parent structure by walking through the parent hierarchy of each GO entry.

    Parameters
    ----------
    GOdict : dict
        A dictionary of GO objects generated by importOBO().
        Keys are of the format `GO-0000001` and map to OBO objects.

    Returns
    -------
    None
        The input GO dictionary will be updated.
        Parent attributes now trace back over the full tree hierarchy.
    """

    # Process each GO term in the GO dictionary
    for GOid, GOobj in GOdict.items():
        # Define new set to store higher order parents
        parentSet = set()
        # Call helper function to propagate through parents
        propagateParents(GOid, GOid, GOdict, parentSet)
        # Update GO term's parents attribute to include all higher order
        # parents
        GOobj.parents.update(parentSet)

    # After all parents have been found, for each ID, add it as a child for
    # all its parents
    completeChildHierarchy(GOdict)

    return None


def propagateParents(currentTerm, baseGOid, GOdict, parentSet):
    """
    Propagates through the parent hierarchy of a provided GO term.

    Parameters
    ----------
    GOdict : dict
        A dictionary of GO objects generated by importOBO().
        Keys are of the format `GO-0000001` and map to OBO objects.

    Returns
    -------
    dict of OBO objects
        Updated GO dict where parent attributes trace back over the full tree hierarchy.
        Keys are of the format `GO-0000001` and map to OBO objects.
    """

    # If current term is not present in GO dictionary, print warning and end
    # recursion
    if currentTerm not in GOdict:
        print('WARNING!\n' + currentTerm, 'was defined as a parent for',
              baseGOid, ', but was not found in the OBO file.')
        parentSet.pop(currentTerm)      # remove missing value
        return

    # If current term has no further parents the recursion will end and move
    # back up the stack, since there are no parents to iterate over
    parents = GOdict.get(currentTerm).parents
    for parent in parents:

        # # Check if parent is present in GO dictionary
        # if parent not in GOdict:
        #     print('WARNING!\n' + parent, 'was defined as a parent for',
        #           baseGOid, ', but was not found in the OBO file.')

        # Add current term's parents to growing set
        parentSet.add(parent)
        # NOTE: better to do this at the start of the function, by adding current term
        #       otherwise a term that is not present in the OBO dict will be added
        #       since the check happens later, i.e. in the next function call

        # and recurse function for each parent
        propagateParents(parent, baseGOid, GOdict, parentSet)

    return None


def completeChildHierarchy(GOdict):
    """
    Generates the entire GO tree's child structure by iterating over the parents
    of each GO object.

    NOTE: completeParentsHierarchy() must be run prior to this function.

    Parameters
    ----------
    GOdict : dict
        A dictionary of GO objects generated by importOBO().
        Keys are of the format `GO-0000001` and map to OBO objects.

    Returns
    -------
    dict of OBO objects
        Updated GO dict where child attributes trace back over the full tree hierarchy.
        Keys are of the format `GO-0000001` and map to OBO objects.
    """
    for GOid, GOobj in GOdict.items():
        [GOdict[parent].childs.add(GOid) for parent in GOobj.parents]

    return None

# Different propagation method
# http://blog.nextgenetics.net/?e=6
# if level is desired, might be better to propagate for each GOid during
# calculation


###############################################################################
################## Hypergeometric testing ###############
###############################################################################

def enrichmentOneSided(GOid, background, subset, GOdict, gafDict, gafSubset):
    """
    Performs a one-sided hypergeometric test for a given GO term.

    Parameters
    ----------
    GOid : str
        A GO identifier (key to the GO dictionary).
    background : set of str
        A set of background uniprot AC's.
    subset : set of str
        A subset of uniprot AC's of interest.
    GOdict : dict
        A dictionary of GO objects generated by importOBO().
        Keys are of the format `GO-0000001` and map to OBO objects.
    gafDict : dict
        A dictionary mapping the background's gene uniprot AC's to GO ID's.        
    gafDict : dict
        A dictionary mapping the subset's gene uniprot AC's to GO ID's.

    Returns
    -------
    float
        The p-value of the one-sided hypergeometric test.
    """

    backgroundTotal = len(background)
    subsetTotal = len(subset)

    validTerms = set([GOid])
    validTerms.update(GOdict['GOid'].childs)

    backgroundGO = countGOassociations(validTerms, gafDict)
    subsetGO = countGOassociations(validTerms, gafDict)

    # k or more successes (= GO associations = subsetGO) in N draws (= subsetTotal)
    # from a population of size M (backgroundTotal) containing n successes (backgroundGO)
    # k or more is the sum of the probability mass functions of k up to N successes
    # since cdf gives the cumulative probability up and including input (less or equal to k successes),
    # and we want P(k or more), we need to calculate 1 - P(less than k) =  1 - P(k-1 or less)
    # .sf is the survival function (1-cdf).
    pVal = hypergeom.sf(subsetGO - 1, backgroundTotal,
                        backgroundGO, subsetTotal)

    return pVal


def countGOassociations(validTerms, gaf):
    """
    """

    GOcounter = 0

    # For each gene:GO id set pair in the GAF dictionary
    for gene, GOids in gaf.items():
        # Increment the GO counter if the valid terms set shares a member
        # with the GO id set of the current gene
        if not validTerms.isdisjoint(GOids):
            GOcounter += 1

    return GOcounter


def enrichmentAnalysis(background, subset, GOdict, gafDict, gafSubset):

    # generate a list of all base GO id's to test
    # i.e. those of all genes in the subset of interest

    GOidsToTest = []


    for gene, goids in gafSubset.items():

    # 


# Run as stand alone script
# https://stackoverflow.com/questions/419163/what-does-if-name-main-do
if __name__ == '__main__':

    # Check provided arguments
    parser = argparse.ArgumentParser(
        description='Script to perform GO enrichment analysis')
    # https://stackoverflow.com/questions/18862836/how-to-open-file-using-argparse
    parser.add_argument('-b', '--background', type=str, required=True, dest='background',
                        help='File containing a list of uniprot accession numbers for the background set of genes')
    parser.add_argument('-s', '--subset', type=str, required=True, dest='subset',
                        help='File containing a list of uniprot accession numbers for the subset of genes of interest')
    parser.add_argument('-o', '--obo', type=str, required=True, dest='obo',
                        help='.obo file containing gene ontology')
    parser.add_argument('-g', '--gaf', type=str, required=True, dest='gaf',
                        help='.gaf file containing GO annotations')
    parser.add_argument('-O', '--output', type=argparse.FileType('w'),
                        default='enrichment-results.csv', dest='output',
                        help='Specifies the name or path of the output csv file')
    parser.add_argument('-m', '--min', type=int, default='3', dest='minGenes',
                        help='Minimum number of genes before considering a GO category')
    parser.add_argument('-t', '--thresh', type=float, default='0.1', dest='threshold',
                        help='P-value cut-off to use for significance testing')
    args = parser.parse_args()

    background = importBackground(args.background)
    interest = importSubset(args.subset)
    isValidSubset(interest, background)

    gafDict = importGAF(args.gaf, background)
    gafSubset = createSubsetGafDict(interest, gafDict)

    GOterms = importOBO(args.obo)

    print('background\n', background)
    print('interest\n', interest)

    for i in ['GO:0060373', 'GO:0048245', 'GO:0044077', 'GO:0042925', 'GO:1902361', 'GO:1902361', 'GO:1902361', 'GO:0000001', 'GO:0000002']:
        print('id', i, 'parents', GOterms[i].parents)

    buildGOtree(GOterms)

    for i in ['GO:0060373', 'GO:0048245', 'GO:0044077', 'GO:0042925', 'GO:1902361', 'GO:1902361', 'GO:1902361', 'GO:0000001', 'GO:0000002']:
        print('id', i, 'parents', GOterms[i].parents)


'''
# need to store iterative parents in list/set and append only after
# propagation is terminated. otherwise set is changed during run

    # Check if parent is present in GO dictionary
    if parent not in GOdict:
        print('WARNING!\n' + parent, 'was defined as a parent for',
              baseGOid, ', but was not found in the OBO file.')

    else:
        storer.add(parent)
        # Add parent to base GO term's parent attribute set
        GOdict.get(baseGOid).parents.add(parent)

        # Propagate further through parent hierarchy
        # by calling function on each of the current parent's parents
        parentParents = GOdict.get(parent).parents
        for parent in parentParents:
            propagateTree(parent, baseGOid, GOdict, storer)

    # # Stop recursion if previous call came from a GO term without parents
    # if not parents:
    #     return False

    # else:
    #     for parent in parents:
    #         # Add current parent term to base GO term's parent set
    #         GOdict.get(baseGOid).parents.add()
    #         # Propagate through parent terms of current parents
    #         propagateTree(GOdict.get(parent).parents, baseGOid, GOdict, storer)
    #         storer.add(parent)

    # # Check if current (parent) term still has parents of its own to process
    # if GOdict.get(GOid).parents:
    #     gdg

    # if GOid not in GOdict[parent].child:
    #     GOdict[parent].child.add(GOid)

    # if GOdict[parent]:
    #     for parentsParent in GOdict[parent].parents:
    #         propagateTree(parentsParent, GOid, GOdict)
    # else:
    #     return False

'''

'''



































































# Import gene list
inputPath = os.path.abspath(args.genes)
with open(inputPath, 'r') as inGenes:
    geneSet = set([line.rstrip() for line in inGenes][1:])

def importGAF(path, geneSet):
    """
    Imports a GAF file and generates a dictionary mapping the gene ID to the GO ID.

    Information on the GAF 2.1 format can be found at 
        http://geneontology.org/page/go-annotation-file-gaf-format-21 

    Parameters
    ----------
    path : str
        The path to the file.
    geneSet : set
        A set containing the ID's of the genes under consideration.

    Returns
    -------
    dict of str
        A dictionary mapping gene ID's to GO ID's.

    Possible improvements:
        Check for `is_obsolete` and `replaced_by`, although the replacement term should be in OBO file as an entry.

        Check for inclusion in provided geneset afterwards using:
            gafDict = {key: value for key, value in gafDict.items() if key in geneSet }

        To do: double dictionary? already filter based on gene list? what to do with NOT?
    """
    gafPath = os.path.abspath(path)

    if not geneSet:
        print('Empty gene set provided during creation of GAF dictionary.')
        exit()

    with open(gafPath, 'r') as gafFile:
        gafDict = {}
        for line in gafFile:
            # ignore comments in gaf file
            if not line.startswith('!'):
                splitLine = line.split('\t')                # split column-wise
                uniprotAC = splitLine[1]
                goTerm = splitLine[4]
                goQualifier = splitLine[3]
                if 'NOT' not in goQualifier:                # ignore annotations with "NOT"
                    if uniprotAC in geneSet:                # only keep genes present in input gene set
                        if uniprotAC not in gafDict:        # Create new key if AC does not already appear in dictionary
                            # make dictionary containing uniprot AC as key and
                            # set of GO's
                            gafDict[uniprotAC] = {goTerm}
                        else:
                            gafDict[uniprotAC].add(goTerm)
    return gafDict

# https://stackoverflow.com/questions/1336791/dictionary-vs-object-which-is-more-efficient-and-why

"""
"""https://stackoverflow.com/questions/3489071/in-python-when-to-use-a-dictionary-list-or-set
When you want to store some values which you'll be iterating over, 
Python's list constructs are slightly faster. 
However, if you'll be storing (unique) values in order to check for their existence, 
then sets are significantly faster.
"""


class goTerm:
    """
    GO term object.

    Stores the ID, name and domain of the GO term and contains dictionaries for child and parent nodes.

    Attributes
    ----------
    id : str
        The identifier of the GO term.
    altid : str
        Optional tag for an alternative id.
    name : str
        The GO term name.
    namespace : str
        The domain of the GO term (Cellular Component, Molecular Function or Biological Process).
    parents : set of str
        The parent terms of the GO term, as indicated by the `is_a` relationship.
    childs : set of str
        The child terms of the GO term, derived from other GO terms after a complete OBO file is processed initially.

not necessary... Methods
    -------
    returnID
        Returns the ID of the GO term.
    gamma(n=1.0)
        Change the photo's gamma exposure.  
    """

    goCount = 0

    __slots__ = ('id', 'name', 'altid', 'namespace', 'childs', 'parents')

    def __init__(self, GOid):
        self.id = GOid
        self.altid = ''
        self.name = ''
        self.namespace = ''
        self.childs = set()
        self.parents = set()

        goTerm.goCount += 1

    def returnID(self):
        return self.id

def importOBO(path):
    """
    Imports an OBO file and generates a dictionary containing an OBO object for each GO term.

    Information on the OBO 1.2 format can be found at 
        https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_2.html

    Parameters
    ----------
    path : str
        The path to the file.

    Returns
    -------
    dict of OBO objects
        Keys of the format `GO-0000001` map to OBO objects.

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
                # and hierarchy from subsequent lines
                if line.startswith('id'):
                    # Store ID for lookup of other attributes in next lines
                    GOid = line.split(': ')[1].rstrip()

                    if not GOid in GOdict:               # check if id is already stored as a key in dictionary
                        # if not, create new GOid object as the value for this
                        # key
                        GOdict[GOid] = goTerm(GOid)

                elif line.startswith('name'):
                    GOdict[GOid].name = line.split(': ')[1].rstrip()
                elif line.startswith('namespace'):
                    GOdict[GOid].namespace = line.split(': ')[1].rstrip()
                elif line.startswith('alt_id'):
                    GOdict[GOid].altid = line.split(': ')[1].rstrip()
                elif line.startswith('is_a'):
                    GOdict[GOid].parents = line.split()[1].rstrip()
    return GOdict

'''
