# -*- coding: utf-8 -*-
'''
Created on Wed January 11 2017
@author: Pieter Moris
'''

import pandas as pd
import argparse
import os

# Check provided arguments
parser = argparse.ArgumentParser(
    description='Script to perform GO enrichment analysis')

parser.add_argument('-i', '--inputFile', type=str, required=True, dest='genes',
                    help='File containing a list of uniprot accession numbers for the genes of interest')
parser.add_argument('-b', '--obo', type=str, required=True, dest='obo',
                    help='.obo file containing gene ontology')
parser.add_argument('-g', '--gaf', type=str, required=True, dest='gaf',
                    help='.gaf file containing GO annotations')
parser.add_argument('-O', '--output', type=argparse.FileType('w'),
                    default='enrichment-results.csv', dest='output',
                    help='Specifies the name or path of the output csv file')
parser.add_argument('-m', '--min', type=int, default='3', dest='minGenes',
                    help='Minimum number of genes before considering a GO category')
parser.add_argument('-t', '--thresh', type=int, default='0.1', dest='threshold',
                    help='P-value cut-off to use for significance testing')

args = parser.parse_args()

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

'''
https://stackoverflow.com/questions/3489071/in-python-when-to-use-a-dictionary-list-or-set
When you want to store some values which you'll be iterating over, 
Python's list constructs are slightly faster. 
However, if you'll be storing (unique) values in order to check for their existence, 
then sets are significantly faster.
'''


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

    def __init__(self, goID):
        self.id = goID
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
                    goID = line.split(': ')[1].rstrip()

                    if not goID in GOdict:               # check if id is already stored as a key in dictionary
                        # if not, create new goID object as the value for this
                        # key
                        GOdict[goID] = goTerm(goID)

                elif line.startswith('name'):
                    GOdict[goID].name = line.split(': ')[1].rstrip()
                elif line.startswith('namespace'):
                    GOdict[goID].namespace = line.split(': ')[1].rstrip()
                elif line.startswith('alt_id'):
                    GOdict[goID].altid = line.split(': ')[1].rstrip()
                elif line.startswith('is_a'):
                    GOdict[goID].parents = line.split()[1].rstrip()
    return GOdict

