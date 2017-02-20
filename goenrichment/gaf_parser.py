#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@author: Pieter Moris
'''

import os


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
        A dictionary mapping gene uniprot AC's (str) to a set GO ID's.


    Possible improvements:
        Check for `is_obsolete` and `replaced_by`, although the replacement term should be in OBO file as an entry.

        Check for inclusion in provided gene set afterwards using:
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

    if len(gafDict) != len(geneSet):
        print('WARNING!\nNot every uniprot AC that was provided in the background set was found in the GAF file:')
        print([AC for AC in geneSet if AC not in gafDict])
        print('WARNING!\n')

    print('Retrieved', len(gafDict),
          'annotated (background filtered) uniprot AC\'s from', gafPath + '\n')

    return gafDict

def importFullGAF(path):
    """
    Imports a GAF file (gene association format) and generates a dictionary mapping the gene uniprot AC to the GO ID.
    Imports the full set of genes in the file.

    Information on the GAF 2.1 format can be found at
        http://geneontology.org/page/go-annotation-file-gaf-format-21

    Parameters
    ----------
    path : str
        The path to the file.

    Returns
    -------
    dict
        A dictionary mapping gene uniprot AC's (str) to a set of GO ID's.
   """

    gafPath = os.path.abspath(path)

    with open(gafPath, 'r') as gafFile:
        gafDict = {}
        for line in gafFile:
            if line.startswith('UniProtKB'):
                splitLine = line.split('\t')                # split column-wise
                uniprotAC = splitLine[1]
                goTerm = splitLine[4]
                goQualifier = splitLine[3]
                if 'NOT' not in goQualifier:                # ignore annotations with "NOT"
                    if uniprotAC not in gafDict:        # Create new key if AC does not already appear in dictionary
                        # make dictionary containing uniprot AC as key and
                        gafDict[uniprotAC] = {goTerm}
                    else:                               # set of GO's as value
                        gafDict[uniprotAC].add(goTerm)

    print('Retrieved', len(gafDict),
          'annotated uniprot AC\'s from', gafPath, 'after filtering on the background set.' + '\n')

    return gafDict


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
