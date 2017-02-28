#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@author: Pieter Moris
'''

import os


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

    backgroundPath = os.path.abspath(path)
    with open(backgroundPath, 'r') as inGenes:
        backgroundSet = set([line.rstrip() for line in inGenes])
    print('Retrieved', len(backgroundSet),
          'background uniprot AC\'s from', backgroundPath)
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

    subsetPath = os.path.abspath(path)
    with open(subsetPath, 'r') as inGenes:
        geneSubset = set([line.rstrip() for line in inGenes])
    print('Retrieved', len(geneSubset), 'subset uniprot AC\'s from', subsetPath)
    return geneSubset


def isValidSubset(subset, background):
    """
    Checks if the gene subset of interest contains genes not present in the background set.
    If there are additional genes they are removed.

    Parameters
    ----------
    subset : set of str
        A subset of uniprot ACs of interest.
    background : set of str
        A set of uniprot ACs to be used as the background.

    Returns
    -------
    set of str
        A cleaned subset of uniprot ACs of interest.
    """

    if subset.issubset(background):
        return subset
    else:
        print('WARNING! Subset contained genes not present in background list...')
        print([AC for AC in subset if AC not in background])
        print('Removing these genes from the set of interest...\n')
        print(subset)
        return subset.difference([AC for AC in subset if AC not in background])

