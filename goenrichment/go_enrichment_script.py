#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@author: Pieter Moris
'''

import argparse

import enrichment_stats
import gaf_parser
import genelist_importer
import obo_tools


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

    background = genelist_importer.importBackground(args.background)
    interest = genelist_importer.importSubset(args.subset)
    if not genelist_importer.isValidSubset(interest, background):
        print("WARNING! Subset contains genes not present in background set!")
        print(interest - background)
        print('Terminating script.')
        exit()

    gafDict = gaf_parser.importGAF(args.gaf, background)
    gafSubset = gaf_parser.createSubsetGafDict(interest, gafDict)

    GOterms = obo_tools.importOBO(args.obo)

    print('background\n', background)
    print('interest\n', interest)

    for i in ['GO:0060373', 'GO:0048245', 'GO:0044077', 'GO:0042925', 'GO:1902361', 'GO:1902361', 'GO:1902361', 'GO:0000001', 'GO:0000002']:
        print('id', i, 'parents', GOterms[i].parents)

    obo_tools.buildGOtree(GOterms)

    for i in ['GO:0060373', 'GO:0048245', 'GO:0044077', 'GO:0042925', 'GO:1902361', 'GO:1902361', 'GO:1902361', 'GO:0000001', 'GO:0000002']:
        print('id', i, 'parents', GOterms[i].parents)

    pValues = enrichment_stats.enrichmentAnalysis(background, interest, GOterms, gafDict, gafSubset,
                       minGenes=3, threshold=0.05)

    print(pValues)