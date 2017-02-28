#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@author: Pieter Moris
'''

import argparse
import pandas as pd

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
    parser.add_argument('-b', '--background', type=str, default='full annotation set', dest='background',
                        help='File containing a list of uniprot accession numbers for the background set of genes')
    parser.add_argument('-s', '--subset', type=str, required=True, dest='subset',
                        help='File containing a list of uniprot accession numbers for the subset of genes of interest')
    parser.add_argument('-o', '--obo', type=str, required=True, dest='obo',
                        help='.obo file containing gene ontology')
    parser.add_argument('-g', '--gaf', type=str, required=True, dest='gaf',
                        help='.gaf file containing GO annotations')
    # parser.add_argument('-O', '--output', type=argparse.FileType('w'),
    # default='enrichment-results.csv', dest='output',
    # help='Specifies the name or path of the output csv file')
    parser.add_argument('-O', '--output', type=str, default='enrichment_results.csv', dest='outputFile',
                        help='Output file or path')
    parser.add_argument('-m', '--min', type=int, default='3', dest='minGenes',
                        help='Minimum number of genes before considering a GO category')
    parser.add_argument('-T', '--testing-thresh', type=float, default='0.05', dest='testing_threshold',
                        help='P-value cut-off to use to stop GO tree propogation during testing')
    parser.add_argument('-t', '--thresh', type=float, default='0.1', dest='threshold',
                        help='P-value cut-off to use for significance testing')
    parser.add_argument('-s', '--suppress', action="store_false",
                        help='Suppress verbose output.')
    args = parser.parse_args()

    # Import the uniprot set of interest
    interest = genelist_importer.importSubset(args.subset)

    # If no background is provided, import the gene association file and retrieve
    # the full background set from it
    if args.background == 'full annotation set':
        print('No background gene set provided, retrieving all genes from the gene annotation file...')
        gafDict = gaf_parser.importFullGAF(args.gaf)
        background = set(gafDict.keys())
    # otherwise import the background set and use this to limit gene association import
    else:
        background = genelist_importer.importBackground(args.background)
        gafDict = gaf_parser.importGAF(args.gaf, background)

    # Generate a gene association file for the subset of interest too
    gafSubset = gaf_parser.createSubsetGafDict(interest, gafDict)

    # Check if the background/interest gene sets contain genes not present in the gene association file
    background = gaf_parser.removeObsoleteGenes(background, gafDict, 'background')
    interest = gaf_parser.removeObsoleteGenes(interest, gafSubset, 'interest')

    # Check if the set of interest contains uniprot ac's not present in the background set
    interest = genelist_importer.isValidSubset(interest, background)

    # import gene ontology file
    GOterms = obo_tools.importOBO(args.obo)

    # for i in ['GO:0060373', 'GO:0048245', 'GO:0044077', 'GO:0042925', 'GO:1902361', 'GO:1902361', 'GO:1902361', 'GO:0000001', 'GO:0000002']:
    #     print('id', i, 'parents', GOterms[i].parents)

    # Build full GO hierarchy
    print('Propagating through ontology to find all children and parents for each term...')
    obo_tools.buildGOtree(GOterms)

    # print('\n After propagating through parents')
    #
    # for i in ['GO:0060373', 'GO:0048245', 'GO:0044077', 'GO:0042925', 'GO:1902361', 'GO:1902361', 'GO:1902361', 'GO:0000001', 'GO:0000002']:
    #     print('id', i, 'parents', GOterms[i].parents)

    # Perform enrichment test
    pValues = enrichment_stats.enrichmentAnalysis(background, interest, GOterms, gafDict, gafSubset, minGenes=args.minGenes, threshold=args.testing_threshold)

    # print('\nTested GO-terms and their p-values:')
    # print(pValues)
    # print('\nOrdered p-values:')
    # print(sorted(pValues.values()))

    # Perform multiple testing correction
    correctedPvalues = enrichment_stats.multipleTestingCorrection(pValues, threshold = args.threshold)

    if not args.suppress:
        print(correctedPvalues)

    # Create dataframe with tested GO terms and results
    output = enrichment_stats.annotateOutput(correctedPvalues, GOterms)

    if not args.suppress:
        print(correctedPvalues)
        print('\nGO term, uncorrected and FDR-corrected p-values and description of GO term:')
        print(output)

    print('Saving output to', args.outputFile)
    output.to_csv(args.outputFile)