#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@author: Pieter Moris

Main script to perform a gene enrichment analysis.

# Method

Performs one-sided hypergeometric tests, starting from the most specific (child) GO terms associated with the
genes in the set of interest. If the p-value of the test does not fall below the specified significance level alpha,
the test wll be carried out for all of the term's parent terms, otherwise the process will terminate. This method
attempts to limit the total number of tests that need to be carried out, since a term that is enriched will likely
also have enriched parent terms. Furthermore, GO terms associated with a small number of genes are skipped.
Next, the Benjamini-Hochberg FDR or Bonferroni multiple testing correction are applied to the test results.
Finally, a csv file containing all the GO terms that were evaluated and their p-values are returned.

# Requires the following input files:
    - gene association file in .gaf format
    - gene ontology file in .obo format
    - text file with genes of interest separated by new lines
        expects uniprot AC's
    - optional text file with background set
        if not provided, the entire gene annotation found in the association file will be used

# The following options can be specified:
    - The output file in which results will be stored
    - The minimum number of genes required to be associated with a term before it will be tested.
    - The FDR cut-off to utilise.
    - The p-value threshold to limit the recursive testing of the GO tree.
    - The namespace to test: all, biological_process, molecular_function or cellular_component
        Uses the naming convention found in the .obo file.

# Dependencies:
    numpy, pandas, scipy.stats, statsmodels.sandbox.stats.multicomp

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
        description='Script to perform GO enrichment analysis', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    parser.add_argument('-n', '--namespace', type=str, default='all',
                        help='Select the GO namespace to limit the enrichment test to')
    parser.add_argument('-m', '--min', type=int, default='3', dest='minGenes',
                        help='Minimum number of genes before considering a GO category')
    parser.add_argument('-T', '--testing-thresh', type=float, default='0.05', dest='testing_threshold',
                        help='P-value cut-off to use to stop GO tree propagation during testing')
    parser.add_argument('-t', '--thresh', type=float, default='0.1', dest='threshold',
                        help='P-value cut-off to use for significance testing')
    parser.add_argument('-c', '--condense', action="store_true",
                        help='Suppress verbose output.')
    args = parser.parse_args()

    # Import the gene set of interest (uniprot AC format)
    interest = genelist_importer.importGeneList(args.subset)

    # If no background is provided, import the gene association file and
    # retrieve the full background set from it
    if args.background == 'full annotation set':
        print('No background gene set provided, retrieving all genes from the gene annotation file...\n')
        gafDict = gaf_parser.importGAF(args.gaf, args.background)
        background = set(gafDict.keys())
    # otherwise, import the background set first and use this to limit gene association import
    else:
        background = genelist_importer.importGeneList(args.background)
        gafDict = gaf_parser.importGAF(args.gaf, background)

    # Check if the set of interest contains uniprot AC's not present in the background set and
    # report and remove them
    interest = genelist_importer.isValidSubset(interest, background)

    # Generate a gene association file for the (pruned) subset of interest too.
    gafSubset = gaf_parser.createSubsetGafDict(interest, gafDict)

    # Check if the background/interest gene sets contain genes not present in the gene association file
    # and report that they won't be considered during the analysis
    background = genelist_importer.reportMissingGenes(background, gafDict, 'background')
    interest = genelist_importer.reportMissingGenes(interest, gafSubset, 'interest')

    # Import gene ontology file
    GOterms = obo_tools.importOBO(args.obo)

    print('lenGOterms',len(GOterms))
    print('lengafdict',len(gafDict))
    print('lensubsetgafdict',len(gafSubset))


    # Reduce gene ontology file to selected namespace
    if not args.namespace == 'all':
        GOterms = obo_tools.filterOnNamespace(GOterms, args.namespace)
        # and reduce the gene association files to these GO terms
        gafDict = gaf_parser.cleanGafTerms(gafDict, GOterms)
        gafSubset = gaf_parser.cleanGafTerms(gafSubset, GOterms)
        print(len(gafDict), 'out of', len(background), 'genes are annotated for', args.namespace)
        print(len(gafSubset), 'out of', len(interest), 'genes are annotated for', args.namespace)

    print('lenFilteredGOterms', len(GOterms))
    print('lenFilteredgafDict', len(gafDict))
    print('lenbackground',len(background))
    print('lensubsetgafdict',len(gafSubset))
    print('leninterest',len(interest))

    background = genelist_importer.reportMissingGenes(background, gafDict, 'background')
    interest = genelist_importer.reportMissingGenes(interest, gafSubset, 'interest')

    print('lenFilteredGOterms', len(GOterms))
    print('lenFilteredgafDict', len(gafDict))
    print('lenbackground',len(background))
    print('lensubsetgafdict',len(gafSubset))
    print('leninterest',len(interest))

    # TODO: remove interest/background set and get lengths from gaf/gafsubset dict?

    # for i in ['GO:0060373', 'GO:0048245', 'GO:0044077', 'GO:0042925', 'GO:1902361', 'GO:1902361', 'GO:1902361', 'GO:0000001', 'GO:0000002']:
    #     print('id', i, 'parents', GOterms[i].parents)

    # Build full GO hierarchy
    print('Propagating through ontology to find all children and parents for each term...\n')
    obo_tools.buildGOtree(GOterms)

    # print('\n After propagating through parents')
    #
    # for i in ['GO:0060373', 'GO:0048245', 'GO:0044077', 'GO:0042925', 'GO:1902361', 'GO:1902361', 'GO:1902361', 'GO:0000001', 'GO:0000002']:
    #     print('id', i, 'parents', GOterms[i].parents)

    # Perform enrichment test
    print('Finished completing ontology...proceeding with enrichment tests...\n')
    enrichmentResults = enrichment_stats.enrichmentAnalysis(background, interest, GOterms, gafDict, gafSubset, minGenes=args.minGenes, threshold=args.testing_threshold)

    # Update results with multiple testing correction
    enrichment_stats.multipleTestingCorrection(enrichmentResults, testType='fdr', threshold=args.threshold)

    # Create dataframe with tested GO terms and results
    output = enrichment_stats.annotateOutput(enrichmentResults, GOterms, background, interest)

    # Print intermediate output
    if not args.condense:
        print('\nGO term, uncorrected and FDR-corrected p-values and description of GO term:\m')
        print(output)

    # Save results
    print('Saving output to', args.outputFile)
    output.to_csv(args.outputFile)
