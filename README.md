# goscripts - GO enrichment analysis tools

A ready to use `python` script to perform GO enrichment tests by inputting a list of uniprot_kb accession numbers, an ontology (`.obo`) file and a gene association (`.gaf`) file.

The `goscripts` module contains further functionality to parse and manipulate `.obo` and `.gaf` files; e.g. retrieving parent/child terms by traversing the hierarchy, remapping a given set of GO terms to a specified depth, etc.

## Installation

### Using the main enrichment test script without installing additional packages

1) Download or clone this repository: `git clone git@github.com:pmoris/goscripts.git`
2) Run the `go_enrichment_script.py`.

If you wish to move the script to another location, be sure to also copy the `goscripts` directory with it.

### Installing all `goscripts` functionality as a python package

1) Optionally: create a new virtual or conda environment: `conda create -n goenv`
2) Install the goscripts package:
    * **Directly**:

            pip install git+https://github.com/pmoris/goscripts.git

    * **Manually**: Download or clone the repository and from within the main directory (where `setup.py` resides), let pip install the package:

            git clone git@github.com:pmoris/goscripts.git
            cd goscripts/
            pip install . # don't forget the dot

3) You can now use the stand-alone go enrichment testing script (which can be moved to any location) by invoking:

        python go_enrichment_script.py

4) All additional `goscripts` functionality can be used in your own python scripts via an import statement:

        import goscripts

## Usage

    usage: go_enrichment_script.py [-h] [-b BACKGROUND] -s SUBSET -o OBO -g GAF
                                [-O OUTPUTFILE]
                                [-n {all,biological_process,molecular_function,cellular_component}]
                                [-m MINGENES] [-l TESTING_LIMIT] [-p THRESHOLD]
                                [--mult-test MULT_TEST] [-v] [--no-propagation]
                                [--no-part-of]

    Script to perform GO enrichment analysis

    optional arguments:
    -h, --help            show this help message and exit
    -b BACKGROUND, --background BACKGROUND
                            File containing a list of Uniprot accession numbers
                            for the background set of genes. If omitted, the full
                            list of genes in the .gaf file will be used. (default:
                            set())
    -s SUBSET, --subset SUBSET
                            File containing a list of Uniprot accession numbers
                            for the subset of genes of interest (default: None)
    -o OBO, --obo OBO     .obo file containing gene ontology (default: None)
    -g GAF, --gaf GAF     .gaf file containing GO associations (default: None)
    -O OUTPUTFILE, --output OUTPUTFILE
                            Output file or path (default: enrichment_results.csv)
    -n {all,biological_process,molecular_function,cellular_component}, --namespace {all,biological_process,molecular_function,cellular_component}
                            Select the GO namespace to limit the enrichment test
                            to (default: all)
    -m MINGENES, --min MINGENES
                            Minimum number of genes before considering a GO
                            category (default: 3)
    -l TESTING_LIMIT, --limit-tests TESTING_LIMIT
                            P-value cut-off to use to stop GO tree propagation
                            during enchrichment tests (default: 0.05)
    -p THRESHOLD, --pval-thresh THRESHOLD
                            Significant p-value threshold to use for significance
                            testing (default: 0.1)
    --mult-test MULT_TEST
                            The type of multiple testing correction to use. Either
                            "fdr" or "bonferroni". (default: fdr)
    -v, --verbose         Verbose output. (default: False)
    --no-propagation      Disables propagation during testing. Use if only
                            strictly associated terms should be tested. (default:
                            True)
    --no-part-of          Ignore part_of relations between GO terms during
                            traversal. (default: False)

## Input files

* Ontology .obo files are described and available at the [Gene Ontology Consortium](http://www.geneontology.org/page/download-ontology).
* The gene association file format is described at the [Gene Ontology Consortium](http://www.geneontology.org/page/go-annotation-file-formats) and made available by EBI at the [GOA ftp site](https://www.ebi.ac.uk/GOA/downloads).
* The `background` and `subset` files should be plain text files containing a single Uniprot accession number per line.

    P00750
    A2BC19
    P12345
    A0A022YWF9

## Details

Performs one-sided hypergeometric tests, starting from the most specific (child) GO terms associated with the genes in the set of interest. If the p-value of the test does not fall below the specified significance level alpha, the test will be carried out for all of the term's parent terms, otherwise the process will terminate. This method attempts to limit the total number of tests that need to be carried out, since a term that is enriched will likely also have enriched parent terms. Furthermore, GO terms associated with a small number of genes are skipped. Next, the Benjamini-Hochberg FDR or Bonferroni multiple testing correction are applied to the test results. Finally, a `.csv` file containing all the GO terms that were evaluated and their p-values are returned. More information is available in the docstrings.

---

## Dependencies

    numpy
    pandas
    scipy.stats
    statsmodels.sandbox.stats.multicomp

---

Copyright (c) 2018 Pieter Moris
Adrem Data Lab - biomina - UAntwerpen
