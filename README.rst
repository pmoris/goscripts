goscripts - Python script and package for Gene Ontology enrichment analysis
===========================================================================

The full documentation for this project is available at: https://pmoris.github.io/goscripts/

--------------

A ready to use ``python`` script to perform GO enrichment tests by
inputting a list of uniprot\_kb accession numbers, an ontology
(``.obo``) file and a gene association (``.gaf``) file. No install required!

The ``goscripts`` package provides further functionality to parse and
manipulate ``.obo`` and ``.gaf`` files; e.g.

-  Parsing ``.obo`` `gene ontology files <http://www.geneontology.org/page/download-ontology>`_
   in order to retrieve child/parent terms: ``goscripts.obo_tools``
-  Remapping gene ontology terms to a specified depth: ``goscripts.obo_tools``
-  Parse ``.gaf`` `gene association files <http://www.geneontology.org/page/go-annotation-file-formats>`_:
   ``goscripts.gaf_parser``
-  Performing an enrichment test using various multiple testing correction
   procedures (provided by `statsmodels <http://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html>`_): ``goscripts.enrichment_stats``

How to get started
------------------

Without installation
~~~~~~~~~~~~~~~~~~~~

1) Download or clone the repository:
   ``git clone git@github.com:pmoris/goscripts.git``
2) Run the script:
   ::

       python go_enrichment_script.py

The script requires functionality stored inside the ``goscripts`` directory
and expects to find this directory. Consequently, if you wish to move the
script to a different location, be sure to also copy this directory with it. Moreover, all :ref:`dependencies <dep-label>` should be installed.

Installing the ``goscripts`` package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1) Optionally: create a new virtual or conda environment:
   ``conda create -n goenv``
2) Install the goscripts package:

   -  **Directly**:

      ::

          pip install git+https://github.com/pmoris/goscripts.git

   -  **Manually**: Download or clone the repository and from within the
      main directory (where ``setup.py`` resides), let pip install the
      package:

      ::

          git clone git@github.com:pmoris/goscripts.git
          cd goscripts/
          pip install . # don't forget the dot

3) You can now use the ``go_enrichment_script.py`` from any location:

   ::

       python go_enrichment_script.py

4) All additional ``goscripts`` functionality can be used in your own
   python scripts via an import statement:

   ::

       import goscripts
       # or
       from goscripts import gaf_parser

Downloading required ontologies and annotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Example for human data:

::

    # Ontology file
    wget http://purl.obolibrary.org/obo/go.obo destination/directory
    # Annotation file
    wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz destination/directory
    gunzip destination/directory/goa_human.gaf.gz


Using the GO enrichment test script
-----------------------------------

::

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
                            "fdr_bh" (default), "bonferroni" or any other method
                            offered by
                            statsmodels.stats.multitest.multipletests(). (default:
                            fdr_bh)
    -v, --verbose         Verbose output. (default: False)
    --no-propagation      Disables propagation during testing. Use if only
                            strictly associated terms should be tested. (default:
                            True)
    --no-part-of          Ignore part_of relations between GO terms during
                            traversal. (default: False)

See the statsmodels documentation for an overview of all available
multiple testing correction procedures:
http://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html.

Input files
~~~~~~~~~~~

-  Ontology .obo files are described and available at the `Gene Ontology
   Consortium <http://www.geneontology.org/page/download-ontology>`__.
-  The gene association file format is described at the `Gene Ontology
   Consortium <http://www.geneontology.org/page/go-annotation-file-formats>`__
   and made available by EBI at the `GOA ftp
   site <https://www.ebi.ac.uk/GOA/downloads>`__.
-  The ``background`` and ``subset`` files should be plain text files
   containing a single Uniprot accession number per line.

   P00750 A2BC19 P12345 A0A022YWF9

Details
~~~~~~~

Performs one-sided hypergeometric tests, starting from the most specific
(child) GO terms associated with the genes in the set of interest. If
the p-value of the test does not fall below the specified significance
level alpha, the test will be carried out for all of the term's parent
terms, otherwise the process will terminate. This method attempts to
limit the total number of tests that need to be carried out, since a
term that is enriched will likely also have enriched parent terms.
Furthermore, GO terms associated with a small number of genes are
skipped. Next, the Benjamini-Hochberg FDR or Bonferroni multiple testing
correction are applied to the test results. Finally, a ``.csv`` file
containing all the GO terms that were evaluated and their p-values are
returned. More information is available in the docstrings.

--------------

.. _dep-label:

Dependencies
------------

::

    numpy
    pandas
    scipy.stats
    statsmodels.stats.multitest

--------------

Similar projects
----------------

-   https://github.com/tanghaibao/goatools
-   https://github.com/jdrudolph/goenrich
-   https://www.psb.ugent.be/cbd/papers/BiNGO/Home.html

--------------

Copyright (c) 2018 Pieter Moris Adrem Data Lab - biomina - UAntwerpen
