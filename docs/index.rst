.. goscripts documentation master file, created by
   sphinx-quickstart on Sun May 13 15:17:19 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


goscripts - Python script and package for Gene Ontology enrichment analysis
===========================================================================

Welcome to goscripts's documentation!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
   procedures (provided by `statsmodels <http://www.statsmodels.org/stable/generated/statsmodels.sandbox.stats.multicomp.multipletests.html#statsmodels.sandbox.stats.multicomp.multipletests>`_): ``goscripts.enrichment_stats``

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   source/README
   source/goscripts

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
