.. topicpy documentation master file, created by
   sphinx-quickstart on Mon Aug 17 16:28:20 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to topicpy's documentation!
===================================

Installation
------------

Install topicpy by running:

   pip install topicpy

This package consists of multiple modules:

- converter: handels gene name conversions
- geneontology: uses GSEA to perform gene ontologies
- gtex: handle GTEx data
- hsbmpy: handle ouput of hierarchical stochastic block model
- hypergeom: perform hypergeometric tests
- lda: Run sklearn LatenDirichletAllocation with some more parameters
- tableanalyser: study distributions of a RNA-Seq data
- TCGA_files: handle TCGA metadata

Resources
---------

* Online documentation: https://topicpy.readthedocs.io/en/latest/
* Paper: in preparation

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. toctree::
  :maxdepth: 2

  converter
  geneontology
  gtex
  hsbmpy
  hypergeom
  lda
  tableanalyser
  tacos_plot
  TCGA_files
