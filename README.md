![Python application](https://github.com/fvalle1/topicpy/workflows/Python%20application/badge.svg)
![CodeQL](https://github.com/fvalle1/topicpy/workflows/CodeQL/badge.svg)
[![made-with-sphinx-doc](https://img.shields.io/badge/Made%20with-Sphinx-1f425f.svg)](https://topicpy.readthedocs.io/en/latest/)
[![GitHub](https://img.shields.io/github/license/fvalle1/topicpy)](LICENSE)

# topicpy
Python package for topic modelling pre and post processing

# Modules

This package consists of multiple modules:

- converter: handels gene name conversions
- geneontology: uses GSEA to perform gene ontologies
- gtex: handle GTEx data
- hsbmpy: handle ouput of hierarchical stochastic block model
- hypergeom: perform hypergeometric tests
- lda: Run sklearn LatenDirichletAllocation with some more parameters
- tableanalyser: study distributions of a RNA-Seq data
- TCGA_files: handle TCGA metadata

# Documentation
[https://topicpy.readthedocs.io/en/latest/](https://topicpy.readthedocs.io/en/latest/)

# Test
Run unittest with

``` python
python test/test.py
```

# License
See [LICENSE](LICENSE)
