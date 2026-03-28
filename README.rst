recount3
========

.. image:: https://img.shields.io/pypi/l/recount3.svg
   :target: https://github.com/MoseleyBioinformaticsLab/recount3/blob/main/LICENSE
   :alt: Clear BSD License with extra clauses

.. image:: https://img.shields.io/pypi/v/recount3.svg
   :target: https://pypi.org/project/recount3
   :alt: Current library version

.. image:: https://img.shields.io/pypi/pyversions/recount3.svg
   :target: https://pypi.org/project/recount3
   :alt: Supported Python versions

.. image:: https://codecov.io/gh/MoseleyBioinformaticsLab/recount3/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/MoseleyBioinformaticsLab/recount3
   :alt: CodeCov

.. image:: https://img.shields.io/github/stars/MoseleyBioinformaticsLab/recount3.svg?style=social&label=Star
   :target: https://github.com/MoseleyBioinformaticsLab/recount3
   :alt: GitHub project

|

The ``recount3`` package is a Python library and command-line tool for
interacting with the `recount3`_ data repository, a uniformly processed
collection of RNA-seq studies covering human and mouse samples from SRA,
GTEx, and TCGA.

The package provides a typed API for discovering, downloading, and assembling
recount3 resources into analysis-ready objects.  A companion CLI implements a
discover -> manifest -> materialize workflow suitable for use in scripts and
pipelines.

The ``recount3`` package can be used in two ways:

   * As a Python library for searching, downloading, and assembling
     recount3 gene/exon/junction count matrices, sample metadata, genome
     annotations, and BigWig coverage files into
     ``SummarizedExperiment`` / ``RangedSummarizedExperiment`` objects via
     the BiocPy stack.
   * As a command-line tool (``recount3``) to search the mirror, produce
     JSONL manifests, and materialize resources to disk or into a ``.zip``
     archive.


Links
~~~~~

   * recount3 @ GitHub_
   * recount3 @ PyPI_
   * Documentation @ Pages_
   * Questions & bug reports @ Issues_


Installation
~~~~~~~~~~~~

The core package requires Python 3.10 or newer and depends on NumPy, pandas,
and SciPy.  Optional extras unlock BiocPy integration and BigWig support.

Install on Linux, Mac OS X
--------------------------

.. code:: bash

   python3 -m pip install recount3

Install on Windows
------------------

.. code:: bash

   py -3 -m pip install recount3

Optional extras
---------------

Install with BiocPy support (``SummarizedExperiment``, ``RangedSummarizedExperiment``,
``GenomicRanges``):

.. code:: bash

   python3 -m pip install "recount3[biocpy]"

Install with BigWig support (``pybigwig``):

.. code:: bash

   python3 -m pip install "recount3[bigwig]"

Install all optional extras:

.. code:: bash

   python3 -m pip install "recount3[biocpy,bigwig]"

Upgrade on Linux, Mac OS X
---------------------------

.. code:: bash

   python3 -m pip install recount3 --upgrade

Upgrade on Windows
------------------

.. code:: bash

   py -3 -m pip install recount3 --upgrade


Quickstart
~~~~~~~~~~

Python API
----------

Discover all resources for a project, stack the gene-level count matrices
across samples, and build a ``RangedSummarizedExperiment``:

.. code:: python

   >>> from recount3 import R3ResourceBundle
   >>>
   >>> # Discover every resource for a human SRA project.
   >>> bundle = R3ResourceBundle.discover(
   ...     organism="human",
   ...     data_source="sra",
   ...     project="SRP009615",
   ... )
   >>>
   >>> # Stack raw gene-count matrices across all samples.
   >>> counts = bundle.only_counts().stack_count_matrices(genomic_unit="gene")
   >>>
   >>> # Build a RangedSummarizedExperiment (requires recount3[biocpy]).
   >>> rse = bundle.to_ranged_summarized_experiment(genomic_unit="gene")

Command-line tool
-----------------

Discover resources, save a JSONL manifest, and download in parallel:

.. code:: bash

   # Search for gene-level count files and write a manifest.
   recount3 search gene-exon \
       organism=human data_source=sra genomic_unit=gene project=SRP009615 \
       --format=jsonl > manifest.jsonl

   # Materialize all resources from the manifest (8 parallel jobs).
   recount3 download --from=manifest.jsonl --dest=./downloads --jobs=8

Stream search output directly into download without an intermediate file:

.. code:: bash

   recount3 search annotations \
       organism=human genomic_unit=gene annotation_extension=G026 \
       --format=jsonl | \
   recount3 download --from=- --dest=./annotations

.. note:: Read the full documentation on Pages_ for the complete API reference,
          CLI guide, and worked examples.


Dependencies
~~~~~~~~~~~~

Core (installed automatically):

.. code:: text

   numpy>=2.0
   pandas>=2.2
   scipy>=1.13

Optional: BiocPy integration (``recount3[biocpy]``):

.. code:: text

   biocframe>=0.7
   genomicranges>=0.8
   summarizedexperiment>=0.6

Optional: BigWig support (``recount3[bigwig]``):

.. code:: text

   pybigwig>=0.3.18


Questions, Feature Requests, and Bug Reports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please submit questions, feature requests, and bug reports on Issues_.


License
~~~~~~~

This package is distributed under the BSD_ license.


.. _recount3: https://rna.recount.bio
.. _BiocPy: https://github.com/BiocPy
.. _GitHub: https://github.com/MoseleyBioinformaticsLab/recount3
.. _PyPI: https://pypi.org/project/recount3
.. _Pages: https://moseleybioinformaticslab.github.io/recount3/
.. _Issues: https://github.com/MoseleyBioinformaticsLab/recount3/issues
.. _BSD: https://github.com/MoseleyBioinformaticsLab/recount3/blob/main/LICENSE
