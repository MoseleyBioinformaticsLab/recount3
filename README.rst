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

``recount3`` is a typed Python library and command-line tool for the `recount3`_
data repository, a uniformly processed collection of RNA-seq studies spanning
tens of thousands of human and mouse samples from SRA, GTEx, and TCGA. It
discovers, downloads, and assembles recount3 resources into analysis-ready
objects. These resources include gene, exon, and junction count matrices,
sample metadata, genome annotations, and BigWig coverage files.

The package provides two interfaces:

* A Python library that assembles count matrices, sample metadata, and
  genomic coordinates into BiocPy ``SummarizedExperiment`` and
  ``RangedSummarizedExperiment`` objects, with recount3-compatible scaling and
  normalization utilities (approximate read counts, AUC- and
  mapped-reads-based scaling, and TPM).
* A command-line tool (``recount3``) that implements a
  discover -> manifest -> materialize workflow for scripts and pipelines. It
  emits JSONL/TSV manifests and materializes resources to a directory or a
  ``.zip`` archive with parallel downloads.


Links
~~~~~

   * recount3 @ GitHub_
   * recount3 @ PyPI_
   * Documentation @ Pages_
   * Questions & bug reports @ Issues_


Installation
~~~~~~~~~~~~

The core package requires Python 3.10 or newer and depends on NumPy, pandas,
and SciPy. Two optional extras enable additional features:

.. code:: bash

   python3 -m pip install recount3                   # core
   python3 -m pip install "recount3[biocpy]"         # + SummarizedExperiment builders
   python3 -m pip install "recount3[bigwig]"         # + BigWig coverage access
   python3 -m pip install "recount3[biocpy,bigwig]"  # everything

* ``biocpy`` (``biocframe``, ``genomicranges``, ``summarizedexperiment``) is
  required for ``create_rse`` and every helper that returns or operates on a
  BiocPy object.
* ``bigwig`` (``pyBigWig``) is required only for BigWig coverage access.

On Windows, substitute ``py -3 -m pip install ...``. Upgrade an existing
installation with ``python3 -m pip install --upgrade recount3``.

**Note:** The optional extras have platform constraints. The ``bigwig`` extra
can be difficult or impossible to install on Windows and macOS. The ``biocpy``
extra can be difficult or impossible to install on Windows. The core package and
the command-line workflow do not depend on either extra and work on all
supported platforms.


The three-layer API
~~~~~~~~~~~~~~~~~~~~~

``recount3`` exposes the same workflow at three levels of abstraction, so a
single project can be assembled in one call while multi-project or custom
workflows retain full control:

* **High level.** ``create_rse()`` builds one project into a
  ``RangedSummarizedExperiment`` and performs discovery, downloading, metadata
  merging, and range assembly in a single call.
* **Mid level.** ``R3ResourceBundle`` is a filterable container of resources
  for combining multiple projects, selecting subsets, and stacking matrices.
* **Low level.** ``R3Resource`` represents a single file and manages its URL,
  cache entry, and parser.

See the Tutorial_ for a complete walkthrough.


Quickstart
~~~~~~~~~~

Python API
----------

Assemble a project into a ``RangedSummarizedExperiment`` (requires the
``recount3[biocpy]`` extra):

.. code:: python

   >>> from recount3 import create_rse
   >>>
   >>> rse = create_rse(
   ...     project="SRP009615",
   ...     organism="human",
   ...     annotation_label="gencode_v26",
   ... )
   >>> rse.shape
   (63856, 12)

For multi-project or custom workflows, use the bundle layer to filter
resources and stack matrices directly:

.. code:: python

   >>> from recount3 import R3ResourceBundle
   >>>
   >>> bundle = R3ResourceBundle.discover(
   ...     organism="human",
   ...     data_source="sra",
   ...     project="SRP009615",
   ... )
   >>> print(f"Found {len(bundle.resources)} resources.")
   Found 10 resources.
   >>>
   >>> gene_counts = bundle.filter(
   ...     resource_type="count_files_gene_or_exon",
   ...     genomic_unit="gene",
   ... ).stack_count_matrices(compat="feature")
   >>> gene_counts.shape
   (63856, 12)

Command-line tool
-----------------

Discover resources, write a JSONL manifest, and download in parallel:

.. code:: bash

   # Search for gene-level count files and write a manifest.
   recount3 search gene-exon \
       organism=human data_source=sra genomic_unit=gene project=SRP009615 \
       --format=jsonl > manifest.jsonl

   # Materialize all resources from the manifest (8 parallel jobs).
   recount3 download --from=manifest.jsonl --dest=./downloads --jobs=8

Because both subcommands operate on JSONL via standard streams, search and
download compose into a single pipeline without an intermediate file:

.. code:: bash

   recount3 search annotations \
       organism=human genomic_unit=gene annotation_extension=G026 \
       --format=jsonl | \
   recount3 download --from=- --dest=./annotations

The ``bundle`` subcommands assemble analysis-ready outputs without a Python
session. Supported outputs are a stacked count matrix (TSV, gzip-compressed
TSV, or Parquet) and a pickled ``SummarizedExperiment`` or
``RangedSummarizedExperiment``:

.. code:: bash

   recount3 bundle rse --from=manifest.jsonl --genomic-unit=gene --out=rse.pkl

**Note:** Read the full documentation on Pages_ for the complete API reference,
the CLI guide, and worked examples.


Data mirrors
~~~~~~~~~~~~

recount3 publishes the same relative file layout on several interchangeable
public mirrors. ``recount3`` targets this layout rather than any single host, so
selecting a different mirror requires only a change to the base URL (the
``RECOUNT3_URL`` environment variable, the ``--base-url`` CLI flag, or the
``base_url`` field of ``recount3.config.Config``):

==============================  =======================================================
Mirror                          Base URL
==============================  =======================================================
Duffel load balancer (default)  ``http://duffel.rail.bio/recount3/``
AWS Open Data                   ``https://recount-opendata.s3.amazonaws.com/recount3/release/``
JHU IDIES (Dataverse)           ``https://data.idies.jhu.edu/recount3/data/``
==============================  =======================================================


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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please submit questions, feature requests, and bug reports on Issues_.


License
~~~~~~~

This package is distributed under the BSD_ license.


.. _recount3: https://rna.recount.bio
.. _BiocPy: https://github.com/BiocPy
.. _GitHub: https://github.com/MoseleyBioinformaticsLab/recount3
.. _PyPI: https://pypi.org/project/recount3
.. _Pages: https://moseleybioinformaticslab.github.io/recount3/
.. _Tutorial: https://moseleybioinformaticslab.github.io/recount3/tutorial.html
.. _Issues: https://github.com/MoseleyBioinformaticsLab/recount3/issues
.. _BSD: https://github.com/MoseleyBioinformaticsLab/recount3/blob/main/LICENSE
