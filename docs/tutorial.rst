Tutorial
========

This tutorial walks through common workflows using the ``recount3`` Python
package.  All examples assume the package is installed::

   pip install "recount3[biocpy]"

For BigWig support add the ``bigwig`` extra as well::

   pip install "recount3[biocpy,bigwig]"


Discovering Resources for a Project
------------------------------------

:class:`~recount3.bundle.R3ResourceBundle` is the main entry point.
Call :meth:`~recount3.bundle.R3ResourceBundle.discover` to fetch every
resource associated with a project:

.. code:: python

   from recount3 import R3ResourceBundle

   bundle = R3ResourceBundle.discover(
       organism="human",
       data_source="sra",
       project="SRP009615",
   )
   print(f"Found {len(bundle.resources)} resources.")

Supported organisms are ``"human"`` and ``"mouse"``.  Supported data sources
are ``"sra"``, ``"gtex"``, and ``"tcga"``.


Stacking Count Matrices
-----------------------

Use :meth:`~recount3.bundle.R3ResourceBundle.stack_count_matrices` to
concatenate all gene-level (or exon-level, junction-level) count files
from a bundle into a single :class:`pandas.DataFrame`:

.. code:: python

   counts = bundle.only_counts().stack_count_matrices(genomic_unit="gene")
   print(counts.shape)   # (n_genes, n_samples)


Building a SummarizedExperiment
--------------------------------

Requires the ``biocpy`` extra (``biocframe``, ``summarizedexperiment``).

.. code:: python

   se = bundle.to_summarized_experiment(genomic_unit="gene")
   print(se)


Building a RangedSummarizedExperiment
---------------------------------------

A :class:`~summarizedexperiment.RangedSummarizedExperiment` attaches
genomic coordinate ranges to each feature row.  A GTF annotation file is
required; include it in the bundle via ``search_project_all`` (done
automatically by :meth:`~recount3.bundle.R3ResourceBundle.discover`) or
add an annotation resource manually.

.. code:: python

   rse = bundle.to_ranged_summarized_experiment(genomic_unit="gene")
   print(rse)
   print(rse.shape)          # (n_genes, n_samples)
   print(rse.colnames[:5])   # first five sample IDs


Accessing Metadata
------------------

Column metadata (sample annotations) is merged from every available
metadata table and aligned to the count columns automatically:

.. code:: python

   col_data = rse.coldata.to_pandas()
   print(col_data.columns.tolist())


Filtering Bundles
-----------------

Bundles are immutable snapshots.  Use :meth:`~recount3.bundle.R3ResourceBundle.filter`
to narrow down by any description field:

.. code:: python

   gene_only = bundle.filter(genomic_unit="gene")
   annots    = bundle.filter(resource_type="annotations")
   metadata  = bundle.only_metadata()


Working with Individual Resources
-----------------------------------

Each :class:`~recount3.resource.R3Resource` knows its URL and can be
downloaded and loaded independently:

.. code:: python

   from recount3 import R3GeneOrExonCounts, R3Resource

   desc = R3GeneOrExonCounts(
       organism="human",
       data_source="sra",
       genomic_unit="gene",
       project="SRP009615",
       sample="SRR387777",
   )
   res = R3Resource(desc)
   res.download(path=None, cache_mode="enable")  # saves to local cache
   df = res.load()                               # returns pd.DataFrame
   print(df.head())


Searching Without a Bundle
---------------------------

The :mod:`recount3.search` helpers return lists of
:class:`~recount3.resource.R3Resource` objects directly:

.. code:: python

   from recount3 import search_count_files_gene_or_exon

   resources = search_count_files_gene_or_exon(
       organism="human",
       data_source="sra",
       genomic_unit="gene",
       project="SRP009615",
   )
   print(f"{len(resources)} count files found.")


Cache Management
----------------

Downloaded files are cached under ``~/.cache/recount3/files`` by default.
Use the :mod:`recount3.config` helpers to inspect or clear the cache:

.. code:: python

   from recount3 import recount3_cache, recount3_cache_files, recount3_cache_rm

   print(recount3_cache())          # Path to the cache directory
   print(recount3_cache_files())    # List all cached files

   # Remove only junction files (dry run first):
   removed = recount3_cache_rm(
       predicate=lambda p: ".junctions." in p.name,
       dry_run=True,
   )
   print(f"Would remove {len(removed)} files.")

The cache directory and other settings can be overridden via environment
variables (``RECOUNT3_CACHE_DIR``, ``RECOUNT3_URL``, etc.) or by
constructing a custom :class:`~recount3.config.Config` and passing it to
any resource or search function.


SRA Attribute Expansion
-----------------------

recount3 stores SRA sample attributes as a single pipe-delimited string
column.  :func:`recount3.se.build_ranged_summarized_experiment` expands
them automatically, but you can also call the helper directly:

.. code:: python

   from recount3.se import expand_sra_attributes

   rse_expanded = expand_sra_attributes(rse)
   col_data = rse_expanded.coldata.to_pandas()
   # Columns like "sra_attribute.age", "sra_attribute.disease" are now present.
