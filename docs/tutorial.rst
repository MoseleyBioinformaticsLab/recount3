Tutorial
========

This tutorial walks through the ``recount3`` Python API end-to-end: resource
discovery, downloading, assembly of
:class:`~summarizedexperiment.SummarizedExperiment` /
:class:`~summarizedexperiment.RangedSummarizedExperiment` objects, sample
metadata merging, count normalization and scaling, BigWig coverage access, and
management of the on-disk cache.

For the ``recount3`` command-line tool, see :doc:`cli`. For full per-symbol
documentation, see :doc:`api`.


Installation
------------

The core package depends only on NumPy, pandas, and SciPy. The optional
extras enable features used throughout this tutorial:

.. code:: bash

   python3 -m pip install recount3               # core only
   python3 -m pip install "recount3[biocpy]"     # + SummarizedExperiment
   python3 -m pip install "recount3[bigwig]"     # + pyBigWig
   python3 -m pip install "recount3[parquet]"    # + .parquet output
   python3 -m pip install "recount3[anndata]"    # + .h5ad output
   python3 -m pip install "recount3[all]"        # everything

What each extra enables:

- ``biocpy`` is required for :func:`~recount3.create_rse`,
  :meth:`~recount3.bundle.R3ResourceBundle.to_summarized_experiment`,
  :meth:`~recount3.bundle.R3ResourceBundle.to_ranged_summarized_experiment`,
  and every helper in :mod:`recount3.se` that returns or operates on a
  BiocPy object.
- ``bigwig`` is required only when you call
  :meth:`~recount3.resource.R3Resource.load` on a BigWig resource or use the
  :class:`~recount3._bigwig.BigWigFile` reader directly.
- ``parquet`` installs ``pyarrow`` and is required to write ``.parquet``, both
  from :meth:`pandas.DataFrame.to_parquet` on a stacked matrix and from
  ``recount3 bundle stack-counts --out=counts.parquet``. pandas accepts either
  ``pyarrow`` or ``fastparquet``; an existing ``fastparquet`` installation is
  used as-is, and the ``io.parquet.engine`` option is honoured.
- ``anndata`` installs ``anndata`` and ``delayedarray`` (and implies
  ``biocpy``). ``SummarizedExperiment.to_anndata()`` imports both, and
  ``summarizedexperiment`` declares neither as a required dependency, so the
  extra is needed to write ``.h5ad`` from ``recount3 bundle se`` /
  ``recount3 bundle rse``.

If an optional dependency is missing, the affected function raises
:exc:`ImportError` on first use; the remainder of the package stays importable
and functional. The two output extras are additionally checked up front by the
CLI, before any download runs, so an unusable output format is reported
immediately rather than after the data has been fetched and assembled.


Quick start
-----------

.. note::

   Every example in this tutorial retrieves data from a live recount3 mirror
   and therefore requires network access. Downloaded files are cached under
   ``~/.cache/recount3/files`` (see :ref:`cache-and-configuration`), so
   re-running an example reuses the local copy rather than downloading again.

The most direct path from a project identifier to an analysis-ready BiocPy
object is :func:`~recount3.create_rse`. It requires the ``biocpy`` extra and
performs discovery, downloads, metadata merging, and range assembly in a single
call:

.. code:: python

   from recount3 import create_rse

   rse = create_rse(
       project="SRP009615",
       organism="human",
       annotation_label="gencode_v26",
   )

   print(rse.shape)             # (n_features, n_samples)
   print(rse.get_column_names()[:5])

This single call is sufficient for the most common workflow; it is expanded in
:ref:`layer-1` below. The remainder of this tutorial describes the steps that
``create_rse`` performs internally and the lower-level components to use when
finer control is required.


The three layers of the API
---------------------------

``recount3`` exposes the same workflow at three levels of abstraction:

.. list-table::
   :header-rows: 1
   :widths: 22 26 52

   * - Layer
     - Primary entry point
     - Recommended when
   * - High-level: BiocPy builders
     - :func:`~recount3.create_rse`
     - You want one project as a ``RangedSummarizedExperiment``.
   * - Mid-level: bundles
     - :class:`~recount3.R3ResourceBundle`
     - You combine multiple projects, filter resources, or stack matrices yourself.
   * - Low-level: resources
     - :class:`~recount3.R3Resource`
     - You want fine-grained control over a single file's URL, download, or parser.

Each layer is a thin wrapper around the next. ``create_rse`` calls
``R3ResourceBundle.discover`` internally; ``R3ResourceBundle`` aggregates
``R3Resource`` objects. Because the layers share a common set of types, they
interoperate freely: a bundle obtained from ``discover`` can be filtered at
Layer 2 and then handed to the same builders that ``create_rse`` invokes.


.. _layer-1:

Layer 1: Building experiments with ``create_rse``
-------------------------------------------------

:func:`~recount3.create_rse` is the recommended entry point for the most
common workflow: one project, one organism, one annotation, one assembled
:class:`~summarizedexperiment.RangedSummarizedExperiment`. Requires the
``biocpy`` extra.

**Use this layer when** you have a study accession and want its expression
data ready to analyze. You read a paper that used ``SRP009615`` and want to
re-run the differential expression yourself; you need gene-level counts for
one GTEx tissue to test a hypothesis; you are checking whether a gene of
interest is expressed in a public dataset before designing an experiment. In
each case one accession goes in and one object comes out, with counts,
sample metadata, and genomic coordinates already aligned to each other.

**Use a different layer when** one project is not the unit of work: reach
for Layer 2 if you need several studies in one matrix or want to inspect the
files before committing to a download, and Layer 3 if you want one specific
file and nothing else.

Gene-level RSE (default)
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   from recount3 import create_rse

   rse = create_rse(
       project="SRP009615",
       organism="human",
       annotation_label="gencode_v26",   # or "gencode_v29", "fantom6_cat", "refseq", "ercc", "sirv"
   )

You may pass the raw extension code instead of a label:

.. code:: python

   rse = create_rse(
       project="SRP009615",
       organism="human",
       annotation_extension="G026",
   )

When both are supplied, ``annotation_extension`` takes precedence. Discover the
available labels with :func:`~recount3.annotation_options`:

.. code:: python

   from recount3 import annotation_options

   annotation_options("human")
   # {'gencode_v26': 'G026', 'gencode_v29': 'G029', 'fantom6_cat': 'F006',
   #  'refseq': 'R109', 'ercc': 'ERCC', 'sirv': 'SIRV'}

   annotation_options("mouse")
   # {'gencode_v23': 'M023'}

Exon-level and junction-level
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   exon_rse = create_rse(
       project="SRP009615",
       organism="human",
       genomic_unit="exon",
       annotation_label="gencode_v26",
   )

   junction_rse = create_rse(
       project="SRP009615",
       organism="human",
       genomic_unit="junction",
   )

For junctions, ``recount3`` prefers the RR sidecar for genomic
coordinates; pass ``prefer_rr_junction_coordinates=False`` to disable
this.

Falling back to a plain ``SummarizedExperiment``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every gene or exon row range comes from an annotation GTF, and every
junction range from an RR sidecar. When that file cannot be turned into
ranges, ``create_rse`` raises :exc:`~recount3.RangesError` (a
:exc:`ValueError`) naming which of three things went wrong:

``the annotation could not be retrieved``
   The download failed. The HTTP layer has already retried it
   (``RECOUNT3_MAX_RETRIES``, default 3), so this points at a mirror or
   network that is actually down rather than a momentary blip. Try again
   later, or point ``RECOUNT3_URL`` at another mirror.

``the annotation could not be parsed``
   The file arrived but is not readable as a GTF, most often a truncated
   cache entry from an interrupted download. Drop it with
   :func:`~recount3.recount3_cache_rm` and let it download again.

``the annotation does not cover every counted feature``
   The annotation parses cleanly but describes a different feature set
   than the counts -- a GENCODE 26 GTF against GENCODE 29 counts, say.
   Pass the matching ``annotation_extension`` or ``annotation_label``;
   :func:`~recount3.annotation_options` lists what is available.

You do not have to predict which of these will happen before you call.
Ask for the RSE; you either get one, or you get a message naming the
cause. Nothing is silently degraded in between.

``allow_fallback_to_se=True`` changes only what happens in that failure
case: rather than raising, you get a plain
:class:`~summarizedexperiment.SummarizedExperiment` and a logged warning
explaining why:

.. code:: python

   experiment = create_rse(
       project="SRP009615",
       organism="human",
       allow_fallback_to_se=True,
   )

Be deliberate about that flag:

- The returned object **has no genomic ranges**. Counts, sample metadata,
  and the operations driven by column data still work --
  :func:`recount3.se.transform_counts`,
  :func:`recount3.se.compute_scale_factors`,
  :func:`recount3.se.expand_sra_attributes`. Anything needing coordinates
  does not: range queries, overlaps, subsetting by region, and
  :func:`recount3.se.compute_tpm`, which needs feature widths and raises on
  a plain SE. Code written for an RSE will still fail, just later and
  further from the cause.
- It is **not a retry**. The download is not attempted again, and the flag
  does nothing about whatever made retrieval fail.
- It does **not repair an annotation mismatch**. Accepting an SE is how a
  wrong ``annotation_extension`` goes unnoticed.

It earns its place in count-only work, and in batch jobs that should not
abort on one bad project. When you do want ranges, fix the cause the
message names instead of passing the flag.

Operations performed by ``create_rse``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For one ``(organism, data_source, project)`` triple, it:

1. discovers gene/exon/junction counts, the matching annotation GTF, and
   the five project metadata tables;
2. downloads everything into the on-disk cache;
3. stacks the count matrix into a feature × sample DataFrame;
4. merges, namespaces, and aligns the metadata tables to the count columns
   (including a ``BigWigURL`` column constructed per sample);
5. parses the GTF (or RR file, for junctions) to attach genomic ranges;
6. constructs the BiocPy object.

To deviate from any of those steps (multiple projects, custom metadata
filtering, stacking only some matrices, or a different join policy), use
Layer 2.


Layer 2: Resource bundles
-------------------------

:class:`~recount3.R3ResourceBundle` is a container of
:class:`~recount3.R3Resource` objects with helpers for filtering,
loading, stacking, and converting to BiocPy objects.

**Use this layer when** ``create_rse``'s one-project, everything-at-once
shape does not fit. Typical cases:

- *Several studies, one matrix.* You want a gene-level matrix spanning
  ``SRP009615`` and ``SRP001558`` to look for an effect that holds across
  both, so the counts have to be stacked on a shared feature axis before
  any analysis starts.
- *Look before you download.* Junction and BigWig files are large. A bundle
  lists what exists for a project -- URLs, sizes, annotation codes -- as
  plain objects, so you can decide what is worth fetching, or write out a
  manifest for someone else to fetch.
- *Only part of what discovery returns.* You need the four gene-count files
  and the QC table, not the exon counts, the junctions, and the other four
  metadata tables that ``create_rse`` would download alongside them.
- *Your own assembly.* You want the stacked counts as a
  :class:`pandas.DataFrame` to merge with clinical data of your own, and
  will build the BiocPy object yourself (or skip it entirely).

A bundle is just a list of resources plus filters, so nothing is downloaded
until you ask for it. That is the practical difference from Layer 1:
``create_rse`` decides what to fetch for you, a bundle lets you decide.

Discovering resources for one or more projects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   from recount3 import R3ResourceBundle

   bundle = R3ResourceBundle.discover(
       organism="human",
       data_source="sra",
       project="SRP009615",
   )
   print(f"Found {len(bundle.resources)} resources.")
   # Found 10 resources.

By default, ``discover`` returns gene + exon counts, the default annotation
GTF for each of those units (gene and exon), the five metadata tables, and the
default junction artifact (``MM``). For ``SRP009615`` this is the ten resources
counted above: 2 counts + 2 annotation GTFs + 5 metadata tables + 1 junction
file. Override with:

.. code:: python

   bundle = R3ResourceBundle.discover(
       organism="human",
       data_source="sra",
       project="SRP009615",
       genomic_units=("gene",),
       annotations=("G026", "G029"),       # or "default" / "all"
       junction_exts=("MM", "RR"),
       include_metadata=True,
       include_bigwig=False,
   )

Multi-project bundles
~~~~~~~~~~~~~~~~~~~~~

Pass an iterable for any of ``organism``, ``data_source``, or
``project``. ``discover`` computes the Cartesian product and produces a
single combined bundle:

.. code:: python

   multi = R3ResourceBundle.discover(
       organism="human",
       data_source="sra",
       project=["SRP009615", "SRP001558"],
       genomic_units=("gene",),
   )
   print(f"Combined: {len(multi.resources)} resources, 2 projects.")
   # Combined: 15 resources, 2 projects.

The count is 15, not 14: each project contributes 7 project-specific resources
(1 gene count + 1 junction file + 5 metadata tables), and the gene annotation
GTF is shared across both projects, so it is deduplicated to a single resource
(7 × 2 + 1 = 15). Note that the junction artifact is included by default
regardless of ``genomic_units``; pass ``junction_exts=()`` to omit it.

When a bundle spans more than one ``(organism, data_source, project)``
triple, its ``organism``/``data_source``/``project`` attributes are left
as ``None`` to avoid misrepresenting its identity; per-resource fields
remain authoritative.

Supported values:

- ``organism``: ``"human"``, ``"mouse"``
- ``data_source``: ``"sra"``, ``"gtex"``, ``"tcga"``

Filtering bundles
~~~~~~~~~~~~~~~~~

Bundles are returned by-value from :meth:`~recount3.R3ResourceBundle.filter`;
the original is not mutated. Each keyword maps to a field on the
underlying :class:`~recount3._descriptions.R3ResourceDescription`, and
accepts any :data:`~recount3.types.FieldSpec`:

- a single string: exact match
- an iterable of strings: membership test
- a callable ``(value) -> bool``: predicate

.. code:: python

   gene_counts = bundle.filter(
       resource_type="count_files_gene_or_exon",
       genomic_unit="gene",
   )

   gene_or_exon = bundle.filter(genomic_unit=["gene", "exon"])

   gencode_only = bundle.filter(
       annotation_extension=lambda ext: ext and ext.startswith("G"),
   )

   no_metadata = bundle.filter(resource_type="metadata_files", invert=True)

Convenience aliases provide shortcuts for the most common filters:
:meth:`~recount3.R3ResourceBundle.only_counts`,
:meth:`~recount3.R3ResourceBundle.only_metadata`,
:meth:`~recount3.R3ResourceBundle.bigwigs`,
:meth:`~recount3.R3ResourceBundle.exclude_metadata`.

.. note::

   Filtering on a field that a resource does not have (for example,
   filtering on ``genomic_unit`` when metadata files have no genomic unit)
   excludes those resources from the result. Combine filters explicitly
   when this matters: ``bundle.filter(resource_type=..., genomic_unit=...)``.

Stacking count matrices
~~~~~~~~~~~~~~~~~~~~~~~

:meth:`~recount3.R3ResourceBundle.stack_count_matrices` concatenates count
DataFrames. It does not take a ``genomic_unit`` argument, so filter the
bundle first to choose which family you want:

.. code:: python

   gene_counts_df = (
       bundle
       .filter(resource_type="count_files_gene_or_exon", genomic_unit="gene")
       .stack_count_matrices(compat="feature")
   )
   print(gene_counts_df.shape)                        # (n_features, n_samples)
   # (63856, 12)

   junction_counts_df = (
       bundle
       .filter(resource_type="count_files_junctions", junction_extension="MM")
       .stack_count_matrices()
   )

Compatibility checking is controlled by ``compat``:

- ``compat="family"`` (default): gene/exon may mix with gene/exon;
  junctions stay with junctions.
- ``compat="feature"``: stricter; the feature space must match exactly: the
  same genomic unit (gene versus exon) for gene/exon counts, or the same
  junction subtype for junctions. (The annotation build is not constrained.)

Mixing incompatible resources raises
:exc:`~recount3.errors.CompatibilityError`.

Building SummarizedExperiment / RangedSummarizedExperiment from a bundle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The bundle methods below are what ``create_rse`` calls internally:

.. code:: python

   se = bundle.to_summarized_experiment(genomic_unit="gene")

   rse = bundle.to_ranged_summarized_experiment(
       genomic_unit="gene",
       annotation_extension="G026",
       allow_fallback_to_se=False,
   )

The same functions are available as standalone wrappers in
:mod:`recount3.se` (:func:`~recount3.build_summarized_experiment`,
:func:`~recount3.build_ranged_summarized_experiment`) for symmetry with
``create_rse``.

Downloading a bundle's files in parallel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:meth:`~recount3.R3ResourceBundle.download` materializes every resource in a
bundle to a local destination. Because retrieval is I/O-bound, resources are
fetched concurrently by a pool of worker threads sized by ``max_workers``
(default 8). This is the same mechanism used by the ``recount3 download``
command-line tool:

.. code:: python

   bundle.download(dest="./downloads", max_workers=8)

``dest`` may be a directory (each resource written as a separate file) or a
path ending in ``.zip`` (resources written into a single archive). The
``cache`` keyword (named ``cache_mode`` on
:meth:`~recount3.R3Resource.download`) accepts the same values: ``"enable"``,
``"update"``, ``"disable"``.


Layer 3: Individual resources
-----------------------------

:class:`~recount3.R3Resource` is the lowest level: one file, one URL, one
cache entry, one parser.

**Use this layer when** the unit of work is a single file. You want the
GENCODE 26 gene annotation itself, not an experiment built from it; you are
mirroring a handful of known URLs into a shared directory for a cluster job;
you want the raw junction MatrixMarket file to feed a tool of your own; or
you are debugging and want to see exactly which URL a description resolves
to before anything is fetched.

A resource is built from a description. Descriptions are typed
dataclasses with field validation; the recommended constructor is the
:class:`~recount3.R3ResourceDescription` factory, which routes to the
appropriate subclass based on ``resource_type``:

.. code:: python

   from recount3 import R3Resource, R3ResourceDescription

   desc = R3ResourceDescription(
       resource_type="count_files_gene_or_exon",
       organism="human",
       data_source="sra",
       genomic_unit="gene",
       project="SRP009615",
       annotation_extension="G026",   # required for gene/exon counts
   )

   res = R3Resource(desc)
   print(res.url)                     # fully-qualified URL on the recount3 mirror
   # http://duffel.rail.bio/recount3/human/data_sources/sra/gene_sums/15/SRP009615/sra.gene_sums.SRP009615.G026.gz

   res.download(path=None, cache_mode="enable")  # cache only, no local copy
   df = res.load()                               # pandas.DataFrame
   print(df.shape)

The full description catalog:

==========================================  =================================================
Resource type                               Description class
==========================================  =================================================
``"annotations"``                           :class:`~recount3.R3Annotations`
``"count_files_gene_or_exon"``              :class:`~recount3.R3GeneOrExonCounts`
``"count_files_junctions"``                 :class:`~recount3.R3JunctionCounts`
``"metadata_files"``                        :class:`~recount3.R3ProjectMetadata`
``"bigwig_files"``                          :class:`~recount3.R3BigWig`
``"data_sources"``                          :class:`~recount3.R3DataSources`
``"data_source_metadata"``                  :class:`~recount3.R3DataSourceMetadata`
==========================================  =================================================

Downloading
~~~~~~~~~~~

:meth:`~recount3.R3Resource.download` has three forms, controlled by
``path``:

.. code:: python

   res.download(path=None)                       # cache only
   res.download(path="/data/recount3")           # copy into a directory
   res.download(path="/data/recount3.zip")       # append to a ZIP archive

``cache_mode`` controls cache interaction:

- ``"enable"`` (default): use cached copy if present; download if not.
- ``"update"``: force a fresh download, then overwrite the cache.
- ``"disable"``: bypass the cache entirely (only valid when ``path`` is
  a directory or ``.zip``).

Loading
~~~~~~~

:meth:`~recount3.R3Resource.load` parses the cached file. The return type
depends on the resource:

================================================  =========================================
Resource type                                     ``load()`` returns
================================================  =========================================
Gene/exon counts                                  :class:`pandas.DataFrame`
Junction MM (with ID sidecar)                     :class:`pandas.DataFrame` (sparse-backed)
Junction ID / RR                                  :class:`pandas.DataFrame`
Metadata tables / source listings                 :class:`pandas.DataFrame`
BigWig                                            :class:`~recount3._bigwig.BigWigFile`
================================================  =========================================

The parsed object is cached on the resource; subsequent ``load()`` calls
return the same instance until you call
:meth:`~recount3.R3Resource.clear_loaded` (or pass ``force=True``).

Searching without a bundle
~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want a flat list of resources rather than a bundle, the
:mod:`recount3.search` helpers return ``list[R3Resource]`` directly. Each
takes :data:`~recount3.types.StringOrIterable` for every parameter and
returns one resource per Cartesian-product combination:

.. code:: python

   from recount3 import (
       search_count_files_gene_or_exon,
       search_metadata_files,
       search_bigwig_files,
   )

   counts = search_count_files_gene_or_exon(
       organism="human",
       data_source="sra",
       genomic_unit="gene",
       project="SRP009615",
       annotation_extension="G026",
   )

   meta = search_metadata_files(
       organism="human",
       data_source="sra",
       project="SRP009615",
       table_name=("recount_project", "recount_qc", "recount_seq_qc",
                   "recount_pred", "sra"),
   )

   bigwigs = search_bigwig_files(
       organism="human",
       data_source="sra",
       project="SRP009615",
       sample=["SRR387777", "SRR387778"],
   )

The single-call equivalent is :func:`~recount3.search_project_all`
(used internally by ``R3ResourceBundle.discover``).


Working with sample metadata
----------------------------

When ``create_rse`` or
:meth:`~recount3.R3ResourceBundle.to_ranged_summarized_experiment`
assembles an RSE, it merges all available per-project metadata tables
into ``column_data``, namespacing non-key columns by their table of
origin (e.g. ``recount_qc__star__all_mapped_reads``).

Access it as a pandas DataFrame:

.. code:: python

   col_df = rse.get_column_data().to_pandas()
   col_df.columns[:10]

Expanding SRA sample attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In an assembled RSE, SRA samples carry an ``sra__sample_attributes`` column
(the ``sra`` metadata table namespaced with ``__`` as described above) that
encodes key–value pairs in the form ``"age;;67.78|disease;;Control|..."``.
:func:`recount3.se.expand_sra_attributes` parses these into separate columns.
It accepts either a DataFrame or an SE/RSE object, and recognizes both the
namespaced ``sra__sample_attributes`` name and the R-style
``sra.sample_attributes`` spelling. Each parsed attribute becomes a new column
named ``sra_attribute.<key>`` (for example, ``sra_attribute.disease``):

.. code:: python

   from recount3.se import expand_sra_attributes

   rse2 = expand_sra_attributes(rse)
   col_df = rse2.get_column_data().to_pandas()
   sra_cols = [c for c in col_df.columns if c.startswith("sra_attribute.")]


Normalization and scaling
-------------------------

recount3 distributes coverage-sum counts ("``raw_counts``" assay), not
read counts. :mod:`recount3.se` provides recount3-compatible helpers to
convert and normalize them. All require a
:class:`~summarizedexperiment.RangedSummarizedExperiment`.

Approximate read counts
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   from recount3.se import compute_read_counts

   reads = compute_read_counts(rse)        # pandas DataFrame, integer-rounded

Values are rounded to whole reads by default; pass ``round_to_integers=False``
to retain the fractional estimates.

Per-sample scale factors
~~~~~~~~~~~~~~~~~~~~~~~~

Two methods are supported, matching the R ``recount3`` reference:

.. code:: python

   from recount3.se import compute_scale_factors, transform_counts

   sf_auc      = compute_scale_factors(rse, by="auc")
   sf_mapreads = compute_scale_factors(rse, by="mapped_reads")

Apply scale factors to the assay:

.. code:: python

   scaled = transform_counts(rse, by="auc")          # default
   scaled = transform_counts(rse, by="mapped_reads", target_read_count=4e7)

TPM (gene/exon only, needs feature widths)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   from recount3.se import compute_tpm

   tpm = compute_tpm(rse)                  # requires rowRanges with widths

:func:`recount3.se.is_paired_end` and the other helpers documented in
:mod:`recount3.se` accept either a DataFrame of metadata or an SE/RSE
object. See the API reference for full signatures.


BigWig coverage
---------------

Per-sample BigWig coverage files are not included by default; pass
``include_bigwig=True`` (or use the ``search_bigwig_files`` helper)
to add them. Requires the ``bigwig`` extra.

.. code:: python

   bundle = R3ResourceBundle.discover(
       organism="human",
       data_source="sra",
       project="SRP009615",
       include_bigwig=True,
   )

   for res, bw in bundle.iter_bigwig():
       with bw:
           print(res.description.sample, bw.chroms("chr1"))
           mean_chr1 = bw.stats("chr1", 0, 1_000_000, type="mean")[0]

:class:`~recount3._bigwig.BigWigFile` is a thin wrapper around
``pyBigWig``. Its main methods are ``chroms()``, ``header()``, ``values()``,
``stats()``, ``intervals()``, and ``close()``:

.. code:: python

   bw_res = bundle.bigwigs().resources[0]
   bw = bw_res.load()                       # BigWigFile
   with bw:                                 # closes the handle on exit
       values = bw.values("chr1", 0, 1000, numpy=True)


.. _cache-and-configuration:

Cache and configuration
-----------------------

Downloaded files are stored under ``~/.cache/recount3/files`` by default.
The :mod:`recount3.config` helpers let you inspect and prune the cache:

.. code:: python

   from recount3 import (
       recount3_cache,
       recount3_cache_files,
       recount3_cache_rm,
   )

   print(recount3_cache())                 # cache directory Path
   files = recount3_cache_files(pattern="*.gtf.gz")

   # Dry-run a deletion first:
   to_remove = recount3_cache_rm(
       predicate=lambda p: ".junctions." in p.name,
       dry_run=True,
   )
   recount3_cache_rm(predicate=lambda p: ".junctions." in p.name)

Configuration precedence is, from lowest to highest: library defaults,
environment variables, an explicit :class:`~recount3.Config` passed to a
resource or search function. The supported environment variables are:

==============================  ========================================
Variable                        Effect
==============================  ========================================
``RECOUNT3_URL``                Base URL of the recount3 mirror
``RECOUNT3_CACHE_DIR``          On-disk cache directory
``RECOUNT3_CACHE_DISABLE``      ``"1"`` to disable caching
``RECOUNT3_HTTP_TIMEOUT``       Network timeout (seconds)
``RECOUNT3_MAX_RETRIES``        Transient-error retry attempts
``RECOUNT3_INSECURE_SSL``       ``"1"`` to skip TLS verification (unsafe)
``RECOUNT3_USER_AGENT``         Custom ``User-Agent`` header
``RECOUNT3_CHUNK_SIZE``         Streaming chunk size (bytes)
==============================  ========================================

.. note::

   recount3 publishes the same file layout on several interchangeable public
   mirrors, so ``RECOUNT3_URL`` may point at any of them: the Duffel load
   balancer (``http://duffel.rail.bio/recount3/``, the default), AWS Open Data
   (``https://recount-opendata.s3.amazonaws.com/recount3/release/``), or JHU
   IDIES (``https://data.idies.jhu.edu/recount3/data/``). ``RECOUNT3_INSECURE_SSL``
   affects only ``https`` mirrors; it is a no-op for the default ``http`` mirror.

For programmatic use, construct a :class:`~recount3.Config` and pass it
explicitly:

.. code:: python

   from pathlib import Path
   from recount3 import Config, R3Resource, R3GeneOrExonCounts

   cfg = Config(
       base_url="http://duffel.rail.bio/recount3/",
       timeout=60,
       insecure_ssl=False,
       max_retries=5,
       user_agent="my-pipeline/0.1",
       cache_dir=Path("/scratch/recount3_cache"),
       cache_disabled=False,
       chunk_size=1024 * 1024,
   )

   res = R3Resource(
       R3GeneOrExonCounts(
           organism="human", data_source="sra", genomic_unit="gene",
           project="SRP009615", annotation_extension="G026",
       ),
       config=cfg,
   )


Errors and troubleshooting
--------------------------

All ``recount3`` exceptions derive from :exc:`~recount3.Recount3Error`,
so a single ``except`` clause catches every package-specific failure:

================================================  ======================================================
Exception                                         Raised when
================================================  ======================================================
:exc:`~recount3.ConfigurationError`               Bad config (env var, cache dir, option combinations)
:exc:`~recount3.DownloadError`                    Network/I-O failure during download
:exc:`~recount3.LoadError`                        Cached file parsed empty, malformed, or shape-mismatched
:exc:`~recount3.CompatibilityError`               Incompatible resources combined in a stack/build
:exc:`~recount3.RangesError`                      Genomic ranges could not be derived for an RSE
:exc:`~recount3.MissingRangesError`               Nothing in the bundle can supply genomic ranges
:exc:`~recount3.RangesCoverageError`              Ranges source omits some counted features
================================================  ======================================================

The three ranges errors also subclass :exc:`ValueError`, which is what
``create_rse`` has always raised on this failure, so existing
``except ValueError`` handlers keep working. ``MissingRangesError`` and
``RangesCoverageError`` are what you catch to tell "there was no
annotation to read" from "the annotation was the wrong one".

Common pitfalls
~~~~~~~~~~~~~~~

``ImportError: summarizedexperiment is required``
   Install the BiocPy extra: ``pip install "recount3[biocpy]"``.

``Writing Parquet requires a Parquet engine``
   No ``pyarrow`` or ``fastparquet`` is installed. Run
   ``pip install "recount3[parquet]"``, or write ``.tsv``, ``.tsv.gz``, or
   ``.csv`` instead.

``Cannot write .h5ad: Optional dependency 'anndata' is required``
   ``SummarizedExperiment.to_anndata()`` needs ``anndata`` and
   ``delayedarray``. Run ``pip install "recount3[anndata]"``, or write a
   ``.pkl`` instead.

``Cannot write .h5ad: N column name(s) contain a forward slash``
   HDF5 reads ``/`` as a path separator, and recount3 STAR QC fields are
   named after splice motifs (``..._gt/ag``). Pass ``--sanitize-columns``
   to rename them to ``..._gt_ag``, or write a ``.pkl`` instead, which
   keeps the names verbatim.

``Cannot write Parquet: N columns use a pandas sparse dtype``
   Junction count matrices load sparse-backed and no Parquet engine accepts
   :class:`pandas.SparseDtype`. Write a text format, or pass ``--densify``
   to materialize every implicit zero first. Densifying a junction matrix
   can need far more memory than the sparse form.

``KeyError: Missing required field: annotation_extension``
   Gene and exon descriptions need an annotation code. Pass it
   explicitly (``annotation_extension="G026"``) or use ``create_rse``,
   which resolves a default for you.

``TypeError: stack_count_matrices() got an unexpected keyword 'genomic_unit'``
   Filter the bundle before calling stack:
   ``bundle.filter(genomic_unit="gene").stack_count_matrices()``.

``RangesError: Could not derive genomic ranges …``
   The rest of the message names the cause: the annotation (or, for
   junctions, the RR coordinate file) could not be retrieved, could not be
   parsed, or does not cover every counted feature -- a mismatch, fixed by
   passing the matching ``annotation_extension``. A fourth variant, "no
   annotation providing genomic ranges was in the bundle", means there is
   no GTF or RR file to work from at all; include one at discovery time.
   ``allow_fallback_to_se=True`` returns a range-less
   :class:`~summarizedexperiment.SummarizedExperiment` instead of raising;
   it does not retry or repair anything.

``CompatibilityError: Incompatible count families …``
   You tried to stack gene/exon counts together with junctions. Filter
   to one family first, or stack each family separately.


Where to go next
----------------

- :doc:`api`: full per-symbol reference for all public modules.
- :doc:`cli`: the ``recount3`` command-line tool, which mirrors this API
  as a discover -> manifest -> materialize workflow.
- The `recount3 raw-files documentation
  <https://rna.recount.bio/docs/raw-files.html>`_ describes the underlying file
  layout (URLs, sharding, annotation codes). Note that this upstream page (not
  this tutorial) contains several inaccuracies.
