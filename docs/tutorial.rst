Tutorial
========

This tutorial walks through the ``recount3`` Python API end-to-end: finding
projects and samples, downloading the files behind them, assembling
:class:`~summarizedexperiment.SummarizedExperiment` and
:class:`~summarizedexperiment.RangedSummarizedExperiment` objects, merging
sample metadata, normalizing and scaling counts, reading BigWig coverage, and
managing the on-disk cache.

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
   python3 -m pip install "recount3[pybiocfilecache]"  # R cache sharing
   python3 -m pip install "recount3[all]"        # every optional feature

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
- ``pybiocfilecache`` enables a shared R/Python registry. The faster default
  filesystem cache needs no extra; see :ref:`cache-and-configuration`.
- ``all`` installs every optional feature above. The ``dev`` and ``docs``
  extras hold the test and documentation toolchains and are installed
  separately, so ``all`` does not pull them in.

If an optional dependency is missing, the affected function raises
:exc:`ImportError` on first use; the remainder of the package stays importable
and functional. The two output extras are additionally checked up front by the
CLI, before any download runs, so an unusable output format is reported
immediately rather than after the data has been fetched and assembled.


Quick start
-----------

.. note::

   Run each section's blocks in order, since later blocks reuse names bound by
   earlier ones. Examples that download data require network access unless
   their files are already cached. Annotation lookup, filtering, and analysis
   of loaded objects run offline. Downloaded files are cached under
   ``~/.cache/recount3/files`` (see :ref:`cache-and-configuration`), so
   re-running an example reuses the local copy rather than downloading again.

The most direct path from a project identifier to an analysis-ready BiocPy
object is :func:`~recount3.create_rse`. It requires the ``biocpy`` extra and
performs discovery, downloads, metadata merging, and range assembly in a single
call:

.. code:: python

   import recount3 as r3

   rse = r3.create_rse(
       project="SRP009615",
       organism="human",
       annotation_label="gencode_v26",
   )

   print("Features x samples:", rse.shape)
   print("Gene assay:", type(rse.get_assay("raw_counts")).__name__)
   print("First samples:", list(rse.get_column_names()[:3]))

Output::

   Features x samples: (63856, 12)
   Gene assay: ndarray
   First samples: ['SRR389077', 'SRR387777', 'SRR387778']

63,856 gene features by 12 samples. The counts are a plain
:class:`numpy.ndarray`; the BiocPy container supplies the feature and sample
labels around it. Exact counts and identifiers depend on the study and on what
the mirror currently serves.

Examples here import the package under the short alias ``r3``, so each call
shows where it comes from. ``from recount3 import create_rse`` works the same
way if you prefer it.

This single call is sufficient for the most common workflow; it is expanded in
:ref:`layer-1` below. The remainder of this tutorial describes the steps that
``create_rse`` performs internally and the lower-level components to use when
finer control is required.


Finding projects and samples before choosing a study
----------------------------------------------------

The quick start above assumes you already know which study you want. If you
do not, start from the source-level metadata: these calls read and parse the
mirror's own project and sample tables, so what they return is a record of
what exists rather than a guess at a filename. Limit ``data_sources`` to the
source you need. The human SRA table alone covers several thousand
projects:

.. code:: python

   import recount3 as r3

   projects = r3.available_projects(organism="human", data_sources="sra")
   print(projects[["project", "n_samples"]].head())

   samples = r3.available_samples(organism="human", data_sources="sra")
   study_samples = samples.loc[samples["project"].eq("SRP009615")]
   print(study_samples[["project", "external_id"]].head())

Both return a :class:`pandas.DataFrame`. The human SRA project table has 8,677
rows; ``n_samples`` there matches the column count you will get from
``create_rse``::

        project  n_samples file_source      project_home
   831  SRP009615         12         sra  data_sources/sra

These tables support selecting projects before downloading counts.
``r3.samples_for_project`` returns one project's sample identifiers as a list;
``r3.project_homes`` returns a table of project locations; and
``r3.create_sample_project_lists(organism="human")`` returns a
``(samples, projects)`` pair of sorted identifier lists for a whole organism,
which is what the CLI's ``recount3 ids`` writes out.
``r3.annotation_label("human", "G026")`` converts an extension back to its
human-readable label.

Keep the distinction in mind while reading the rest of this tutorial: the
discovery calls on this page read metadata, whereas the resource and bundle
layers below mostly *construct* candidate URLs from the parameters you give
them. A constructed URL is not evidence that the file exists.


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

   import recount3 as r3

   rse = r3.create_rse(
       project="SRP009615",
       organism="human",
       annotation_label="gencode_v26",   # or "gencode_v29", "fantom6_cat", "refseq", "ercc", "sirv"
   )

You may pass the raw extension code instead of a label:

.. code:: python

   rse = r3.create_rse(
       project="SRP009615",
       organism="human",
       annotation_extension="G026",
   )

When both are supplied, ``annotation_extension`` takes precedence. Discover the
available labels with :func:`~recount3.annotation_options`:

.. code:: python

   import recount3 as r3

   r3.annotation_options("human")
   r3.annotation_options("mouse")

Output::

   {'gencode_v26': 'G026', 'gencode_v29': 'G029', 'fantom6_cat': 'F006',
    'refseq': 'R109', 'ercc': 'ERCC', 'sirv': 'SIRV'}
   {'gencode_v23': 'M023'}

Exon-level and junction-level
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   exon_rse = r3.create_rse(
       project="SRP009615",
       organism="human",
       genomic_unit="exon",
       annotation_label="gencode_v26",
   )

   junction_rse = r3.create_rse(
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
   than the counts, for example a GENCODE 26 GTF against GENCODE 29 counts.
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

   experiment = r3.create_rse(
       project="SRP009615",
       organism="human",
       allow_fallback_to_se=True,
   )

Be deliberate about that flag:

- A fallback object **has no genomic ranges**. Counts, sample metadata,
  and the operations driven by column data still work, namely
  :func:`recount3.se.compute_scale_factors`,
  :func:`recount3.se.expand_sra_attributes`, and
  :func:`recount3.se.is_paired_end`. The helpers that need an RSE,
  ``compute_read_counts``, ``transform_counts``, and ``compute_tpm``, raise
  ``TypeError`` on a plain SE, as do range queries. The flag returns an RSE
  as usual whenever ranges are available.
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
  describes candidate files for a project, their URLs and annotation codes,
  as plain objects, so you can decide what is worth fetching, or write out a
  manifest for someone else to fetch.
- *Only part of what discovery returns.* You need gene counts and the QC
  table, without the other resources in a default bundle. ``create_rse``
  already limits counts to the requested genomic unit and annotation.
- *Your own assembly.* You want the stacked counts as a
  :class:`pandas.DataFrame` to merge with clinical data of your own, and
  will build the BiocPy object yourself (or skip it entirely).

A bundle is just a list of resources plus filters, so nothing is downloaded
by default during discovery without BigWigs. BigWig discovery reads a sample
index. Candidate URLs do not verify that files exist and do not report their
sizes. That is the practical difference from Layer 1:
``create_rse`` decides what to fetch for you, a bundle lets you decide.

Discovering resources for one or more projects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   import recount3 as r3

   bundle = r3.R3ResourceBundle.discover(
       organism="human",
       data_source="sra",
       project="SRP009615",
   )
   print(f"Found {len(bundle.resources)} resources.")

Output::

   Found 10 resources.

By default, ``discover`` returns gene + exon counts, the default annotation
GTF for each of those units (gene and exon), the five metadata tables, and the
default junction artifact (``MM``). For ``SRP009615`` this is the ten resources
counted above: 2 counts + 2 annotation GTFs + 5 metadata tables + 1 junction
file. Override with:

.. code:: python

   custom_bundle = r3.R3ResourceBundle.discover(
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
``project``. ``discover`` uses all combinations of the supplied values and
produces a single combined bundle:

.. code:: python

   multi = r3.R3ResourceBundle.discover(
       organism="human",
       data_source="sra",
       project=["SRP009615", "SRP001558"],
       genomic_units=("gene",),
   )
   print(f"Combined: {len(multi.resources)} resources, 2 projects.")

Output::

   Combined: 15 resources, 2 projects.

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
       annotation_extension=lambda ext: bool(ext) and ext.startswith("G"),
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
       .filter(annotation_extension="G026")
       .stack_count_matrices(compat="feature")
   )
   print(type(gene_counts_df).__name__, gene_counts_df.shape)

Output::

   DataFrame (63856, 12)

A plain :class:`pandas.DataFrame`, feature IDs on the index and sample IDs on
the columns, identical to what ``create_rse`` puts in its ``raw_counts``
assay. Junctions stack the same way, and stay sparse-backed:

.. code:: python

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

``compat="feature"`` does not verify an annotation build or align metadata.
Select one annotation explicitly before stacking. Use
``axis=1, join_policy="inner"`` to combine samples on shared feature IDs, and
check that the resulting sample names are unique. For experiment construction,
use the builders below; they also validate annotations and sample metadata.
For multi-project junctions, include MM, ID, and RR files and use the builders
to align junctions by genomic coordinates rather than project-local row numbers.

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

   import recount3 as r3

   desc = r3.R3ResourceDescription(
       resource_type="count_files_gene_or_exon",
       organism="human",
       data_source="sra",
       genomic_unit="gene",
       project="SRP009615",
       annotation_extension="G026",   # required for gene/exon counts
   )

   res = r3.R3Resource(desc)
   print(res.url)                     # fully-qualified URL on the recount3 mirror
   # http://duffel.rail.bio/recount3/human/data_sources/sra/gene_sums/15/SRP009615/sra.gene_sums.SRP009615.G026.gz

   res.download(path=None, cache_mode="enable")  # cache only, no local copy
   df = res.load()                               # pandas.DataFrame
   local_path = res.ensure_cached()              # path for another file reader
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
   res.download(path="./downloads")              # copy into a directory
   res.download(path="./recount3.zip")           # append to a ZIP archive

``path`` accepts a string or a :class:`pathlib.Path`
(:data:`~recount3.types.StrPath`), so a directory built with ``/`` works as
well, as does a ``Path`` handed back by
:meth:`~recount3.R3Resource.ensure_cached` or
:func:`~recount3.recount3_cache`:

.. code:: python

   from pathlib import Path

   out = Path("results-output")
   res.download(path=out / "downloads")

The same applies to ``dest`` on
:meth:`~recount3.R3ResourceBundle.download`.

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

   import recount3 as r3
   counts = r3.search_count_files_gene_or_exon(
       organism="human",
       data_source="sra",
       genomic_unit="gene",
       project="SRP009615",
       annotation_extension="G026",
   )

   meta = r3.search_metadata_files(
       organism="human",
       data_source="sra",
       project="SRP009615",
       table_name=(
           "recount_project", "recount_qc", "recount_seq_qc",
           "recount_pred", "sra",
       ),
   )

   bigwigs = r3.search_bigwig_files(
       organism="human",
       data_source="sra",
       project="SRP009615",
       sample=["SRR387777", "SRR387778"],
   )

The single-call equivalent is :func:`~recount3.search_project_all`
(used internally by ``R3ResourceBundle.discover``). The remaining helpers
follow the same shape:

==========================================  ==========================================
Function                                    Returns
==========================================  ==========================================
``r3.search_count_files_gene_or_exon``      Gene or exon count files
``r3.search_count_files_junctions``         Junction ``MM``/``ID``/``RR`` files
``r3.search_metadata_files``                Per-project metadata tables
``r3.search_bigwig_files``                  Per-sample BigWig coverage files
``r3.search_annotations``                   Annotation GTFs for an organism
``r3.search_data_sources``                  The ``homes_index`` for an organism
``r3.search_data_source_metadata``          Source-level (not project) metadata
``r3.search_project_all``                   Everything above for one project
==========================================  ==========================================

To go the other way, from a manifest line back to a resource, use
:meth:`~recount3.R3Resource.from_mapping`, which rehydrates one JSONL record
written by ``recount3 search``:

.. code:: python

   import json

   with open("manifest.jsonl", encoding="utf-8") as handle:
       resources = [r3.R3Resource.from_mapping(json.loads(line))
                    for line in handle if line.strip()]

The ``url`` and ``arcname`` keys in the record are recomputed from the
description and the active configuration, so a manifest written against one
mirror can be replayed against another.


Working with sample metadata
----------------------------

When ``create_rse`` or
:meth:`~recount3.R3ResourceBundle.to_ranged_summarized_experiment`
assembles an RSE, it merges all available per-project metadata tables
into ``column_data``, namespacing non-key columns by their table of
origin (e.g. ``recount_qc__star.all_mapped_reads``).

Construction uses ``metadata_join="inner"`` by default, matching R's
intersection of the nonempty metadata tables within each project. Samples
absent from that intersection are excluded from both the assay and
``column_data``. Use ``metadata_join="outer"`` on ``create_rse`` or either
bundle builder to retain every count sample with missing metadata values.
This option is independent of ``join_policy``, which joins feature rows
across projects. With no metadata resources, all count samples are retained.

Requested files that fail to load, conflicting sample identifiers, and mixed
annotations raise errors. Select ``annotation_extension`` explicitly when a
bundle contains several annotations. For junctions from multiple projects,
each MM matrix needs its matching ID and RR sidecars: junction rows are
matched by chromosome, inclusive coordinates, and strand. An outer feature
join introduces zeros only for features absent from a project; missing values
inside an input count file are errors.

Gene and exon assays preserve their numeric types. Junction assays remain
SciPy CSC sparse matrices, and count-transform helpers return sparse-backed
DataFrames for sparse inputs. Standard BiocPy assay access, slicing, copying,
and range operations work directly on the constructed objects. Converting a
large junction assay with ``toarray()`` explicitly allocates its dense form.

Repeated exon IDs retain their individual transcript annotations in
``row_data`` and ``row_ranges``. Python gives repeated rows unique names while
preserving the original IDs in ``row_data["feature_id"]``. Compatible projects
with identical repeated feature ordering can be combined; ambiguous repeated
feature alignments raise an error. The GTF score is preserved as ``bp_length``
(covered exonic length), which can differ from genomic span. TPM uses this
annotated length when present.

Experiment metadata records the project, organism, annotation, source URLs,
construction options, creation time, package version, and the mapping of
metadata columns to their source tables. Access individual fields through
``rse.metadata["project"]``. A bundle retains at most one aligned annotation
cache for repeated RSE construction; changing its annotation file invalidates
the cache. With ``autoload=False``, counts and sample metadata must already be
loaded, while genomic range files must be cached locally.

Access it as a pandas DataFrame:

.. code:: python

   col_df = rse.get_column_data().to_pandas()
   col_df.columns[:10]

Expanding SRA sample attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In an assembled RSE, SRA samples carry an ``sra__sample_attributes`` column
(the ``sra`` metadata table namespaced with ``__`` as described above) that
encodes key-value pairs in the form ``"age;;67.78|disease;;Control|..."``.
:func:`recount3.se.expand_sra_attributes` parses these into separate columns.
It accepts either a DataFrame or an SE/RSE object, and recognizes both the
namespaced ``sra__sample_attributes`` name and the R-style
``sra.sample_attributes`` spelling. Each parsed attribute becomes a new column
named ``sra_attribute.<key>`` (for example, ``sra_attribute.disease``):

.. code:: python

   import recount3 as r3

   rse2 = r3.se.expand_sra_attributes(rse)
   col_df = rse2.get_column_data().to_pandas()
   sra_cols = [c for c in col_df.columns if c.startswith("sra_attribute.")]
   print(sra_cols)
   print(col_df[["sra_attribute.cells", "sra_attribute.cell_line"]].head(4))

Output::

   ['sra_attribute.cell_line', 'sra_attribute.shRNA_expression',
    'sra_attribute.source_name', 'sra_attribute.treatment',
    'sra_attribute.cells']
             sra_attribute.cells sra_attribute.cell_line
   SRR389077                 NaN                    K562
   SRR387777                K562                     NaN
   SRR387778                K562                     NaN
   SRR389078                 NaN                    K562

Worth dwelling on: the submitters of this study recorded the same fact under
two different attribute names, so neither column alone describes every sample.
``expand_sra_attributes`` reports what the submitters wrote and does not
reconcile it for you. Inspect the expanded columns before defining analysis
groups; here you would combine them, for example with
``col_df["sra_attribute.cell_line"].combine_first(col_df["sra_attribute.cells"])``.


Normalization and scaling
-------------------------

Gene and exon matrices contain base-pair coverage sums (the ``raw_counts``
assay). Junction matrices contain junction-supporting counts (the ``counts``
assay); do not apply the gene coverage-to-read or TPM formulas to junctions.
:mod:`recount3.se` provides helpers consistent with the R implementation.
These reach the ``se`` submodule rather than the package root, so they are
called as ``r3.se.compute_tpm(rse)``, whereas the builders and discovery
helpers used so far are available directly as ``r3.create_rse(...)``.

``compute_read_counts``, ``transform_counts``, and ``compute_tpm`` require an
RSE. ``compute_scale_factors`` and ``is_paired_end`` also accept sample
metadata DataFrames and plain SE objects.

Approximate read counts
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   import recount3 as r3

   reads = r3.se.compute_read_counts(rse)        # pandas DataFrame, integer-rounded

Values are rounded to whole reads by default; pass ``round_to_integers=False``
to retain the fractional estimates.

Per-sample scale factors
~~~~~~~~~~~~~~~~~~~~~~~~

Two methods are supported, matching the R ``recount3`` reference:

.. code:: python

   import recount3 as r3

   sf_auc = r3.se.compute_scale_factors(rse, by="auc")
   sf_mapreads = r3.se.compute_scale_factors(rse, by="mapped_reads")
   print(pd.DataFrame({"auc": sf_auc, "mapped_reads": sf_mapreads}).head(3))

Both return a :class:`pandas.Series` indexed by ``external_id``::

                     auc  mapped_reads
   external_id
   SRR389077    0.045855      0.129044
   SRR387777    0.039961      0.112357
   SRR387778    0.034571      0.097252

Apply scale factors to the assay:

.. code:: python

   scaled = r3.se.transform_counts(rse, by="auc")          # default
   scaled = r3.se.transform_counts(rse, by="mapped_reads", target_read_count=4e7)

TPM (gene/exon only, needs feature widths)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   import recount3 as r3

   tpm = r3.se.compute_tpm(rse)
   print(type(tpm).__name__, tpm.shape)
   print(tpm.sum(axis=0).round().unique())

Output::

   DataFrame (63856, 12)
   [1000000.]

Every sample sums to one million across the full feature set, which is the
quickest check that the normalization ran over the whole matrix rather than a
subset. The five highest-expressed genes in the first three samples::

                          gene_name  SRR389077  SRR387777  SRR387778
   ENSG00000210082.2        MT-RNR2   36013.14   16179.52   14534.20
   ENSG00000281383.1  CH507-513H4.5   67606.31    7793.08   10009.38
   ENSG00000213934.6           HBG1    7827.91   12294.58   10770.67
   ENSG00000198712.1         MT-CO2    6718.71    9145.23    6520.00
   ENSG00000228253.1        MT-ATP8    4585.29    6566.27    8629.65

TPM uses annotated ``bp_length`` (covered exonic bases) when present, falling
back to range widths. A gene's genomic span can include introns and is not an
interchangeable length. Normalize the full gene matrix before selecting genes
for display; each nonzero sample should sum to approximately one million.
Check for missing or non-finite results before downstream analysis. The helpers
return new DataFrames and do not replace the RSE assay; rounded estimates have
integer-like values but may retain a floating dtype.

:func:`recount3.se.is_paired_end`, ``compute_scale_factors``, and
``expand_sra_attributes`` accept a metadata DataFrame or an SE/RSE object.
See the API reference for full signatures.


From package objects to downstream analysis
-------------------------------------------

The gene assay is a NumPy array, metadata and normalized counts are pandas
objects, and the junction assay is a SciPy sparse matrix. Inspect the real
objects and preserve their labels when moving between representations:

.. code:: python

   import numpy as np
   import pandas as pd
   from scipy import sparse

   raw = rse.get_assay("raw_counts")
   print(type(raw), raw.shape)
   print(type(tpm), tpm.shape)
   assert np.isfinite(tpm.to_numpy()).all()
   nonzero = tpm.sum(axis=0) > 0
   np.testing.assert_allclose(tpm.sum(axis=0)[nonzero], 1_000_000)

   # Filter only after normalizing the complete gene matrix.
   expressed = (tpm >= 1).sum(axis=1) >= 3
   log_tpm = np.log2(tpm.loc[expressed] + 1)
   variable_ids = log_tpm.var(axis=1).nlargest(2000).index
   X = log_tpm.loc[variable_ids].T.to_numpy()  # samples × genes
   sample_metadata = rse.get_column_data().to_pandas().reindex(tpm.columns)
   assert sample_metadata.index.tolist() == tpm.columns.tolist()
   correlation = pd.DataFrame(
       np.corrcoef(X), index=tpm.columns, columns=tpm.columns,
   )
   print("Analysis array (samples x genes):", X.shape, type(X).__name__)
   print(correlation.iloc[:4, :4].round(3))

Output::

   Analysis array (samples x genes): (12, 2000) ndarray
              SRR389077  SRR387777  SRR387778  SRR389078
   SRR389077      1.000      0.731      0.736      0.902
   SRR387777      0.731      1.000      0.955      0.723
   SRR387778      0.736      0.955      1.000      0.724
   SRR389078      0.902      0.723      0.724      1.000

``X`` is a samples-by-genes :class:`numpy.ndarray`, the orientation
scikit-learn and most Python machine-learning tooling expect, and it can be
handed straight to an estimator.

The same array summarizes by principal component analysis (PCA). NumPy's SVD
is enough; no extra dependency is needed:

.. code:: python

   X_centered = X - X.mean(axis=0, keepdims=True)
   U, singular_values, _ = np.linalg.svd(X_centered, full_matrices=False)
   scores = U[:, :2] * singular_values[:2]
   variance_fraction = singular_values**2 / np.sum(singular_values**2)
   pc_df = pd.DataFrame(scores, index=tpm.columns, columns=["PC1", "PC2"])
   print(pc_df.head(4).round(3))
   print("PC1/PC2 explained variance (%):",
         np.round(100 * variance_fraction[:2], 2).tolist())

Output::

                 PC1     PC2
   SRR389077 -56.396   5.020
   SRR387777  -6.816 -23.199
   SRR387778  -4.931 -20.632
   SRR389078 -67.903   7.778
   PC1/PC2 explained variance (%): [42.19, 14.57]

This is an exploratory expression-profile comparison. The thresholds are
explicit choices for SRP009615, not universal defaults; component signs can
flip between numerical libraries without changing the result. For prediction,
fit gene selection and centering within training folds. This example does not
estimate treatment effects or remove batch effects.

For a junction RSE built earlier, summarize sparse counts without allocating
the entire dense matrix:

.. code:: python

   junctions = junction_rse.get_assay("counts")
   assert sparse.issparse(junctions)
   totals = np.asarray(junctions.sum(axis=0)).ravel()
   detected = np.asarray((junctions > 0).sum(axis=0)).ravel()
   junction_summary = pd.DataFrame(
       {"count_sum": totals, "detected_junctions": detected},
       index=junction_rse.get_column_names(),
   )
   print(type(junctions).__module__ + "." + type(junctions).__name__)
   print("Features x samples:", junctions.shape, "stored entries:", junctions.nnz)
   print(junction_summary.head(3))

Output::

   scipy.sparse._csc.csc_matrix
   Features x samples: (281448, 12) stored entries: 1341130
              count_sum  detected_junctions
   SRR389079    1732182              142791
   SRR389080    1315344              117890
   SRR389081     844075              106754

281,448 junctions by 12 samples, of which 1,341,130 cells are nonzero. The
three SciPy buffers hold about 16 MB against roughly 27 MB for the equivalent
dense array, and that gap widens sharply for a multi-project assembly.

The default junction assay is named ``counts``, not ``raw_counts``. Summing or
filtering it with SciPy keeps the full matrix sparse; ``toarray()`` allocates
every cell. For Parquet output, pandas sparse columns must be densified
explicitly (CLI ``--densify``), so estimate the memory requirement first.

Exporting to AnnData
~~~~~~~~~~~~~~~~~~~~

AnnData is the container most Python single-cell and machine-learning tooling
reads. It is transposed relative to BiocPy: samples are rows (``obs``),
features are columns (``var``), and each recount3 assay becomes a layer.
:func:`recount3.se.to_anndata` performs the conversion and needs the
``anndata`` extra:

.. code:: python

   adata = r3.se.to_anndata(rse)
   print("samples x genes:", adata.shape)
   print("uns provenance:", adata.uns["project"], adata.uns["annotation"],
         len(adata.uns["resource_urls"]), "resource URLs")
   print("layers:", list(adata.layers), "| obs columns:", adata.obs.shape[1])

Output::

   samples x genes: (12, 63856)
   uns provenance: SRP009615 gencode_v26 7 resource URLs
   layers: ['raw_counts'] | obs columns: 176

Sample metadata travels as ``obs``, gene annotations as ``var``, and the
experiment's provenance as ``uns``.

Call :func:`~recount3.se.to_anndata` rather than the BiocPy
``rse.to_anndata()`` method. BiocPy holds experiment metadata as a
``NamedList`` and hands it to ``AnnData(uns=...)`` unchanged, which AnnData
rejects with ``Only mutable mapping types (e.g. dict) are allowed for
`.uns`.``. Since every experiment built here carries provenance metadata, that
applies to all of them. :func:`~recount3.se.to_anndata` converts the
provenance to plain dicts and lists first, so ``project``, ``annotation``,
``resource_urls``, and ``metadata_columns`` all survive the conversion.

Writing the result to HDF5 has a second requirement. recount3 STAR QC fields
are named after splice motifs and contain ``/``, which HDF5 reads as a path
separator, both in the sample column names and in the ``uns`` provenance map
keyed by them. Pass ``sanitize_for_hdf5=True`` to rename them:

.. code:: python

   adata = r3.se.to_anndata(rse, sanitize_for_hdf5=True)
   adata.write_h5ad("rse.h5ad")

That is the same fixup ``recount3 bundle rse --sanitize-columns`` applies.
Renaming is opt-in because it changes what an analysis indexes by; each
renamed provenance entry keeps its original name in its value, so nothing is
lost. If you would rather not rename anything, write ``.pkl`` instead, which
preserves the RSE itself, including its genomic ranges, which AnnData has no
place for.


BigWig coverage
---------------

Per-sample BigWig coverage files are not included by default; pass
``include_bigwig=True`` (or use the ``search_bigwig_files`` helper)
to add them. Requires the ``bigwig`` extra.

.. code:: python

   bundle = r3.R3ResourceBundle.discover(
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
   bw = bw_res.load()                       # a BigWigFile wrapper
   with bw:                                 # closes the handle on exit
       values = bw.values("chr1", 0, 1000, numpy=True)

Note which object each spelling gives you. ``load()`` returns the
:class:`~recount3._bigwig.BigWigFile` wrapper, so ``bw`` above is the wrapper
and reopens its handle automatically on the next read after a ``close()``.
Entering the wrapper as a context manager instead yields the live ``pyBigWig``
handle, not the wrapper:

.. code:: python

   bw_res = r3.search_bigwig_files(
       organism="human", data_source="sra", project="SRP009615",
       sample="SRR387777",
   )[0]
   with bw_res.load() as bw:                # bw is the pyBigWig handle here
       coverage = bw.values("chr1", 10_000, 11_000, numpy=True)
       mean_signal = bw.stats("chr1", 10_000, 11_000, type="mean", exact=True)[0]
   print(type(coverage).__name__, coverage.shape, "| mean:", mean_signal)

Output::

   ndarray (1000,) | mean: 0.036

Per-base coverage comes back as a :class:`numpy.ndarray`, one value per base
over the requested interval.

Both spellings close the file on exit and expose the same ``values()``,
``stats()``, and ``intervals()`` calls, so either spelling works; just do not
expect wrapper-only behaviour from the second.

BigWig intervals use zero-based, half-open coordinates; uncovered positions can
be NaN. ``load()`` first caches the full BigWig file. Reading a small interval
limits the array returned, not the initial download. To limit downloads to one
sample, use ``r3.search_bigwig_files(..., sample="SRR387777")`` as shown above
rather than discovering a whole project's coverage files.


.. _cache-and-configuration:

Cache and configuration
-----------------------

The default filesystem backend stores downloads under
``~/.cache/recount3/files``. The optional ``pybiocfilecache`` backend defaults
to R's recount3 cache directory, as described below.
The :mod:`recount3.config` helpers let you inspect and prune the cache:

.. code:: python

   import recount3 as r3
   print(r3.recount3_cache())                 # cache directory Path
   files = r3.recount3_cache_files(pattern="*.gtf.gz")

   # Dry-run a deletion first:
   to_remove = r3.recount3_cache_rm(
       predicate=lambda p: ".junctions." in p.name,
       dry_run=True,
   )
   r3.recount3_cache_rm(predicate=lambda p: ".junctions." in p.name)


Threaded operations
~~~~~~~~~~~~~~~~~~~

``R3ResourceBundle.download(max_workers=8)`` and ``recount3 download --jobs 8``
use worker threads. A single ``R3Resource.download()`` performs one operation.
``create_rse()`` and the bundle's experiment constructors do not automatically
call the parallel bundle downloader. Prefetch a bundle explicitly before
constructing an experiment when overlapping transfers is useful::

    bundle.download(dest="downloads", cache="enable", max_workers=4)
    rse = bundle.to_ranged_summarized_experiment(genomic_unit="gene")

The worker count bounds concurrency, not guaranteed throughput. Network
capacity, server policies, dataset count, disk I/O, and ZIP writes can limit
scaling. Experiment construction and parsing have separate memory costs.

Cache destinations and refreshes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cache transfers are coordinated by the canonical destination path, resolving
relative paths and symlink aliases. Different cache entries may transfer
concurrently. Requests for the same missing entry check existence under its
lock: after one succeeds, waiting ``enable`` requests reuse that file. If a
transfer fails, its lock is released and a waiter can attempt the download.
Locks are weakly retained, so an unbounded list of historical URLs is not kept.

Each ``update`` request performs its own transfer, including simultaneous
requests. Updates to the same destination run sequentially; lock acquisition
order is unspecified. An ``enable`` request waits for a transfer already
holding that destination's lock, then uses the available file. It does not
force a further refresh. The last successful update determines the cache
contents. Returned paths are not immutable snapshots of that version.

Transfers write a temporary sibling and atomically replace the payload only
on success. A failed transfer leaves any previous payload intact and removes
its temporary file. Already-open readers and hard-linked materializations can
continue to refer to the previous version after replacement. A later refresh
does not silently update those materializations or parsed objects in memory.

ZIP writes retain a lock per canonical archive path. Downloads for different
archive members can overlap, but adding or replacing members is serialized.
With ``cache="disable"``, ZIP downloads use temporary files before inserting
complete members; directory downloads stream to their destinations without
cache deduplication. Cache-disabled calls do not load pybiocfilecache.

All these locks coordinate threads **within one Python process**. They do not
coordinate separate CLI jobs, Python processes, R sessions, or cluster nodes.
Use separate cache and output paths per concurrent process, or provide external
coordination. Run cache removal, registry maintenance, and external file edits
only when the cache is idle. Atomic payload replacement alone does not provide
cross-process deduplication or transactional registry and payload updates.


Optional shared R/Python cache
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

recount3 already includes a fast URL-hash filesystem cache. It is the default
and requires no cache extra. Use the optional ``pybiocfilecache`` backend when
you need to share cached resources with R's BiocFileCache. Its SQLite lookups,
metadata writes, and local-file checksums add measurable overhead; it is not
a performance upgrade. Installing the extra alone does not change the default.

In a controlled WSL2 benchmark (Python 3.14.7, three repetitions, four 1 MiB
files, 25 ms per loopback HTTP request), median download times were:

.. list-table:: Native cache versus optional shared cache
   :header-rows: 1
   :widths: 30 15 20 20

   * - Cache state / workers
     - Native (s)
     - Shared (s)
     - Shared/native
   * - Empty / 1
     - 0.131
     - 0.429
     - 3.3x
   * - Empty / 4
     - 0.054
     - 0.365
     - 6.8x
   * - Warm / 1
     - 0.0015
     - 0.0382
     - 26.1x
   * - Warm / 4
     - 0.0091
     - 0.0622
     - 6.9x

The optional backend added roughly 22--25 MiB to peak client RSS in these
cases. Empty-cache timing includes the optional import and registration;
warm timing excludes prefill but includes lookups and local-file checksums.
These figures describe these small controlled workloads, not public-mirror
throughput or a universal slowdown ratio. Transfer size, network latency,
registry size, and filesystem performance affect the overhead. The native
backend is the faster default; shared-cache reuse can still avoid a much more
expensive download when the file was already obtained from R.

Install and explicitly select it with::

    python -m pip install "recount3[pybiocfilecache]"
    export RECOUNT3_CACHE_BACKEND=pybiocfilecache
    recount3 download --from resources.jsonl --dest downloads

Or select it using the global CLI flag, before the subcommand::

    recount3 --cache-backend pybiocfilecache download \
        --from resources.jsonl --dest downloads

Both examples use R's default recount3 cache directory. No R installation
is needed to determine the path. Python follows the rules used by
``tools::R_user_dir("recount3", "cache")``:

- A nonempty ``R_USER_CACHE_DIR`` takes precedence over ``XDG_CACHE_HOME``.
  With either variable, the directory is ``<value>/R/recount3``.
- Otherwise, Linux/WSL uses ``~/.cache/R/recount3``; macOS uses
  ``~/Library/Caches/org.R-project.R/R/recount3``; Windows uses
  ``%LOCALAPPDATA%/R/cache/R/recount3``.

An explicit ``--cache-dir`` overrides ``RECOUNT3_CACHE_DIR``, which overrides
these defaults. The native ``filesystem`` backend keeps its existing
``~/.cache/recount3/files`` default. Existing files are not moved when the
backend changes. If R uses a custom ``options(recount3_cache=...)`` setting,
set the same directory explicitly in Python; Python cannot read R session
options.

In Python, select the backend before its default directory is resolved::

    import recount3 as r3

    cfg = r3.default_config(cache_backend="pybiocfilecache")

To use a custom shared directory instead::

    from dataclasses import replace
    from pathlib import Path
    import recount3 as r3

    cfg = replace(
        r3.default_config(cache_backend="pybiocfilecache"),
        cache_dir=Path("/data/recount3-shared"),
    )

Pass ``config=cfg`` to :class:`~recount3.R3Resource`, the only entry point
that retrieves data and takes this keyword. The cache helpers
:func:`~recount3.recount3_cache`, :func:`~recount3.recount3_cache_files`, and
:func:`~recount3.recount3_cache_rm` accept it too, so that cache inspection
and pruning address the same backend. The search functions,
``R3ResourceBundle.discover``, and ``create_rse`` do not take it; configure
their defaults through the environment variables below before calling them.
The supported backend values are ``filesystem`` and ``pybiocfilecache``.
``recount3[all]`` includes this extra.
The CLI flag overrides the environment; the default remains ``filesystem``.
Changing an existing configuration with ``dataclasses.replace`` preserves
its ``cache_dir`` unless that field is also replaced.

The optional backend opens ``BiocFileCache.sqlite`` directly in ``cache_dir``.
With the default location, R's ``recount3::recount3_cache()`` uses the same
directory. For a custom location, point R at that **same directory**::

    bfc <- recount3::recount3_cache("/data/recount3-shared")
    info <- BiocFileCache::bfcinfo(bfc)
    BiocFileCache::bfcpath(bfc, rids = info$rid)

The supported pybiocfilecache 0.7.x series maintains R database compatibility.
Both interfaces use schema version ``0.99.4`` and the same resource columns.
Sequential interoperability was tested with the supplied R recount3 source
and BiocFileCache 3.0.0: R can query and use Python-created entries, Python can
reuse R-created entries, and resource IDs remain unique after interleaved
additions and deletion. SQL defaults can differ without preventing these reads.

The exact remote URL is the lookup key (``rname``). An existing R record is
used at its registered path, including relative and web-resource paths.
Python does not download another URL-hash copy when that registered file is
present. A missing registered payload is repaired at the same destination.
Multiple records with the same URL name are ambiguous and raise ``ValueError``;
resolve those duplicates with BiocFileCache before retrying.

New Python downloads retain the native URL-hash filenames and are registered
as local files without an extra copy. Existing native files are adopted on
cached download. R can find them by URL through ``bfcquery`` or ``bfcrpath``.
R's ``recount3::file_retrieve`` also makes a remote availability check before
consulting its cache; use ``bfcpath`` to access cached data without that check.

recount3 continues to handle network timeouts, TLS, retries, and atomic
replacement. It does not delegate transfers to the optional package's web
downloader. Existing R resource types and relative paths are preserved.
A Python refresh of an R web record clears the old HTTP validators, because
the unconditional transfer does not collect replacements. Local records use
pybiocfilecache's checksum; even warm local-cache downloads recompute it.

Registry operations are serialized per shared database within this Python
process. Network transfers occur outside that lock. A registry error propagates
while leaving any completed payload available. Retrying a local cache hit
repairs its checksum; after an interrupted web refresh, repeat ``update`` to
reset its validators. Payload replacement and SQLite commits are not one
transaction, so a process crash may require registry maintenance.

``recount3_cache_files()`` excludes SQLite databases, journals, and R lock
files. With this backend it includes registered paths outside the cache root and
registered entries whose payload was removed externally. Patterns filter these
paths. ``recount3_cache_rm()`` deletes selected payloads inside the cache root
and removes their registry rows, including stale rows for missing payloads.
For registered files outside the cache root, including R's ``action="asis"``
references, it removes only the registry rows and leaves the files untouched.
Containment is checked after resolving symlinks. The returned list includes
external paths whose registrations were removed; it does not mean those files
were deleted. The database itself is preserved, and ``dry_run=True`` changes
neither files nor registry rows. Removal can still affect cached files used
by R and other applications; restrict it with a predicate as appropriate.
Run maintenance while all users of the shared cache are idle.

Switching to the native backend avoids SQLite lookups and checksum overhead.
It does not consult R's registry or maintain its rows. Use the optional backend
for removal of registered files. Absolute-path local entries require repair if
the cache is moved; relative R entries remain relative to the shared root.

Inspect the same registry directly from Python with::

    from pybiocfilecache import BiocFileCache

    with BiocFileCache(cfg.cache_dir) as registry:
        records = registry.list_resources()

Shared format does not imply simultaneous-transfer coordination. The thread
locks in recount3 do not coordinate R sessions or other Python processes.
Use the shared directory sequentially, or provide external coordination for
writers, refreshes, and removal.


Configuration precedence is, from lowest to highest: library defaults,
environment variables, an explicit :class:`~recount3.Config` passed to a
resource or search function. The supported environment variables are:

==============================  ========================================
Variable                        Effect
==============================  ========================================
``RECOUNT3_URL``                Base URL of the recount3 mirror
``RECOUNT3_CACHE_DIR``          On-disk cache directory
``RECOUNT3_CACHE_BACKEND``      filesystem or pybiocfilecache
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
   import recount3 as r3

   cfg = r3.Config(
       base_url="http://duffel.rail.bio/recount3/",
       timeout=60,
       insecure_ssl=False,
       max_retries=5,
       user_agent="my-pipeline/0.1",
       cache_dir=Path("/scratch/recount3_cache"),
       cache_disabled=False,
       chunk_size=1024 * 1024,
   )

   res = r3.R3Resource(
       r3.R3GeneOrExonCounts(
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
   AnnData export needs ``anndata`` and ``delayedarray``. Run
   ``pip install "recount3[anndata]"``, or write a ``.pkl`` instead.

``Cannot write .h5ad: N name(s) contain a forward slash``
   HDF5 reads ``/`` as a path separator, and recount3 STAR QC fields are
   named after splice motifs (``..._gt/ag``). This affects the sample
   columns and the ``uns`` provenance map keyed by them. Pass
   ``--sanitize-columns`` to rename both to ``..._gt_ag`` (in Python,
   ``r3.se.to_anndata(rse, sanitize_for_hdf5=True)``), or write a ``.pkl``
   instead, which keeps the names verbatim.

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
   parsed, or does not cover every counted feature. This is a mismatch, fixed by
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
