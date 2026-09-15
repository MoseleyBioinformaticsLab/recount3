Changelog
=========

1.2.0 (2026-09-15)
------------------

Breaking changes
~~~~~~~~~~~~~~~~

- Experiment construction intersects nonempty metadata tables by default,
  matching R. Pass ``metadata_join="outer"`` to retain count samples with
  incomplete metadata, which is the previous behavior.
- Construction raises ``ValueError`` instead of returning an experiment with
  no samples when sample alignment discards every count sample. Callers that
  checked the result for emptiness must catch the error instead.
- Sparse junction inputs produce sparse (SciPy CSC) assays instead of being
  densified, and dense assays keep their own numeric dtype instead of being
  cast to ``float64``.
- Counts that are missing, non-finite, or non-numeric raise instead of
  becoming zeros; duplicate or blank sample identifiers raise instead of
  reaching the assay; and a requested count matrix that cannot be loaded
  raises instead of being left out of the assembled experiment.
- A bundle offering more than one annotation for the requested genomic unit
  raises ``CompatibilityError`` instead of choosing one of them by heuristic.
  Filter the bundle, or name the one to use with ``annotation_extension``.
- ``search_project_all()`` validates the project against the data source's
  sample index only when ``include_bigwig=True``. Callers that relied on
  discovery raising for an unknown project without requesting BigWig files
  must check the project themselves.
- The ``biocpy`` extra requires ``summarizedexperiment>=0.7.1``, up from
  ``>=0.6``.
- ``recount3[all]`` installs only the optional runtime features (``bigwig``,
  ``biocpy``, ``parquet``, ``anndata``, and ``pybiocfilecache``). It
  previously also pulled in the ``dev`` and ``docs`` toolchains, so ``pip
  install "recount3[all]"`` installed pytest, pylint, build, twine, and Sphinx
  alongside the features it advertised. Contributors now install those two
  extras explicitly.

Added
~~~~~

- ``recount3.se.to_anndata(experiment, sanitize_for_hdf5=False)`` converts a
  ``SummarizedExperiment`` to an ``anndata.AnnData``, with samples in rows,
  features in columns, each assay a layer, and provenance in ``uns``. BiocPy's
  own ``experiment.to_anndata()`` cannot be used on a recount3 experiment: it
  forwards metadata straight to ``AnnData(uns=...)``, and BiocPy stores that
  metadata as a ``NamedList``, which AnnData rejects. Every experiment the
  package builds carries provenance, so that method failed for all of them.
  ``sanitize_for_hdf5=True`` additionally makes the object writable by
  ``write_h5ad``.
- ``parquet`` and ``anndata`` extras, both included in ``all``, and pre-flight
  dependency checks for the outputs that need them. ``parquet`` installs
  ``pyarrow``; ``anndata`` installs ``anndata`` and ``delayedarray`` and
  implies ``biocpy``. ``bundle stack-counts`` resolves a Parquet engine, and
  ``bundle se``/``bundle rse`` verify anndata support, before anything is
  downloaded, so a missing dependency fails in a second rather than after the
  matrix has been assembled. (#8)
- An opt-in ``pybiocfilecache`` extra and a second cache backend that shares
  R's BiocFileCache registry, selected with ``cache_backend="pybiocfilecache"``
  in Python, the global ``--cache-backend`` flag, or the
  ``RECOUNT3_CACHE_BACKEND`` environment variable. It defaults to R's recount3
  cache directory, honoring ``R_USER_CACHE_DIR``, ``XDG_CACHE_HOME``, and the
  platform-specific locations, so a file an R session already downloaded is
  reused where it lies instead of being fetched again; an explicit cache
  directory still takes precedence. The default filesystem backend is
  unchanged, requires no extra, and remains the faster of the two.
- ``bundle stack-counts --densify`` materializes implicit zeros so that
  sparse-backed count matrices can be written to Parquet. Junction matrices
  load sparse and no Parquet engine accepts pandas sparse dtypes; without the
  flag, ``.parquet`` output now reports which columns are sparse and names
  the alternatives instead of failing inside the engine. Text formats, which
  already write dense rows, are unaffected.
- ``bundle se --sanitize-columns`` and ``bundle rse --sanitize-columns``
  replace ``/`` with ``_`` in ``obs``/``var`` column names before writing
  ``.h5ad``. HDF5 reads a forward slash as a path separator and recount3's
  STAR QC fields are named after splice motifs (for example,
  ``recount_qc__star.number_of_splices:_gt/ag``), so most projects carrying
  SRA metadata need it. Columns that are entirely missing are cast for HDF5
  automatically. ``.pkl`` output keeps every name verbatim.
- ``metadata_join`` on the SE/RSE builders and their wrappers. The default
  ``"inner"`` intersects the nonempty metadata tables; ``"outer"`` explicitly
  retains count samples with incomplete metadata.
- ``RangesError`` and its subclasses ``MissingRangesError`` and
  ``RangesCoverageError`` in ``recount3.errors``, also exported at the top
  level. They make the reason ranges could not be derived catchable rather
  than only readable: ``MissingRangesError`` means nothing in the bundle can
  supply ranges, ``RangesCoverageError`` that the file does not describe every
  counted feature. All three also subclass ``ValueError``, which is what
  ``to_ranged_summarized_experiment`` has always raised, so existing handlers
  keep working.
- ``R3Resource.ensure_cached(download=...)`` returns a resource's local path,
  fetching it first when absent. ``_cached_path()`` only computes a path and
  never reports a cache miss, so callers that read the cached file themselves
  had no way to prepare it.
- ``annotation_label()``, the inverse of ``annotation_ext()``, exported at the
  top level. It names an annotation the way users and R do
  (``"gencode_v26"``) given the file-naming code (``"G026"``). Unlike
  ``annotation_ext()``, an organism or extension this release does not
  recognize is returned unchanged rather than raising, because the label is
  only ever read by a human.
- Experiment provenance. Built experiments now carry the creation time,
  the ``recount3`` version, project, organism, project home, genomic unit,
  annotation name and extension, the URLs every contributing resource came
  from, the metadata columns and their source tables, both join policies,
  and the mirror base URL.

Changed
~~~~~~~

- Cached downloads are coordinated per canonical destination rather than by a
  single process-wide lock, so ``R3ResourceBundle.download(max_workers=...)``
  and ``recount3 download --jobs`` transfer different files at the same time.
  The previous module-level lock was held for the length of each transfer, so
  worker threads queued behind one another and the configured concurrency
  bought little. Requests for the same missing file still deduplicate, each
  forced refresh (``cache="update"``) takes its destination's lock in turn,
  and writes into one ZIP archive remain serialized. These locks coordinate
  threads within one Python process; separate CLI jobs, Python processes, R
  sessions, and cluster nodes are not coordinated.
- Count transforms return sparse-backed DataFrames for sparse inputs instead
  of densifying them.
- A metadata table that recount3 never published for a project is dropped
  with a warning as long as one usable table remains, and a table that was
  retrieved but cannot be parsed raises, being a damaged file rather than an
  absent one. Every kind of failure was previously skipped in silence. A
  bundle whose metadata tables were all empty or unretrievable now raises
  rather than producing an experiment with no metadata.
- ``compute_scale_factors(by="auc")``, and ``transform_counts`` through it, no
  longer resolves paired-end status, which cancels out of the AUC formula.
  Scaling by AUC therefore works on projects whose ``recount_seq_qc`` table is
  empty or was never published, where inferring paired-end status raised over
  the missing ``recount_seq_qc.avg_len`` column. ``by="mapped_reads"``, which
  does use paired-end status, is unchanged.
- MatrixMarket junction files are parsed into a ``scipy.sparse.csr_array``
  instead of a ``scipy.sparse.csr_matrix``. SciPy 1.18 deprecated the
  implicit sparse-matrix return of ``scipy.io.mmread`` ahead of flipping the
  default in SciPy 1.20; ``recount3`` now requests the new behavior
  explicitly where the installed SciPy supports it and normalizes the result,
  so the parsed type is identical across the whole supported SciPy range. The
  keyword is feature-detected, keeping SciPy releases older than 1.18 (the
  newest installable on Python 3.10 and 3.11) working unchanged.
- When genomic ranges cannot be derived, the reported reason distinguishes a
  failure to retrieve the annotation, an annotation that cannot be parsed,
  and an annotation that does not cover every counted feature. Only the last
  is fixed by choosing a different ``annotation_extension``.
- ``to_ranged_summarized_experiment`` raises ``RangesError`` instead of a bare
  ``ValueError``. ``RangesError`` is a ``ValueError``, so this is
  source-compatible.
- The "falling back" warning is emitted only when a plain
  ``SummarizedExperiment`` is actually returned, instead of whenever range
  derivation failed.
- ``autoload`` now reaches annotation selection, so ``autoload=False``
  inspects only already-cached annotations and never downloads.
- Annotation selection skips the GTF peek when exactly one annotation matches
  the requested genomic unit, and the ranges derived from an annotation are
  kept on the bundle, so building an SE and then an RSE from the same bundle
  parses the GTF once.

Fixed
~~~~~

- ``search`` subcommands reject ``key=value`` selectors the chosen mode does
  not read, instead of ignoring them. ``recount3 search project ...
  genomic_units=gene annotations=G026`` previously emitted the full
  ten-resource default set, because those three plural names belong to the
  Python ``R3ResourceBundle.discover`` API and not to the CLI, so copying
  Python keyword names to the shell produced a manifest that was wrong without
  any warning. The error names each unrecognized selector, the CLI spelling
  when the name is a known Python-API alias, and the selectors the mode does
  read.
- ``bundle se``/``bundle rse`` can write ``.h5ad`` again. The export called
  BiocPy's ``to_anndata()``, which hands experiment metadata to
  ``AnnData(uns=...)`` as a ``NamedList`` and fails with ``Only mutable
  mapping types (e.g. dict) are allowed for `.uns`.`` for every experiment the
  package builds. The export now goes through
  ``recount3.se.to_anndata``. ``--sanitize-columns`` also covers nested
  ``uns`` keys, which are the sample-metadata column names and carry the same
  forward slashes that HDF5 reads as path separators; each renamed provenance
  entry keeps its original name in its value.
- ``R3Resource.download(path=...)`` and
  ``R3ResourceBundle.download(dest=...)`` are annotated to accept any
  ``os.PathLike``, not only ``str``. Both already handled a ``pathlib.Path``
  at runtime, because each normalizes its argument with ``Path()`` before
  doing anything else. The narrow annotation nonetheless meant that feeding a
  destination straight back from ``ensure_cached()``, ``recount3_cache()``, or
  ``recount3_cache_files()``, all of which return ``Path``, was an error under
  mypy and pyright for anyone type-checking against the shipped ``py.typed``
  marker. The new ``recount3.StrPath`` alias names the accepted type. Return
  types are unchanged: ``download()`` still returns ``str | None``.
- ``R3Resource.filepath`` is annotated ``StrPath | None``, matching the
  ``os.PathLike`` the constructor already accepted and normalizes.
- ``R3Resource(filepath=...)`` normalizes an ``os.PathLike`` to ``str``.
  ``download()`` always stored a ``str``, so the attribute previously held
  either type depending on how it was populated: ``repr`` rendered
  ``PosixPath('...')`` for one and ``'...'`` for the other, and two resources
  naming the same file compared unequal.
- GTF strand values of ``.`` are normalized to ``*``, allowing unstranded
  annotations such as SIRV gene sums to build ranged experiments. Existing
  code remains source-compatible.
- Sample metadata whose tables share no samples no longer produces an
  experiment with zero samples. The inner-join check compared the number of
  matched samples with the number of merged rows, which are both zero when
  the join matches nothing, so every sample could be dropped and the caller
  handed a valid-looking object with an empty assay. Such a join now raises
  and names ``metadata_join="outer"`` as the alternative, and construction
  additionally refuses to return an experiment with no samples.
- Metadata tables are joined on identifier values rather than on a
  dtype-dependent rendering of them. A single blank cell widens a numeric
  ``rail_id`` column to ``float64``, so one table's ``"123488"`` met another
  table's ``"123488.0"`` and the two could not be joined; depending on the
  tables involved this surfaced either as a spurious "Conflicting
  external_id/rail_id mappings" error or as a silently emptied join. Junction
  ``.ID`` rail IDs, which become the count matrix's column labels, use the
  same canonical rendering.
- GTF attribute values are parsed with quoting rules. A quoted value ending
  at a ``;`` inside the quotes truncated the value and left its remainder to
  be rescanned as further ``key value`` pairs, so ``note "a, b; c (d)"``
  yielded ``note="a, b"`` plus a fabricated ``c`` column in ``rowData`` and
  ``rowRanges`` metadata. An empty ``key ""`` is now an empty value rather
  than an absent attribute. A repeated key still takes its last value.
- Metadata TSVs are parsed without quote processing, so literal ``NULL``,
  ``None``, quotes, and empty character fields survive as themselves instead
  of becoming missing values. Numeric blanks still become missing values.
- Metadata tables are merged within each project before samples are aligned,
  so a multi-project bundle no longer mixes one project's tables into
  another's samples.
- ``BigWigURL`` is built from each sample's own study and from the mirror its
  resource was configured with. Every sample previously received the URL of
  the first count resource's project, which is wrong for multi-project
  bundles, and the column always pointed at the default mirror even when the
  bundle was built against another one.
- Discovery no longer requests ``recount_pred`` metadata for GTEx and TCGA,
  which do not publish it. ``search_project_all(include_metadata=True)``
  enumerates five metadata tables for SRA and four for GTEx and TCGA.
- RSE annotation selection filters count matrices consistently with SE
  construction. Multi-project junctions align through each project's RR
  coordinates instead of local MM row numbers.
- GTF scores are preserved as ``bp_length``, and repeated exon occurrences
  retain their individual transcript annotations.
- Generated unique feature names cannot collide with a name already present
  in the annotation.
- Genomic ranges are built from the coordinate columns alone, with the
  remaining annotation columns attached positionally as ``mcols`` and the
  feature names set explicitly. The whole indexed frame was previously handed
  to ``GenomicRanges.from_pandas()``, where a non-default pandas index could
  misalign the metadata columns.
- ``create_rse`` no longer logs a spurious ``Failed to peek GTF features ...
  FileNotFoundError`` warning on the first use of an annotation. Annotation
  selection opened the computed cache path without first ensuring the file
  was there, so a cold cache produced a warning that vanished once the normal
  loading path had downloaded the file. Selection now prepares the cache
  through the existing download implementation, which already retries
  transient network errors. (#6)
- Comparing two loaded ``R3Resource`` objects no longer raises. Equality
  covers the fields that identify the file (``description``, ``url``,
  ``filepath``, and ``config``) and excludes the parsed object held on the
  resource, whose ``==`` is element-wise for the ``pandas.DataFrame`` that
  counts and metadata parse into.
- ``R3Resource.download(path=...)`` creates the destination directory when it
  does not exist and the file is already cached. Only the uncached path
  created it, so materializing a cached file into a directory that did not
  exist yet failed with ``FileNotFoundError``.
- Concurrent writes into one ZIP archive are serialized on Windows.
  ``Path.resolve()`` keeps the extended-length ``\\?\`` prefix whenever it
  cannot confirm that the plain spelling names the same file, which is exactly
  what happens while another thread holds that file open. The prefixed and
  plain spellings keyed separate locks, so one destination could be written
  through two locks at once. Lock keys now drop the prefix and normalize case.
  This has affected ``recount3 download --jobs`` into a ``.zip`` destination
  since 1.1.0.

Documentation
~~~~~~~~~~~~~

- Moved the CLI guide out of the ``recount3.cli`` module docstring and into
  ``docs/cli.rst``, so shell examples render with a single continuation
  backslash. Examples are marked as Bash code blocks, with single-line
  alternatives where useful. (#9)
- Added ``CONTRIBUTING.rst``, covering the development environment,
  repository layout, testing, code style, and how to submit changes.
- Added ``CITATION.cff`` and a README citation section, citing both the
  package preprint and the original recount3 data paper.
- Added a table of contents to the README.
- Documented the ``parquet`` and ``anndata`` extras, ``--densify``, and
  ``--sanitize-columns`` in the README, the tutorial, and the CLI reference.
- Documented threaded downloads, cache destination and refresh semantics,
  cross-process limits, and the optional shared R/Python cache, including its
  configuration precedence, maintenance, and measured overhead, in the
  tutorial's cache chapter.
- Documented ``metadata_join``, assay storage and dtypes, the ranges errors
  and the plain-``SummarizedExperiment`` fallback, and experiment provenance
  in the tutorial.
- Every authored example now uses ``import recount3 as r3`` and calls through
  that prefix, in the README, the tutorial, and the API docstrings. The
  tutorial states that builders and discovery helpers are reached directly
  (``r3.create_rse``) while the normalization helpers live on the ``se``
  submodule (``r3.se.compute_tpm``) and are not top-level exports.
- Added real, executed output to the tutorial's examples, covering assay
  types and shapes, normalized values, a sample correlation matrix, sparse
  junction storage figures, and an AnnData round trip, so the documented calls
  can be checked against what they actually return.
- Added a tutorial section on finding projects and samples with
  ``available_projects`` and ``available_samples`` before an accession is
  known, and a section on moving the returned NumPy, pandas, and SciPy
  objects into a downstream analysis.
- Corrected the README's bundle example, which filtered gene counts without
  naming an annotation. A bundle carrying more than one gene annotation now
  raises ``CompatibilityError``, so the example passes
  ``annotation_extension="G026"``.
- Corrected tutorial claims that did not match the implementation: only
  ``compute_scale_factors``, ``is_paired_end``, and ``expand_sra_attributes``
  accept a plain ``SummarizedExperiment``, while ``compute_read_counts``,
  ``transform_counts``, and ``compute_tpm`` require an RSE; discovery
  describes candidate URLs rather than validating existence or reporting
  sizes; ``compat="feature"`` does not pin an annotation build; and
  ``compute_tpm`` prefers ``bp_length`` over genomic span.
- Documented AnnData export: why the conversion is a package function rather
  than a BiocPy method call, what ``--sanitize-columns`` covers, and how to
  reach it from Python with ``r3.se.to_anndata``.
- Added a tutorial note that ``BigWigFile.load()`` returns the wrapper while
  entering it as a context manager yields the live ``pyBigWig`` handle.
- Corrected the description of ``create_sample_project_lists``: it returns a
  ``(samples, projects)`` pair of identifier lists for an organism; it is the
  CLI's ``recount3 ids`` that writes them out.
- Moved the project and sample discovery section ahead of the three API
  layers, where a reader without an accession needs it, and documented the
  remaining ``search_*`` helpers and ``R3Resource.from_mapping``.

Internal
~~~~~~~~

- Deprecation warnings originating in ``recount3`` modules now fail the test
  suite, so upstream deprecations surface before they become breaking
  changes.
- The coverage workflow enforces 100% statement and branch coverage over
  ``src``; new tests cover the optional CLI output paths (sparse Parquet
  densification, SE and RSE HDF5 column sanitization), cache concurrency and
  lifecycle, and the temp-file cleanup guard on the success path.
- Renamed the CI workflows for clarity (``ci.yml``, ``release.yml``,
  ``run-tests.yml``, ``publish-pypi.yml``, ``publish-docs.yml``,
  ``report-coverage.yml``), routed the per-platform install commands through
  the reusable test workflow, and installed the optional cache dependencies on
  every platform so the cache paths are exercised and counted.
- Development now requires ``pytest-cov>=7``. pytest-cov 5 aborts collection
  on Python 3.14 with "ImportError: numpy: cannot load module more than once
  per process", which broke coverage reporting on the ``3.x`` runner rather
  than degrading it.
- Closed dispatch paths in the core modules and in annotation-extension
  resolution use ``match`` statements.
- Tests build the local mirror URL with ``Path.as_uri()``.

1.1.0 (2026-06-12)
------------------

Added
~~~~~

- Multithreaded downloading in the API via
  ``R3ResourceBundle.download(max_workers=...)``, sharing the same thread-pool
  implementation as the CLI ``download`` command.
- ``available_samples``, ``available_projects``, and ``project_homes`` are now
  part of the public top-level ``recount3`` API.

Changed
~~~~~~~

- Increased the default download concurrency from 4 to 8 worker threads for
  the CLI (``--jobs``), matching the new API's 8 (``max_workers``).

Fixed
~~~~~

- ``create_rse(genomic_unit="junction")`` now succeeds by default. The junction
  extension defaults are unit-aware and include the ``RR`` sidecar required to
  attach genomic coordinates, so ``create_rse()`` returns a
  ``RangedSummarizedExperiment`` instead of raising ``ValueError``.
- SRA sample-attribute expansion now resolves both the R-style
  ``sra.sample_attributes`` column name if used and the namespaced
  ``sra__sample_attributes`` produced by the bundle layer.
- ``recount3 download`` now creates the ``--dest`` directory when it does not
  already exist (matching the existing behavior for ``.zip`` destinations).

Documentation
~~~~~~~~~~~~~

- Rewrote and expanded the tutorial to cover all three API layers, metadata
  merging, normalization and scaling, BigWig access, and cache management.
- Restructured and expanded the README, and corrected numerous docstrings,
  examples, and the CLI reference.

Internal
~~~~~~~~

- Routed resource-description construction through the resource layer,
  decoupling the bundle and CLI from the private ``_descriptions`` module.
- Continuous integration tests now also run on Windows and macOS (the
  ``bigwig`` extra is exercised on Linux only, ``biocpy`` on macOS and Linux
  only).

1.0.0 (2026-03-24)
------------------

Initial public release.
