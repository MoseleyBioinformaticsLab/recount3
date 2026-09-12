Changelog
=========

Unreleased
----------

Added
~~~~~

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

Changed
~~~~~~~

- MatrixMarket junction files are now parsed into a ``scipy.sparse.csr_array``
  instead of a ``scipy.sparse.csr_matrix``. SciPy 1.18 deprecated the implicit
  sparse-matrix return of ``scipy.io.mmread`` ahead of flipping the default in
  SciPy 1.20; ``recount3`` now requests the new behaviour explicitly where the
  installed SciPy supports it and normalises the result, so the parsed type is
  identical across the whole supported SciPy range. The keyword is
  feature-detected, keeping SciPy releases older than 1.18 (the newest
  installable on Python 3.10 and 3.11) working unchanged.
- Deprecation warnings originating in ``recount3`` modules now fail the test
  suite, so upstream deprecations surface before they become breaking changes.
- When genomic ranges cannot be derived, the reported reason now distinguishes
  a failure to retrieve the annotation, an annotation that cannot be parsed,
  and an annotation that does not cover every counted feature. Only the last
  is fixed by choosing a different ``annotation_extension``.
- The "falling back" warning is emitted only when a plain
  ``SummarizedExperiment`` is actually returned, instead of whenever range
  derivation failed.
- ``autoload`` now reaches annotation selection, so ``autoload=False``
  inspects only already-cached annotations and never downloads.
- ``to_ranged_summarized_experiment`` raises ``RangesError`` instead of a bare
  ``ValueError``. ``RangesError`` is a ``ValueError``, so this is
  source-compatible.

Fixed
~~~~~

- ``create_rse`` no longer logs a spurious ``Failed to peek GTF features ...
  FileNotFoundError`` warning on the first use of an annotation. Annotation
  selection opened the computed cache path without first ensuring the file was
  there, so a cold cache produced a warning that vanished once the normal
  loading path had downloaded the file. Selection now prepares the cache
  through the existing download implementation, which already retries
  transient network errors.

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
~~~~~~~~~~~~~~

- Rewrote and expanded the tutorial to cover all three API layers, metadata
  merging, normalization and scaling, BigWig access, and cache management.
- Restructured and expanded the README, and corrected numerous docstrings,
  examples, and the CLI reference.

Internal
~~~~~~~~

- Routed resource-description construction through the resource layer,
  decoupling the bundle and CLI from the private ``_descriptions`` module.
- Continuous integration tests now also runs on Windows and macOS (``bigwig``
  extra is exercised on Linux only, ``biocpy`` on macOS and Linux only).

1.0.0 (2026-03-24)
------------------

Initial public release.