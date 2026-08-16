Changelog
=========

Unreleased
----------

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