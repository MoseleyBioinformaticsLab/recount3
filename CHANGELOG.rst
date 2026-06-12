Changelog
=========

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