TODO
====

1.5.0
-----

Add
~~~

- Add ``compute_fpkm()`` and ``compute_fpkm_uq()`` to ``recount3.compute``.

Change
~~~~~~

- Improve CLI diagnostics with actionable errors, manifest paths and line
  numbers, consistent exit codes, and tracebacks under ``--verbose``. Report
  readable download progress on stderr while preserving stdout formats.

Internal
~~~~~~~~

- Split ``_utils.py`` into ``_cache.py``, ``_filesystem.py``,
  ``_download.py`` (HTTP and ZIP), ``_metadata.py``, ``_optional.py`` (lazy
  imports), and ``_export.py`` (AnnData and Parquet). Move count assay
  selection to ``compute.py``, genomic-unit validation to ``bundle.py``,
  and junction sidecar URL construction to ``resource.py``.

2.0.0
-----

Breaking changes
~~~~~~~~~~~~~~~~

- Remove the five deprecated ``recount3.se`` wrappers for functions moved to
  ``recount3.compute``.

Unscheduled
-----------

Internal
~~~~~~~~

- Verify collection discovery and assembly on the public recount3 mirrors
  when they publish collection files.
