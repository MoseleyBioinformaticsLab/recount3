CLI Reference
=============

The ``recount3`` command-line tool implements a
discover -> manifest -> materialize workflow.

Summary
-------

Use ``recount3`` to:

* ``ids`` - Emit unique sample and project IDs.
* ``search`` - Discover resources and print a machine-readable manifest
  (JSONL or TSV).
* ``download`` - Materialize resources from a manifest (directory or
  ``.zip`` file).
* ``bundle`` - Operate on multiple resources (for example, stack count
  matrices).
* ``smoke-test`` - Run a small connectivity test for CI or local validation.

Quick start
-----------

Discover a handful of gene-level count files, save a manifest, and download
them:

.. code-block:: bash

   recount3 search gene-exon \
       organism=human data_source=sra genomic_unit=gene project=SRP009615 \
       --format=jsonl > manifest.jsonl

   recount3 download --from=manifest.jsonl --dest=./downloads --jobs=8

The search command can also be copied as a single line:

.. code-block:: bash

   recount3 search gene-exon organism=human data_source=sra genomic_unit=gene project=SRP009615 --format=jsonl > manifest.jsonl

Or stream a manifest directly, without an intermediate file:

.. code-block:: bash

   recount3 search annotations \
       organism=human genomic_unit=gene annotation_extension=G026 \
       --format=jsonl | \
       recount3 download --from=- --dest=./annots

The same pipeline on one line is:

.. code-block:: bash

   recount3 search annotations organism=human genomic_unit=gene annotation_extension=G026 --format=jsonl | recount3 download --from=- --dest=./annots

Commands
--------

``ids``
  Emit unique ID lists. By default, output is written to standard output.

  Flags:

  .. code-block:: text

     --organism=human|mouse|""   Empty means all organisms.
     --samples-out=<file>        Write samples to a text file (else stdout).
     --projects-out=<file>       Write projects to a text file (else stdout).

``search``
  Discover resources and print a manifest (JSONL or TSV). Filters are passed
  as space-separated ``key=value`` tokens.

  By default, results are written to standard output for easy piping. Use
  ``--output <file>`` to write a specific file, or ``--outdir <dir>`` to
  create a timestamped filename in that directory.

  Modes and required filters:

  .. code-block:: text

     annotations   organism, genomic_unit, annotation_extension
     gene-exon     organism, data_source, genomic_unit, project
                   (optional: annotation_extension; default G026)
     junctions     organism, data_source, project
                   (optional: junction_type=ALL, junction_extension=MM)
     metadata      organism, data_source, table_name, project
     bigwig        organism, data_source, project, sample
     project       organism, data_source, project
                   (optional: genomic_unit=gene,exon;
                   annotation=default|all|gencode_v26,gencode_v29
                     Human-readable name or convenience alias.
                     'default' -> primary annotation (G026 for human,
                     M023 for mouse). 'all' -> every available annotation.
                     A comma list of names (for example, gencode_v26) or raw
                     extension codes (for example, G026) also works.
                   annotation_extension=G026,G029
                     Raw annotation file-extension codes. When set, this
                     overrides 'annotation' completely. Use it when you
                     already know the exact code(s) you need.
                   junction_type=ALL;
                   junction_extension=MM,RR,ID;
                   include_metadata=true|false;
                   include_bigwig=true|false)
     sources       organism
     source-meta   organism, data_source

  Example:

  .. code-block:: bash

     recount3 search junctions \
         organism=human data_source=sra project=SRP009615 \
         junction_type=ALL junction_extension=MM --format=tsv

``download``
  Materialize resources from a manifest file or one inline JSON object.
  Writes one JSONL progress event per resource to standard output.

  Source:

  .. code-block:: text

     --from=<path>|-       Read JSONL manifest from a file or stdin ('-').
     --inline='<json>'     One JSON object for a single resource.

  Destination:

  .. code-block:: text

     --dest=<dir-or-zip>   Directory or .zip file path.
     --overwrite           Overwrite existing files (directory mode only).

  Behavior:

  .. code-block:: text

     --jobs=<n>            Maximum parallel downloads (default 8).
     --cache=MODE          Cache behavior (default: enable). MODE is one of:
                           enable - use cache; disable - bypass cache;
                           update - force re-download, then cache.

``bundle stack-counts``
  Concatenate compatible count matrices (gene/exon or junctions).

  Required:

  .. code-block:: text

     --from=<manifest>     JSONL manifest (or '-' for stdin).
     --out=<path>          Output file (.csv, .tsv, .tsv.gz, or .parquet).
                           .parquet needs: pip install "recount3[parquet]"

  Options:

  .. code-block:: text

     --compat=family|feature    Compatibility mode (default: family).
     --join=inner|outer         Pandas join type (default: inner).
     --axis=0|1                 Concatenate rows (0) or columns (1).
     --verify-integrity         Fail on duplicate index after concat.
     --densify                  Densify sparse columns for Parquet output.

``bundle se`` / ``bundle rse``
  Assemble a (Ranged)SummarizedExperiment. This requires the
  ``recount3[biocpy]`` extra; ``.h5ad`` output also needs the
  ``recount3[anndata]`` extra.

  Required:

  .. code-block:: text

     --from=<manifest>     JSONL manifest (or '-' for stdin).
     --genomic-unit=gene|exon|junction
     --out=<path>          Output file (.pkl or .h5ad).

  Options:

  .. code-block:: text

     --sanitize-columns    Replace '/' with '_' in obs/var column names,
                           which .h5ad output requires (see below).

  HDF5 reads a forward slash as a path separator, and recount3 STAR QC
  fields are named after splice motifs (for example,
  ``recount_qc__star.number_of_splices:_gt/ag``), so ``.h5ad`` output fails
  for most projects carrying SRA metadata unless those columns are renamed.
  ``--sanitize-columns`` performs the rename and logs every one of them.
  Writing a ``.pkl`` keeps the names verbatim.

``smoke-test``
  Download a few tiny files to verify connectivity and configuration.

  Options:

  .. code-block:: text

     --dest=<dir>          Destination directory (default ./recount3-smoke).
     --limit=<n>           Number of resources to attempt (default 1).

Input and output formats
------------------------

**JSONL** (also called NDJSON) stores one JSON object per line. It is ideal for
streaming, grepping, and piping, and is used for both ``search`` output and
``download`` input.

Each manifest line contains all resource description fields plus two
convenience keys: ``url`` (the fully qualified HTTP URL) and ``arcname`` (the
destination path inside a ``.zip`` archive). For example, this record is
wrapped for readability:

.. code-block:: json

   {"resource_type":"count_files_gene_or_exon","organism":"human",
    "data_source":"sra","genomic_unit":"gene","project":"SRP009615",
    "sample":null,"annotation_extension":"G026","junction_type":null,
    "junction_extension":null,"table_name":null,
    "url":".../sra/gene_sums/15/SRP009615/sra.gene_sums.SRP009615.G026.gz",
    "arcname":"human/data_sources/.../sra.gene_sums.SRP009615.G026.gz"}

``download`` writes one progress event per resource to standard output:

.. code-block:: json

   {"url":"...","status":"ok","dest":"/path/to/file"}
   {"url":"...","status":"skipped","dest":"/existing/file"}
   {"url":"...","status":"error","dest":null,"error":"<repr>"}

**TSV** is tab-separated text for quick human scanning or spreadsheet import.
TSV is available for ``search --format=tsv`` only; ``download`` expects JSONL.

Configuration
-------------

Configuration is centralized in :class:`recount3.config.Config`. Values come
from, in order of decreasing precedence, CLI flags, environment variables,
and library defaults. The relevant environment variables are:

.. code-block:: text

   RECOUNT3_URL               Base URL (trailing slash added automatically)
   RECOUNT3_CACHE_DIR         Directory for the on-disk cache
   RECOUNT3_CACHE_BACKEND     filesystem (default) or pybiocfilecache
   RECOUNT3_CACHE_DISABLE     "1" disables cache, anything else enables
   RECOUNT3_HTTP_TIMEOUT      HTTP timeout in seconds (int)
   RECOUNT3_MAX_RETRIES       Max retry attempts for transient errors (int)
   RECOUNT3_INSECURE_SSL      "1" to disable TLS verification (unsafe; https
                             base URLs only, no-op for default http mirror)
   RECOUNT3_USER_AGENT        Custom HTTP User-Agent string
   RECOUNT3_CHUNK_SIZE        Streaming chunk size in bytes

Global flags mirror these settings: ``--base-url``, ``--cache-dir``,
``--cache-backend``,
``--timeout``, ``--retries``, ``--insecure-ssl``, ``--user-agent``, and
``--chunk-size``.

With ``--cache-backend pybiocfilecache`` (or the equivalent environment
setting), the default directory matches R's recount3 cache: normally
``~/.cache/R/recount3`` on Linux/WSL. R's ``R_USER_CACHE_DIR`` and
``XDG_CACHE_HOME`` overrides and platform-specific defaults are respected.
An explicit ``--cache-dir`` or ``RECOUNT3_CACHE_DIR`` still takes precedence.
The ``filesystem`` default remains ``~/.cache/recount3/files``. See
:ref:`cache-and-configuration` for platform paths and shared-cache examples.

Logging
-------

Logging defaults to INFO. Use ``--quiet`` for WARNING or ``--verbose`` for
DEBUG. Log messages follow pattern-string formatting (not f-strings), per the
Google guide, and include greppable context such as ``url=...`` and
``dest=...``.

Exit codes
----------

.. code-block:: text

   0    Success
   1    Malformed --inline JSON in download
   2    Fatal error (missing filters, I/O failures, bad configuration; also
        argparse validation errors such as unrecognized flags)
   3    Partial failure in download (some items failed)
   130  Interrupted (Ctrl-C)

Security and safety
-------------------

* TLS verification is on by default for ``https://`` mirrors.
  ``--insecure-ssl`` disables it and should only be used to debug certificate
  issues. It applies only to ``https://`` base URLs and is a no-op for the
  default ``http://`` Duffel mirror; the AWS Open Data and JHU IDIES https
  mirrors have valid certificates and need no flag.
* The cache reduces repeated downloads. Choose ``--cache=disable`` to bypass
  it when correctness requires a direct fetch.

Performance tips
----------------

* Increase ``--jobs`` to improve throughput when network-bound.
* Keep the cache enabled for repeated workflows.
* Use streaming pipelines with JSONL and standard tools (``jq``, ``grep``,
  ``head``, and ``tail``) to avoid loading everything into memory.

Example recipes
---------------

List human SRA data sources, then download their metadata:

.. code-block:: bash

   recount3 search sources organism=human --format=jsonl > sources.jsonl
   recount3 search source-meta organism=human data_source=sra \
       --format=jsonl > meta.jsonl
   recount3 download --from=meta.jsonl --dest=./meta

Stack gene-level matrices across samples and write Parquet (first install the
engine with ``pip install "recount3[parquet]"``):

.. code-block:: bash

   recount3 search gene-exon \
       organism=human data_source=sra genomic_unit=gene project=SRP009615 \
       --format=jsonl > counts.jsonl
   recount3 bundle stack-counts --from=counts.jsonl --compat=family \
       --join=inner --axis=1 --out=counts.parquet

Junction matrices load sparse, and no Parquet engine accepts pandas sparse
dtypes. Either keep a text format or densify explicitly:

.. code-block:: bash

   recount3 bundle stack-counts --from=junctions.jsonl \
       --out=junctions.tsv.gz
   recount3 bundle stack-counts --from=junctions.jsonl --densify \
       --out=junctions.parquet

Troubleshooting
---------------

* "Missing required filters": Check the mode-specific filter list above.
* "Writing Parquet requires a Parquet engine": Install one with
  ``pip install "recount3[parquet]"``, or write ``.tsv``, ``.tsv.gz``, or
  ``.csv`` instead.
* "Cannot write .h5ad: Optional dependency ... is required": Install
  ``pip install "recount3[anndata]"``, or write a ``.pkl`` file instead.
* "Cannot write .h5ad: ... column name(s) contain a forward slash": HDF5
  reads ``/`` as a path separator. Add ``--sanitize-columns`` to rename
  those columns, or write a ``.pkl`` file instead.
* "Cannot write Parquet: ... columns use a pandas sparse dtype": You are
  stacking junctions. Add ``--densify`` (which materializes every zero) or
  write a text format.
* ``json.JSONDecodeError``: Ensure your manifest is valid JSONL. Each line
  must be one JSON object.
* Permission or path errors: Verify ``--dest`` exists (or its parent for a
  ``.zip`` file) and is writable; on shared filesystems, reduce ``--jobs`` to
  avoid pressure.
* TLS/SSL errors: Try updating CA certificates, or as a last resort temporarily
  use ``--insecure-ssl`` to isolate the issue.


Full usage
----------

Run any subcommand with ``--help`` for the full option list:

.. code-block:: bash

   recount3 --help
   recount3 search --help
   recount3 download --help
   recount3 bundle stack-counts --help
   recount3 bundle se --help
   recount3 bundle rse --help
   recount3 smoke-test --help
