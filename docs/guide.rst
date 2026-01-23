User Guide
==========

Overview
--------

recount3 provides a typed, minimal API for interacting with the recount3 repository.
The public surface is intentionally flat: core classes and helper functions are
re-exported for discovery and convenience.

Key Features
------------

* **Typed API**: Full type hints for better IDE support and reliability
* **Caching**: Configurable file caching to avoid redundant downloads
* **Resource Management**: Unified interface for all recount3 resource types
* **Bundle Operations**: Work with multiple resources simultaneously

Quick Start
-----------

.. code-block:: python

   from recount3 import (
       R3Resource, R3ResourceDescription, R3Annotations,
       search_count_files_gene_or_exon, R3ResourceBundle,
   )

   # Download and load gene annotations
   desc = R3Annotations(
       organism="human", genomic_unit="gene", annotation_file_extension="G026")
   res = R3Resource(desc)
   res.download(path=None, cache_mode="enable")
   df = res.load()

Configuration
-------------

Environment variables:

* ``RECOUNT3_URL``: Base URL for recount3 mirror (default: http://duffel.rail.bio/recount3/)
* ``RECOUNT3_CACHE_DIR``: Cache directory (default: ~/.cache/recount3/files)
* ``RECOUNT3_CACHE_DISABLE``: Disable caching (0 or 1)
* ``RECOUNT3_HTTP_TIMEOUT``: Network timeout in seconds (default: 60)
* ``RECOUNT3_MAX_RETRIES``: HTTP retry attempts (default: 3)