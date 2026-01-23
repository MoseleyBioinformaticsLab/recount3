Command Line Interface
======================

The recount3 package includes a beta command-line interface for basic tasks.

Basic Usage
-----------

.. code-block:: bash

   # Generate sample and project lists
   python -m recount3

   # This creates:
   # - samplist.txt: List of sample identifiers
   # - projlist.txt: List of project identifiers  
   # - downloads/: Test downloads of example resources

Commands
--------

Currently, the CLI provides:

* **Sample/Project Discovery**: Generate comprehensive lists of available data
* **Smoke Testing**: Download test files to verify connectivity and configuration

Environment Variables
--------------------

The CLI respects the same environment variables as the Python API:

.. code-block:: bash

   export RECOUNT3_CACHE_DIR="/custom/cache/path"
   export RECOUNT3_URL="https://custom.mirror/recount3/"
   python -m recount3