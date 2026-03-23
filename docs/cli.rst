CLI Reference
=============

The ``recount3`` command-line tool implements a
discover -> manifest -> materialize workflow.

.. include:: ../src/recount3/cli.py
   :start-after: """
   :end-before: """


Full usage
----------

Run any subcommand with ``--help`` for the full option list:

.. code:: bash

   recount3 --help
   recount3 search --help
   recount3 download --help
   recount3 bundle stack-counts --help
   recount3 bundle se --help
   recount3 bundle rse --help
   recount3 smoke-test --help
