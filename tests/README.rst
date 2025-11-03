Tests
=====

This folder contains the unit/integration tests for recount3.

How to run
----------

.. code-block:: console

   python -m pip install -e ".[dev]"
   pytest -q

Markers
-------

- ``slow``: opt-in, longer tests.
- ``docs``: opt-in, builds Sphinx docs (requires extras ``docs``).
- ``net``: opt-in, allows network (default is blocked).

Opt-in examples:

.. code-block:: console

   pytest -m "slow or docs" -q
   pytest -m net -q
