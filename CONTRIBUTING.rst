Contributing to recount3
========================

Contributions of bug reports, feature requests, documentation, and code are
welcome. This document specifies the development environment, the conventions
enforced in this repository, and the procedure for submitting changes.

Discussion and issue tracking take place at
https://github.com/MoseleyBioinformaticsLab/recount3/issues.


Contents
~~~~~~~~

* `Reporting bugs`_
* `Requesting features`_
* `Development environment`_
* `Repository layout`_
* `Testing`_
* `Code style`_
* `Documentation`_
* `Changelog`_
* `Commit messages`_
* `Pull requests`_
* `Releases`_
* `License`_


Reporting bugs
~~~~~~~~~~~~~~

Open an issue that includes:

* The package version (``recount3 --version``), the Python version, and the
  operating system.
* The optional extras installed (``bigwig``, ``biocpy``, ``parquet``,
  ``anndata``), since several code paths depend on them.
* A minimal, self-contained reproducer, and the complete traceback or CLI
  output.
* The observed behaviour and the expected behaviour.

Report defects attributable to the upstream recount3 data resource or its
mirrors to the resource maintainers rather than here. Defects in this package's
handling of that data belong in this tracker.


Requesting features
~~~~~~~~~~~~~~~~~~~

Open an issue describing the intended use case, the proposed API or CLI
surface, and the alternatives considered. Features that depend on data not yet
published on the public recount3 mirrors are deferred until that data is
available; see ``TODO.rst``.


Development environment
~~~~~~~~~~~~~~~~~~~~~~~

The package requires Python 3.10 or newer and is tested on Python 3.10 through
3.14. Install in editable mode with the development extras:

.. code:: bash

   git clone https://github.com/MoseleyBioinformaticsLab/recount3.git
   cd recount3
   python3 -m venv .venv
   . .venv/bin/activate
   python3 -m pip install --upgrade pip
   python3 -m pip install -e ".[dev,bigwig,biocpy,parquet,anndata]"

The optional extras are not installable on every platform. Use the command
above on Linux. On macOS, omit ``bigwig``; on Windows, install
``-e ".[dev,parquet]"``. Continuous integration applies the same three
configurations, and tests that require an unavailable extra are skipped
automatically.


Repository layout
~~~~~~~~~~~~~~~~~

==========================  ==================================================
Path                        Contents
==========================  ==================================================
``src/recount3/``           Package source. Modules prefixed with ``_`` are
                            internal.
``tests/``                  Test suite; one ``test_<module>.py`` per source
                            module.
``tests/data/``             Local mirror of a small subset of the recount3
                            raw-files tree.
``docs/``                   Sphinx documentation sources.
``.github/workflows/``      Continuous integration and release workflows.
==========================  ==================================================


Testing
~~~~~~~

Run the suite from the repository root:

.. code:: bash

   python -m pytest              # full suite
   python -m pytest -n auto      # parallel, via pytest-xdist

Configuration resides in ``[tool.pytest.ini_options]`` of ``pyproject.toml``;
no additional flags are required.

**The suite must not access the network.** Tests read from the local mirror in
``tests/data/recount3_mirror/`` by pointing ``RECOUNT3_URL`` at it, construct
fixtures on disk, or substitute test doubles with ``pytest-mock``. A test that
contacts a live mirror is not acceptable.

**Coverage must remain at 100 percent.** CI fails the build otherwise. Verify
locally before submitting:

.. code:: bash

   python -m pytest --cov=src --cov-report=term-missing
   coverage report --fail-under=100

Lines that are genuinely unreachable under test may be excluded with
``# pragma: no cover``; use the exclusion only when the alternative is an
artificial test.

Tests that require an optional dependency must carry the corresponding marker,
which is registered and applied in ``tests/conftest.py``:
``requires_pybigwig``, ``requires_biocpy``, ``requires_parquet``, and
``requires_anndata``.

``DeprecationWarning`` and ``PendingDeprecationWarning`` raised from
``recount3`` modules are configured as errors, so upstream deprecations fail
the suite rather than accumulating silently. Resolve them rather than
suppressing them.


Code style
~~~~~~~~~~

* **Formatting.** Black at a line length of 80, as configured in
  ``[tool.black]``. Black is not part of the ``dev`` extra; install it
  separately and run ``black src tests``. Black does not reflow docstrings or
  comments; wrap those by hand at the same width.
* **Linting.** Pylint, from the ``dev`` extra, is run at its default
  configuration: ``python -m pylint src/recount3``. The existing code is not
  free of messages; those that remain are design-complexity reports
  (``too-many-arguments``, ``too-many-locals``, and similar) and over-length
  lines in prose. Introduce no new messages, and prefer resolving any that
  your change touches.
* **Typing.** Every function and method is fully annotated. Modules begin with
  ``from __future__ import annotations``. The package ships ``py.typed``, so
  annotations are part of the public contract.
* **Docstrings.** Google-style sections (``Args``, ``Returns``, ``Raises``),
  rendered by ``sphinx.ext.napoleon``. Reference other objects with Sphinx
  roles, for example ``:class:`~recount3.resource.R3Resource```. Every public
  module, class, and function carries a docstring.
* **Public surface.** Names intended for public use are re-exported from
  ``src/recount3/__init__.py`` and listed in its ``__all__``.
* **License header.** Every ``.py`` file in ``src/`` and ``tests/`` begins with
  the Clear BSD license header. Copy it verbatim from an existing file when
  adding a new one.
* **Errors.** Raise the typed exceptions defined in ``recount3.errors`` rather
  than bare built-ins, so that failure modes remain catchable by category.


Documentation
~~~~~~~~~~~~~

Build the documentation with:

.. code:: bash

   python3 -m pip install -e ".[docs]"
   sphinx-build docs docs/_build/html

A change to public behaviour requires a corresponding documentation change:

* API additions are documented in their docstrings. A new module additionally
  requires an ``automodule`` entry in ``docs/api.rst``.
* CLI changes are reflected in ``docs/cli.rst``.
* Changes affecting installation, dependencies, or the introductory examples
  are reflected in ``README.rst`` and, where relevant, ``docs/tutorial.rst``.


Changelog
~~~~~~~~~

``CHANGELOG.rst`` follows the Keep a Changelog convention. Add an entry under
the ``Unreleased`` heading, in the appropriate category (``Added``,
``Changed``, ``Fixed``, ``Removed``). State what changed and why it matters to
a user, and note explicitly whether existing code remains source-compatible.
Released sections are headed ``<version> (<YYYY-MM-DD>)`` and are not edited
after release.


Commit messages
~~~~~~~~~~~~~~~

Commits follow the Conventional Commits format:

.. code:: text

   type(scope): imperative summary

   Body explaining the motivation and the consequences of the change.

   Closes #123

Types in use are ``feat``, ``fix``, ``docs``, ``test``, ``refactor``,
``style``, ``build``, and ``chore``. The scope names the affected module or
area, for example ``fix(bundle):``. Keep the summary in the imperative mood and
at or under 72 characters. Reference the issue the commit resolves.


Pull requests
~~~~~~~~~~~~~

1. Branch from ``main``.
2. Confine the branch to one logical change. Submit unrelated changes
   separately.
3. Ensure the full suite passes, coverage remains at 100 percent, the code is
   Black-formatted, Pylint reports no new messages, and the documentation
   builds.
4. Update ``CHANGELOG.rst``.
5. Open the pull request against ``main``, describing the change, its
   motivation, and any change in public behaviour. Link the relevant issue.

CI runs the suite across the supported Python versions on Linux, macOS, and
Windows, and enforces the coverage threshold. All checks must pass before
review concludes.


Releases
~~~~~~~~

Releases are prepared by the maintainers. The package follows semantic
versioning. The procedure is:

1. Update ``__version__`` in ``src/recount3/version.py``.
2. Retitle the ``Unreleased`` section of ``CHANGELOG.rst`` as
   ``<version> (<YYYY-MM-DD>)``.
3. Publish a GitHub release tagged ``v<version>``.

The release workflow then publishes to TestPyPI, runs the suite against that
artifact, publishes to PyPI, runs the suite again against the published
package, and deploys the documentation to GitHub Pages.


License
~~~~~~~

The package is distributed under The Clear BSD License with an additional
citation clause; see ``LICENSE``. Submitting a contribution constitutes
agreement that it is licensed under those terms.
