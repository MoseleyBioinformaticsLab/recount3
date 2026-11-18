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
3.15. Install in editable mode with the development extras:

.. code:: bash

   git clone https://github.com/MoseleyBioinformaticsLab/recount3.git
   cd recount3
   python3 -m venv .venv
   . .venv/bin/activate
   python3 -m pip install --upgrade pip
   python3 -m pip install -e ".[dev,bigwig,biocpy,parquet,anndata]"

The optional extras are not installable on every platform. Use the command
above on Linux. On macOS, omit ``bigwig``; on Windows, install
``-e ".[dev,parquet]"``. On Python 3.15, also omit ``anndata``, whose ``h5py``
dependency publishes no Python 3.15 wheels yet. Continuous integration applies
the same configurations, and tests that require an unavailable extra are
skipped automatically.


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

* **Line endings.** Use LF line endings and end text files with a newline.
* **Formatting.** Black at a line length of 80, as configured in
  ``[tool.black]``. Black is not part of the ``dev`` extra; install it
  separately and run ``black src tests``. Black does not reflow docstrings or
  comments; wrap those by hand at the same width.
* **Linting.** Ruff is included in the ``dev`` extra. From the repository
  root, run ``python -m ruff check .``; CI runs the same command. This checks
  package source, tests, the Sphinx configuration, and ``pyproject.toml``.
  Configuration lives in ``[tool.ruff]`` and ``[tool.ruff.lint]``, with a
  Python 3.10 target, an 80-character line length, and Google-style
  docstrings and import ordering. Black remains the formatter; Ruff handles
  linting and import sorting. To sort imports, run
  ``python -m ruff check --select I --fix .`` before Black.
* **Typing.** Every function and method is fully annotated. Modules begin with
  ``from __future__ import annotations``. The package ships ``py.typed``, so
  annotations are part of the public contract.
* **Imports.** Import modules, not individual names, as the Google Python
  Style Guide specifies: ``from recount3 import resource``, then
  ``resource.R3Resource``. ``typing``, ``typing_extensions``, and
  ``collections.abc`` are exempt, and ``src/recount3/__init__.py`` imports
  names directly because it re-exports them. Standard scientific aliases
  are ``numpy as np`` and ``pandas as pd``; keep ``datetime`` spelled out.
  Other aliases follow Google's name-collision and readability exceptions:
  ``types`` is imported as ``r3_types`` (the standard library module), and
  ``config`` and ``bundle`` as ``r3_config`` and ``r3_bundle`` (public
  parameter names). Use the same alias in every module, and add any new alias
  to ``_MODULE_ALIASES`` in ``docs/conf.py`` so that annotations using it
  still link in the API reference. Ruff's ``I001`` sorts imports by package
  path within each group and keeps module imports on separate lines, with
  exceptions for the typing modules above. ``ICN001`` enforces the configured
  scientific aliases. ``ICN003`` requires module imports for the modules in
  ``banned-from``; review other imports against the same Google convention.
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


The Ruff configuration enables all stable rules (``select = ["ALL"]``), then
ignores individual rules with existing findings. Each ignored rule is named
in ``pyproject.toml``. Additional exceptions specific to tests and Sphinx are
scoped to those directories; for example, assertions are permitted in tests
while unused imports are still checked in package source. This establishes a
passing baseline without changing runtime code. Preview rules are disabled.

To tighten the baseline, remove one ignore, run ``python -m ruff check .``,
resolve its findings, and submit that as a separate change. Keep ignores
specific to individual rule codes, and use a narrow per-file exception when
it applies only to tests or documentation. Avoid blanket ``noqa`` comments.
Ruff is constrained to the audited minor release series in both the ``dev``
dependency and ``required-version``. When upgrading that series, review the
`Ruff versioning policy <https://docs.astral.sh/ruff/versioning/>`_ and audit
newly enabled rules before updating the baseline.


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
   Black-formatted, ``python -m ruff check .`` passes, and the documentation
   builds.
4. Update ``CHANGELOG.rst``.
5. Open the pull request against ``main``, describing the change, its
   motivation, and any change in public behaviour. Link the relevant issue.

CI runs Ruff in a dedicated job and the suite across the supported Python
versions on Linux, macOS, and Windows. It enforces the coverage threshold on
Python 3.14, the newest version on which every extra installs. All checks must
pass before review concludes.


Releases
~~~~~~~~

Releases are prepared by the maintainers. The package follows semantic
versioning. The procedure is:

1. Set the new version in ``src/recount3/version.py`` (``__version__``),
   ``tests/test_version.py`` (the expected version), and ``CITATION.cff``
   (``version``, with ``date-released`` set to the release date). The release
   workflow runs the suite only after publishing to TestPyPI, so a stale
   expected version fails the release there.
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
