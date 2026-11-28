Contributing to recount3
========================

Contributions of bug reports, feature requests, documentation, and code are
welcome. This document specifies the development environment, the conventions
enforced in this repository, and the procedure for submitting changes.

Discussion and issue tracking take place at
https://github.com/MoseleyBioinformaticsLab/recount3/issues.


.. contents:: Contents
   :local:
   :depth: 1
   :backlinks: none


Reporting bugs
~~~~~~~~~~~~~~

Open an issue that includes:

* The package version (``recount3 --version``), the Python version, and the
  operating system.
* The optional extras installed (``anndata``, ``bigwig``, ``biocpy``, ``parquet``,
  ``pybiocfilecache``), since several code paths depend on them.
* A minimal, self-contained reproducer and the complete traceback or CLI
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
3.15. On Linux through Python 3.14, install the package in editable mode with
the development tools and all feature extras.

.. code:: bash

   git clone https://github.com/MoseleyBioinformaticsLab/recount3.git
   cd recount3
   python3 -m venv .venv
   . .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install -e ".[dev,all]"

For other setups, use the extras from the corresponding CI job in the
editable install command.

* Linux on Python 3.15 uses ``.[dev,bigwig,biocpy,parquet,pybiocfilecache]``.
* macOS uses ``.[dev,anndata,biocpy,parquet,pybiocfilecache]`` through Python
  3.14 and omits ``anndata`` on Python 3.15.
* Windows uses ``.[dev,parquet,pybiocfilecache]``.

The test suite skips tests that require unavailable extras.


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

Keep tests, mocks, and module-specific fixtures and helpers in the matching
``test_<module>.py`` file. When two modules need the same mock, define it
locally in each test file.

**The suite must not access the network.** Tests read from the local mirror in
``tests/data/recount3_mirror/recount3/`` by pointing ``RECOUNT3_URL`` at its
``file://`` URL, construct fixtures on disk, or substitute test doubles with
``pytest-mock``.

**Coverage must remain at 100 percent.** CI fails the build otherwise. Verify
locally before submitting:

.. code:: bash

   python -m pytest --cov=src --cov-report=term-missing
   coverage report --fail-under=100

Lines that are genuinely unreachable under test may be excluded with
``# pragma: no cover``; use the exclusion only when the alternative is an
artificial test.

Use the optional-dependency markers defined in ``tests/conftest.py`` where
they apply. The markers are ``requires_pybigwig``, ``requires_biocpy``,
``requires_parquet``, and ``requires_anndata``. Tests that need
``pybiocfilecache`` use ``pytest.importorskip``.

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
  are ``numpy as np`` and ``pandas as pd``. Keep ``datetime`` spelled out.
  Import ``recount3.types`` as ``r3_types`` to distinguish it from the
  standard library's ``types`` module. The ``r3_config`` and ``r3_bundle``
  aliases avoid collisions with public parameter names. Use the same alias
  in every module and add new aliases to ``_MODULE_ALIASES`` in
  ``docs/conf.py`` so that annotations link in the API reference.
  Ruff's ``I001`` sorts imports by package path within each group and keeps
  module imports on separate lines, with exceptions for the typing modules
  above. ``ICN001`` enforces the configured scientific aliases. ``ICN003``
  requires module imports for the modules in ``banned-from``. Review other
  imports against the same Google convention.
* **Docstrings.** Google-style sections (``Args``, ``Returns``, ``Raises``),
  rendered by ``sphinx.ext.napoleon``. Reference other objects with Sphinx
  roles, for example ``:class:`~recount3.resource.R3Resource```. Every public
  module, class, and function carries a docstring.
* **Public surface.** Names intended for public use are re-exported from
  ``src/recount3/__init__.py`` and listed in its ``__all__``. Each commented
  group lists names from shortest to longest, with case-insensitive
  alphabetical tie-breaking.
* **License header.** Every ``.py`` file in ``src/`` and ``tests/`` begins with
  the Clear BSD license header. Copy it verbatim from an existing file when
  adding a new one.
* **Errors.** Raise the typed exceptions defined in ``recount3.errors`` rather
  than bare built-ins, so that failure modes remain catchable by category.


Ruff exceptions
~~~~~~~~~~~~~~~

The Ruff configuration in ``pyproject.toml`` selects all stable rules
applicable under the Google docstring convention. Preview rules remain off
while Ruff develops them. The test and Sphinx exceptions appear separately
under ``per-file-ignores``. The test entries cover fixtures, mocks,
assertions, and test-specific formatting. The docs entries cover Sphinx
configuration conventions. The reasons below apply to package code.

Global ignores
^^^^^^^^^^^^^^

* ``ANN401`` (explicit ``Any``). Optional backends and dynamically loaded
  objects expose several types at the same boundary. For example,
  ``expand_sra_attributes`` accepts a pandas DataFrame or a BiocPy
  experiment, and AnnData metadata can contain heterogeneous values.
  ``Any`` describes those dynamic inputs where a narrower type would
  misstate what the code accepts.
* ``BLE001`` (broad exception handlers). Download boundaries translate
  errors from several backends into ``DownloadError`` with the original
  cause. Best-effort bundle operations can skip a failed resource, and
  the CLI reports per-resource failures while processing the rest.
* ``C901`` (cyclomatic complexity). Range alignment and annotation
  selection branch on resource type, annotation coverage, and
  compatibility policies. Their branch counts reflect supported input
  cases that benefit from seeing the decision order in one place.
* ``COM812`` (missing trailing comma). Black controls trailing commas
  and wrapping. It can split an expression across lines even when its
  surrounding call fits with the last argument and closing parenthesis on
  the same line. Requiring an extra comma in those calls can force another
  layout.
* ``EM101`` (string literal in an exception). Fixed validation messages
  stay beside the condition that raises them, such as the check for
  missing RR columns. A separate name adds indirection for a one-use
  message.
* ``EM102`` (f-string in an exception). Values such as a resource URL,
  annotation code, or observed matrix shape are formatted at the point
  of failure. This keeps the condition and the reported context together.
* ``FBT001`` (positional Boolean annotation). Public functions such as
  ``compute_read_counts`` and ``transform_counts`` accept an existing
  positional ``round_to_integers`` argument. Making it keyword-only
  would break positional calls.
* ``FBT002`` (positional Boolean default). The same public rounding
  arguments default to ``True``. They remain optional so callers can
  omit the argument or pass it positionally.
* ``LOG015`` (root logger calls). The CLI configures the root logger and
  includes the logger name in its diagnostic format. The root calls
  preserve the existing ``DEBUG:root:``, ``WARNING:root:``, and
  ``ERROR:root:`` output, including diagnostics emitted by package helpers.
* ``PERF203`` (exception handler in a loop). Retry handling needs to catch
  an exception for each individual attempt, and best-effort resource
  loading needs to continue after an individual resource fails. Their
  costs are dominated by I/O, and the handler belongs with the operation
  it governs.
* ``PLC0415`` (local import). Parquet engine resolution imports
  ``pandas.io.parquet.get_engine`` when Parquet output is requested.
  This keeps the dependency on that pandas implementation detail out of
  package startup.
* ``PLR0912`` (branch count). Functions such as ``R3Resource.load`` and
  ``available_samples`` dispatch on file formats, source, and available
  metadata. Keeping dispatch branches together makes the priority and
  fallbacks visible.
* ``PLR0913`` (argument count). Public discovery and experiment-building
  methods expose independent choices such as organism, project, annotation,
  genomic unit, and join policy. Most of the wider signatures use named
  arguments, which keep those controls visible at the call site.
* ``PLR0915`` (statement count). CLI parser setup and experiment assembly
  keep related configuration, validation, and construction steps in
  execution order. A fixed statement threshold would split those sequences
  even when the helper would have no other caller.
* ``PLR2004`` (literal comparisons). Small values can directly express
  data invariants. For example, a count matrix has two dimensions and a
  paired-end read-length ratio of two has a domain-specific meaning
  explained in the surrounding code.
* ``SIM117`` (nested context managers). Separate blocks make the resource
  lifetimes clear when one context depends on another. For example,
  the shared cache acquires the R-compatible lock before accessing the SQLite
  registry, and archive copying opens a member through an open archive.
* ``SLF001`` (private member access). Package modules share internal
  implementation details, such as ``R3Resource._cached_path`` and helpers
  in ``_utils``. Those calls let the bundle and resource layers share
  cache and parsing behavior while keeping those helpers internal.
* ``TC002`` (third-party type-only imports). NumPy is a required
  dependency. Keeping its name available at runtime lets introspection
  resolve BigWig method annotations that refer to it.
* ``TC003`` (standard library type-only imports). Names such as
  ``collections.abc.Mapping`` and ``types.TracebackType`` are available
  when Sphinx or other runtime consumers resolve annotations.
* ``TRY003`` (message in a raise statement). Use the established built-in
  and recount3 exception types with context at the failure site. Many
  different validation failures use the same exception class, and the
  message tells the caller which input or resource failed.
* ``TRY004`` (use of ``ValueError`` for a type check). The count-scaling
  API raises ``ValueError`` for invalid values inside otherwise valid
  metadata containers. That includes non-numeric entries in the
  mapped-read-length column. Callers can handle malformed metadata with a
  single exception type.
* ``TRY400`` (error logging in an exception handler). The CLI reports
  expected user-facing failures such as invalid manifests, missing extras,
  and output errors with a short message. ``logging.exception`` would
  attach a traceback to those messages.

Local exceptions in package code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following rules are enabled. Their existing exceptions use a targeted
``noqa`` at the site that needs it.

* ``A002``. ``BigWigFile.stats`` exposes the ``type=`` keyword used by
  pyBigWig and by existing callers.
* ``ARG001``. ``_coerce_gtf_bp_length`` accepts coordinate keywords for
  call compatibility. It derives covered length from the GTF score
  because a feature's genomic span can differ from its covered bases.
* ``DTZ005``. Generated manifest filenames use local wall-clock time.
* ``PD101``. Metadata reconciliation counts distinct partner IDs within
  each sample group. A comparison against one constant value across the
  series would not detect the same conflicts.
* ``PLR0911``. The CLI dispatcher returns each handler's exit status
  directly from its command branch.
* ``PLR0917``. ``compute_scale_factors`` retains its public positional
  signature for existing callers.
* ``PYI034``. ``R3ResourceDescription.__new__`` is a factory that can
  return any registered concrete description. Its return type describes
  that factory behavior and works on Python 3.10.
* ``RUF022``. ``recount3.__all__`` is grouped by API role, then by name
  length within each group as described above.
* ``S110``. Reading a junction RR table can fall back to parsing the
  cached file when the resource's generic loader fails.
* ``S112``. ``R3ResourceBundle.iter_loaded`` yields what is available
  and skips a resource that fails to autoload.
* ``S310``. The download helper intentionally opens ``http://`` and
  ``https://`` mirror URLs as well as ``file://`` URLs for local mirrors.
* ``TRY301``. BigWig validation raises inside the handle's cleanup block
  so a wrong-format file closes the native handle before the error
  propagates.

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

``CHANGELOG.rst`` follows Keep a Changelog by grouping entries by version and
type of change. In reStructuredText, underline version headings with hyphens
and category headings with tildes. Use dashes for bullets.

Add an entry under ``Unreleased`` in the appropriate category, such as
``Added``, ``Changed``, ``Fixed``, or ``Removed``. This project also uses
categories such as ``Breaking changes`` and ``Internal``. State what changed
and why it matters to a user, and note explicitly whether existing code
remains source-compatible.

List released versions newest first with headings in the form
``<version> (<YYYY-MM-DD>)``. Leave released sections as published.


TODO
~~~~

Follow the `Changelog`_ structure, category conventions, and markup in
``TODO.rst``. Use forward-looking wording wherever the changelog uses past
tense. For example, use ``Add`` where the changelog uses ``Added``. Begin each
TODO entry with an imperative verb.

Use undated target versions such as ``1.5.0`` in the TODO and list the nearest
release first. Keep items without a target release under ``Unscheduled`` at
the end.


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
