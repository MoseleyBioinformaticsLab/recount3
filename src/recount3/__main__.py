"""Module entry point.

Delegates to :mod:`recount3.cli`.
"""

from __future__ import annotations

from .cli import main

if __name__ == "__main__":
    main()
 
"""
TODO: 

tests/
  unit/
    test_descriptions.py
    test_resource.py
    test_bigwig.py
    test_search.py
    test_bundle.py
    test_config.py
  integration/
    test_cli_smoketest.py
    test_download_materialize.py
docs/
  api_overview.md
pyproject.toml
pylintrc
"""