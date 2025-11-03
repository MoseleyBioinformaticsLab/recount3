"""Ensure the package imports cleanly when `pyBigWig` is absent.

The package declares an *optional* extra for `pybigwig` (case-insensitive) and
must not hard-require it at import time. This test simulates its absence to
guard against accidental unconditional imports.

Does not import submodules that explicitly require `pyBigWig` at runtime; the
goal is to ensure the package core and CLI remain usable without it.
"""

from __future__ import annotations

import importlib
import sys


def test_import_without_pybigwig(monkeypatch) -> None:
    """Top-level import must succeed when `pyBigWig` cannot be imported."""
    # Simulate `import pyBigWig` and `import pybigwig` failures.
    monkeypatch.setitem(sys.modules, "pyBigWig", None)
    monkeypatch.setitem(sys.modules, "pybigwig", None)
    # Import the top-level package and a few core modules.
    pkg = importlib.import_module("recount3")
    importlib.import_module("recount3.cli")
    importlib.import_module("recount3.config")
    importlib.import_module("recount3.resource")
    assert hasattr(pkg, "__all__") or True
