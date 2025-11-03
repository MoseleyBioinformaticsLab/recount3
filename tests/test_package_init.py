"""Basic import checks for the installed package."""

from __future__ import annotations

import importlib


def test_import_recount3_top_level() -> None:
    """Top-level import should succeed and expose a version."""
    pkg = importlib.import_module("recount3")
    assert hasattr(pkg, "__version__") or hasattr(pkg, "__all__")
