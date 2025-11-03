"""Smoke tests for top-level help and version behaviors.

Validates that:
  * The module entrypoint prints a usage summary with `-h/--help`.
  * `--version` exits with code 0 and prints a SemVer-like version string.
"""

from __future__ import annotations

import re

from .helpers import call_cli_inprocess, call_cli_subprocess


def test_cli_help_inprocess() -> None:
    """`-h` should print usage text and exit cleanly."""
    code, out, _err = call_cli_inprocess(["-h"])
    assert code == 0
    assert "usage:" in out.lower()
    assert "recount3" in out


def test_cli_version_inprocess() -> None:
    """`--version` prints a SemVer-like version and exits 0."""
    code, out, _err = call_cli_inprocess(["--version"])
    assert code == 0
    assert re.search(r"\d+\.\d+\.\d+", out) is not None


def test_cli_module_entrypoint_subprocess() -> None:
    """`python -m recount3 --version` should also work (e2e entrypoint)."""
    proc = call_cli_subprocess(["--version"])
    assert proc.returncode == 0
    assert re.search(r"\d+\.\d+\.\d+", proc.stdout) is not None
