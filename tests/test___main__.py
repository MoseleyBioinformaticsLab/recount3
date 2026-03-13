from __future__ import annotations

import os
import pathlib
import subprocess
import sys

import recount3.__main__ as main_module
import recount3.cli


def test_main_module_imports_without_error() -> None:
    """recount3.__main__ must be importable without raising."""
    assert main_module is not None


def test_main_delegates_to_cli() -> None:
    """The 'main' name in __main__ must be the same object as recount3.cli.main."""
    assert main_module.main is recount3.cli.main


def test_python_m_recount3_help() -> None:
    """python -m recount3 --help must exit 0 and print usage information."""
    src = str(pathlib.Path(__file__).parent.parent / "src")
    env = os.environ.copy()
    env["PYTHONPATH"] = src + os.pathsep + env.get("PYTHONPATH", "")

    result = subprocess.run(
        [sys.executable, "-m", "recount3", "--help"],
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode == 0
    assert "recount3" in result.stdout.lower()
