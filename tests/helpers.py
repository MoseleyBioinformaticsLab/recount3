"""Utility helpers for CLI testing.

This module centralizes test helpers so CLI tests remain concise, deterministic,
and Google-style. Nothing here performs network I/O.

Functions:
  call_cli_inprocess: Run the CLI `main()` in-process and capture stdio.
  call_cli_subprocess: Run `python -m recount3` in a subprocess safely.
  make_resource: Create a minimal R3Resource for mocking search results.
  write_jsonl: Write JSONL lines to a file atomically.
"""

from __future__ import annotations

import io
import json
import runpy
import sys
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path
from typing import Iterable, Sequence

from recount3._descriptions import R3ResourceDescription
from recount3.config import default_config
from recount3.resource import R3Resource


def call_cli_inprocess(argv: Sequence[str]) -> tuple[int, str, str]:
    """Run `recount3.cli.main(argv)` in-process and capture stdio.

    Args:
      argv: Argument vector excluding the program name, e.g. `["--version"]`.

    Returns:
      (exit_code, stdout_text, stderr_text).
    """
    # Import the module fresh to ensure the latest code is used.
    import recount3.cli as cli

    out = io.StringIO()
    err = io.StringIO()
    code = 0
    with redirect_stdout(out), redirect_stderr(err):
        try:
            cli.main(list(argv))  # main() will sys.exit(code)
        except SystemExit as exc:  # expected control path
            code = int(exc.code or 0)
    return code, out.getvalue(), err.getvalue()


def call_cli_subprocess(argv: Sequence[str], input_text: str | None = None):
    """Run `python -m recount3` in a new interpreter process.

    This is closer to real user behavior and lets us pipe data to stdin.

    Args:
      argv: Argument vector for the module entrypoint.
      input_text: Optional text to feed to stdin. If provided, stdin is piped.

    Returns:
      CompletedProcess with `returncode`, `stdout`, `stderr` (all text).
    """
    import subprocess

    cmd = [sys.executable, "-m", "recount3", *argv]
    return subprocess.run(
        cmd,
        input=input_text,
        capture_output=True,
        text=True,
        check=False,
    )


def make_resource(**fields) -> R3Resource:
    """Construct a minimal `R3Resource` for use in mocked search results.

    The `fields` are fed to `R3ResourceDescription(**fields)`. Provide a valid
    combination per resource type (tests use simple “annotations” defaults).

    Returns:
      A concrete `R3Resource` with a default config, safe for URL/arcname use.
    """
    desc = R3ResourceDescription(**fields)
    return R3Resource(description=desc, config=default_config())


def write_jsonl(objs: Iterable[dict], path: Path) -> None:
    """Write JSONL objects to `path` atomically.

    Args:
      objs: Iterable of JSON-serializable mappings; one per line.
      path: Destination file path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as fh:
        for obj in objs:
            fh.write(json.dumps(obj, ensure_ascii=False) + "\n")
    tmp.replace(path)
