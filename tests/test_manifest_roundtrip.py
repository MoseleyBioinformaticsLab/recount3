"""Tests for manifest JSONL round-trip and object rehydration.

Verify that:
  * `_iter_manifest('-')` reads from stdin.
  * `_resource_from_dict` ignores convenience keys (`url`, `arcname`).
"""

from __future__ import annotations

import json
from typing import Iterable

from .helpers import call_cli_subprocess, make_resource


def _jsonl_lines() -> str:
    res = make_resource(
        resource_type="annotations",
        organism="human",
        genomic_unit="gene",
        annotation_file_extension="G026",
    )
    # One well-formed manifest object including convenience keys.
    d = res.description.__dict__.copy()
    d.update({"url": res.url, "arcname": res.arcname})
    return json.dumps(d) + "\n" + json.dumps(d) + "\n"


def test_iter_manifest_from_stdin_round_trip() -> None:
    """`download --from=-` should ingest JSONL from stdin without errors."""
    # Don't need to actually download; pass `--dest` to a tmp zip path so
    # dir creation is exercised while downloads are mocked as skipped by design.
    proc = call_cli_subprocess(["download", "--from=-", "--dest", "out.zip"], input_text=_jsonl_lines())
    # All entries will "ok" or "skipped"; the purpose is decoding & rehydrate.
    assert proc.returncode in (0, 3)
    # Two JSONL events should be printed.
    assert len([l for l in proc.stdout.strip().splitlines() if l.strip()]) == 2
