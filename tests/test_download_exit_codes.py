"""Exit-code behavior for the `download` command under varying outcomes.

Simulates:
  * 0 failures  -> exit 0
  * 1 of 3 fails -> exit 3
  * 3 of 3 fail -> exit 2

No network I/O occurs because `R3Resource.download` is monkeypatched.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from recount3.resource import R3Resource
from .helpers import call_cli_inprocess, write_jsonl


def _manifest_lines() -> list[dict]:
    base = dict(
        resource_type="annotations",
        organism="human",
        genomic_unit="gene",
        annotation_file_extension="G026",
    )
    return [dict(base), dict(base), dict(base)]


def _run(monkeypatch, tmp_path: Path, outcomes: list[str]) -> tuple[int, list[dict]]:
    # outcomes: "ok" or "error" per call order
    calls = {"i": 0}

    def fake_download(self: R3Resource, **_kw):
        i = calls["i"]
        calls["i"] = i + 1
        if outcomes[i] == "error":
            raise RuntimeError("boom")
        return str(tmp_path / f"ok_{i}.dat")

    monkeypatch.setattr(R3Resource, "download", fake_download, raising=True)

    man = tmp_path / "m.jsonl"
    write_jsonl(_manifest_lines(), man)
    code, out, _err = call_cli_inprocess(["download", "--from", str(man), "--dest", str(tmp_path / "dest")])
    events = [json.loads(line) for line in out.strip().splitlines() if line.strip()]
    return code, events


def test_download_all_ok(monkeypatch, tmp_path: Path) -> None:
    """All successful -> exit 0."""
    code, events = _run(monkeypatch, tmp_path, ["ok", "ok", "ok"])
    assert code == 0
    assert all(e["status"] == "ok" for e in events)


def test_download_partial(monkeypatch, tmp_path: Path) -> None:
    """Partial failures -> exit 3."""
    code, events = _run(monkeypatch, tmp_path, ["ok", "error", "ok"])
    assert code == 3
    assert any(e["status"] == "error" for e in events)


def test_download_all_fail(monkeypatch, tmp_path: Path) -> None:
    """All failures -> exit 2."""
    code, events = _run(monkeypatch, tmp_path, ["error", "error", "error"])
    assert code == 2
    assert all(e["status"] == "error" for e in events)
