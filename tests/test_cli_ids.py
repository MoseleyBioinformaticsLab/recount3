"""Tests for the `ids` subcommand writing to files and stdout."""

from __future__ import annotations

from recount3.search import create_sample_project_lists
from .helpers import call_cli_inprocess


def test_ids_stdout(monkeypatch) -> None:
    """`ids` should print two lists when no output files are passed."""
    monkeypatch.setattr(
        "recount3.search.create_sample_project_lists",
        lambda organism="": (["S1", "S2"], ["P1"]),
        raising=True,
    )
    code, out, _err = call_cli_inprocess(["ids", "--organism", "human"])
    # Both lists should be interleaved on stdout (samples then projects).
    assert code == 0
    assert "S1" in out and "S2" in out and "P1" in out


def test_ids_files(monkeypatch, tmp_path) -> None:
    """`ids` writes to the requested files and logs success."""
    monkeypatch.setattr(
        "recount3.search.create_sample_project_lists",
        lambda organism="": (["S1"], ["P1", "P2"]),
        raising=True,
    )
    sf = tmp_path / "samples.txt"
    pf = tmp_path / "projects.txt"
    code, _out, _err = call_cli_inprocess(
        ["ids", "--organism", "mouse", "--samples-out", str(sf), "--projects-out", str(pf)]
    )
    assert code == 0
    assert sf.read_text(encoding="utf-8").strip().splitlines() == ["S1"]
    assert pf.read_text(encoding="utf-8").strip().splitlines() == ["P1", "P2"]
