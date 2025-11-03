"""Logging behavior and verbosity flags."""

from __future__ import annotations

from .helpers import call_cli_inprocess


def test_quiet_suppresses_info(monkeypatch) -> None:
    """`--quiet` should suppress the standard INFO summary line."""
    # Minimal controlled search result.
    monkeypatch.setattr(
        "recount3.search.search_annotations",
        lambda **_: [],
        raising=True,
    )
    code, _out, err = call_cli_inprocess(
        [
            "--quiet",
            "search",
            "annotations",
            "organism=human",
            "genomic_unit=gene",
            "annotation_file_extension=G026",
            "--format=jsonl",
        ]
    )
    assert code == 0
    # With zero results there may be no summary message; ensure no INFO noise.
    assert "INFO:" not in err
