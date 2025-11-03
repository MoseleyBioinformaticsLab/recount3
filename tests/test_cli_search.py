"""CLI tests for `recount3 search`.

These tests validate:
  * Required-filter enforcement per mode.
  * JSONL and TSV emission (structure and stable columns).
  * Destination selection logic (`--output` vs. `--outdir`) and timestamped
    filenames.
  * Logging messages for stdout vs file output.

All network access is mocked by default via the `no_network` autouse fixture.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Any

import pytest

from recount3._descriptions import R3ResourceDescription
from recount3.resource import R3Resource

from .helpers import call_cli_inprocess


def _fake_resources_for_annotations(n: int = 1) -> list[R3Resource]:
    """Create `n` fake annotation resources with deterministic fields.

    Args:
      n: Number of resources to yield.

    Returns:
      A list of :class:`R3Resource` objects with minimal fields set.
    """
    res: list[R3Resource] = []
    for i in range(n):
        desc = R3ResourceDescription(
            resource_type="annotations",
            organism="human",
            genomic_unit="gene",
            annotation_file_extension="G026",
        )
        # Note: config is attached later by the CLI; a None/placeholder is fine.
        res.append(R3Resource(description=desc))  # type: ignore[call-arg]
    return res


@pytest.mark.parametrize(
    ("mode", "missing"),
    [
        ("annotations", "annotation_file_extension"),
        ("gene-exon", "project"),
        ("junctions", "project"),
        ("metadata", "table_name"),
        ("bigwig", "sample"),
        ("sources", "organism"),
        ("source-meta", "data_source"),
    ],
)
def test_search_required_filters_enforced(mode: str, missing: str) -> None:
    """Missing required filters should cause a usage error handled by CLI."""
    base_args = {
        "annotations": ["organism=human", "genomic_unit=gene"],
        "gene-exon": ["organism=human", "data_source=sra", "genomic_unit=gene"],
        "junctions": ["organism=human", "data_source=sra"],
        "metadata": ["organism=human", "data_source=sra", "project=SRP1"],
        "bigwig": ["organism=human", "data_source=sra", "project=SRP1"],
        "sources": [],
        "source-meta": ["organism=human"],
    }[mode]

    code, _out, err = call_cli_inprocess(["search", mode, *base_args, "--format=jsonl"])
    # argparse & main() convert validation exceptions to exit code 2.
    assert code == 2
    assert "Missing required filters" in err


def test_search_annotations_emits_jsonl(monkeypatch) -> None:
    """`search annotations` prints JSONL with augmented url/arcname keys."""
    # Monkeypatch the search function to return deterministic resources.
    import recount3.search as r3_search

    monkeypatch.setattr(
        r3_search, "search_annotations", lambda **_: _fake_resources_for_annotations(2)
    )

    code, out, err = call_cli_inprocess(
        [
            "search",
            "annotations",
            "organism=human",
            "genomic_unit=gene",
            "annotation_file_extension=G026",
            "--format=jsonl",
        ]
    )
    assert code == 0
    lines = [ln for ln in out.splitlines() if ln.strip()]
    assert len(lines) == 2, "One JSON object per resource expected."
    # Each line should contain description fields plus convenience keys.
    assert '"resource_type":"annotations"' in lines[0]
    assert '"url":' in lines[0]
    assert '"arcname":' in lines[0]
    # Logging goes to stderr; ensure the summary is present.
    assert "Emitted 2 resources" in err


def test_search_annotations_emits_tsv(monkeypatch, tmp_path) -> None:
    """`search` can emit TSV with stable column order to a file via --output."""
    import recount3.search as r3_search

    monkeypatch.setattr(
        r3_search, "search_annotations", lambda **_: _fake_resources_for_annotations(1)
    )
    out_file = tmp_path / "manifest.tsv"

    code, _out, err = call_cli_inprocess(
        [
            "search",
            "annotations",
            "organism=human",
            "genomic_unit=gene",
            "annotation_file_extension=G026",
            "--format=tsv",
            "--output",
            str(out_file),
        ]
    )
    assert code == 0
    text = out_file.read_text(encoding="utf-8").splitlines()
    assert text[0].startswith(
        "resource_type\torganism\tdata_source\tgenomic_unit"
    ), "Stable TSV header expected."
    assert len(text) == 2, "Header + one row expected."
    assert "Emitted 1 resources to" in err


def test_search_outdir_timestamped_filename(monkeypatch, tmp_path) -> None:
    """When `--outdir` is provided, a timestamped file name is constructed."""
    import recount3.search as r3_search
    import recount3.cli as cli_mod

    monkeypatch.setattr(
        r3_search, "search_annotations", lambda **_: _fake_resources_for_annotations(1)
    )
    # Patch datetime.now in the CLI module to keep the filename deterministic.
    class _FrozenDT:
        @staticmethod
        def now():
            from datetime import datetime

            return datetime(2025, 1, 2, 3, 4, 5)

        @staticmethod
        def strptime(*_a, **_k):
            raise NotImplementedError

    monkeypatch.setattr(cli_mod, "datetime", _FrozenDT, raising=True)
    outdir = tmp_path / "manifests"

    code, _out, err = call_cli_inprocess(
        [
            "search",
            "annotations",
            "organism=human",
            "genomic_unit=gene",
            "annotation_file_extension=G026",
            "--format=jsonl",
            "--outdir",
            str(outdir),
        ]
    )
    assert code == 0
    # The CLI logs the full path containing the deterministic timestamp.
    assert "annotations-20250102-030405.jsonl" in err
    # File exists:
    files = list(outdir.glob("annotations-*.jsonl"))
    assert len(files) == 1
