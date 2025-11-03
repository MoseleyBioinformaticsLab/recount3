"""CLI tests for `recount3 bundle stack-counts`.

Mocks the heavy frame-building work so the CLI logic (argument parsing,
compatibility flags, writing output) is exercised without requiring real data
files or remote access.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from .helpers import call_cli_inprocess, write_jsonl


def test_bundle_stack_counts_writes_tsv(monkeypatch, tmp_path) -> None:
    """Bundle writes a TSV file and logs success."""
    # Build a tiny manifest with one item (content is irrelevant for this test).
    man = tmp_path / "counts.jsonl"
    write_jsonl([{"resource_type": "count_files_gene_or_exon",
                  "organism": "human", "data_source": "sra",
                  "genomic_unit": "gene", "project": "SRP1",
                  "annotation_file_extension": "G026"}], man)

    # Patch the bundle to avoid pulling real data. Return a 2x2 DataFrame.
    import pandas as pd
    import recount3.bundle as bundle_mod

    class _FakeBundle:
        def __init__(self) -> None:
            pass

        def extend(self, _resources) -> None:
            return None

        def stack_count_matrices(self, **_kwargs):
            return pd.DataFrame({"A": [1, 2]}, index=["g1", "g2"])

    monkeypatch.setattr(bundle_mod, "R3ResourceBundle", _FakeBundle, raising=True)

    out = tmp_path / "out.tsv"
    code, _out, err = call_cli_inprocess(
        [
            "bundle",
            "stack-counts",
            "--from",
            str(man),
            "--out",
            str(out),
            "--compat=family",
            "--join=inner",
            "--axis=1",
        ]
    )
    assert code == 0
    assert out.exists()
    assert "Wrote stacked table to:" in err
