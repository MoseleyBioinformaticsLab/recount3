"""CLI tests for `bundle stack-counts` writing TSV and Parquet.

Avoid pandas/pyarrow by returning a tiny object with `.to_csv` and
`.to_parquet` methods from the mocked `R3ResourceBundle.stack_count_matrices`.
"""

from __future__ import annotations

from pathlib import Path

import recount3.cli as cli
from .helpers import call_cli_inprocess, write_jsonl


class _FakeDF:
    """A minimal 'dataframe-like' object with pandas' writer methods."""

    def __init__(self) -> None:
        self.written_to: list[Path] = []

    def to_csv(self, path: Path, sep: str = "\t") -> None:
        path = Path(path)
        path.write_text("c1{}c2\n1{}2\n".format(sep, sep), encoding="utf-8")
        self.written_to.append(path)

    def to_parquet(self, path: Path) -> None:
        path = Path(path)
        path.write_bytes(b"PAR1fake")  # trivial sentinel
        self.written_to.append(path)


class _FakeBundle:
    """Bundle stub used to bypass heavy dependencies and network I/O."""

    def __init__(self) -> None:  # pragma: no cover - trivial
        self._added = 0

    def extend(self, _resources) -> None:  # pragma: no cover - trivial
        self._added += len(list(_resources))

    def stack_count_matrices(self, **_kw):
        return _FakeDF()


def _manifest_row() -> dict:
    return dict(
        resource_type="count_files_gene_or_exon",
        organism="human",
        data_source="sra",
        genomic_unit="gene",
        project="SRP012345",
        annotation_file_extension="G026",
    )


def test_bundle_writes_tsv_and_parquet(monkeypatch, tmp_path: Path) -> None:
    """Stacking writes the chosen file format without importing pandas."""
    # Replace the bundle class used by the CLI.
    monkeypatch.setattr(cli, "R3ResourceBundle", _FakeBundle, raising=True)

    man = tmp_path / "m.jsonl"
    write_jsonl([_manifest_row(), _manifest_row()], man)

    # TSV
    tsv_out = tmp_path / "out.tsv"
    code, _out, err = call_cli_inprocess([
        "bundle", "stack-counts",
        "--from", str(man),
        "--out", str(tsv_out),
        "--compat", "family",
        "--join", "inner",
        "--axis", "1",
    ])
    assert code == 0, err
    assert tsv_out.exists() and "c1\tc2" in tsv_out.read_text(encoding="utf-8")

    # Parquet (fake)
    pq_out = tmp_path / "out.parquet"
    code, _out, err = call_cli_inprocess([
        "bundle", "stack-counts",
        "--from", str(man),
        "--out", str(pq_out),
    ])
    assert code == 0, err
    assert pq_out.exists() and pq_out.read_bytes().startswith(b"PAR1")
