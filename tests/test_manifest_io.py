"""Unit tests for manifest input/output helpers within the CLI.

Validates the following:
  * JSONL writer adds `url` and `arcname` convenience keys.
  * TSV writer uses a stable, documented column order and blanks for missing
    fields.
  * Manifest rehydration (`_iter_manifest`) strips convenience keys and
    reconstructs :class:`R3Resource` objects with equivalent descriptions.

These tests import the underscore-prefixed helpers from `recount3.cli` since
the CLI is the canonical implementation of manifest I/O for this package.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import pytest

import recount3.cli as cli
from recount3._descriptions import R3ResourceDescription
from recount3.resource import R3Resource


def _mk(res_type: str, **kwargs) -> R3Resource:
    """Small helper to build a resource with a description only."""
    desc = R3ResourceDescription(resource_type=res_type, **kwargs)
    return R3Resource(description=desc)  # type: ignore[call-arg]


def test_write_jsonl_and_iter_manifest_roundtrip(tmp_path) -> None:
    """JSONL writer should be readable by the manifest iterator."""
    res = [
        _mk("annotations", organism="human", genomic_unit="gene", annotation_file_extension="G026"),
        _mk(
            "metadata_files",
            organism="human",
            data_source="sra",
            project="SRP1",
            table_name="recount_pred",
        ),
    ]
    out = tmp_path / "m.jsonl"
    cli._write_jsonl(res, out)

    text = out.read_text(encoding="utf-8").splitlines()
    assert '"url":' in text[0] and '"arcname":' in text[0], "Convenience keys expected."

    cfg = cli._build_config_from_env_and_flags(
        type("Args", (), dict(
            base_url=None, cache_dir=None, timeout=None, retries=None,
            insecure_ssl=False, user_agent=None, chunk_size=None,
        ))()  # type: ignore
    )
    roundtrip = list(cli._iter_manifest(str(out), cfg))
    assert len(roundtrip) == 2
    # Descriptions should match the originals (url/arcname were ignored).
    for a, b in zip(res, roundtrip):
        assert dataclasses.asdict(a.description) == dataclasses.asdict(b.description)  # type: ignore


def test_write_tsv_header_and_row_counts(tmp_path) -> None:
    """TSV writer emits a stable header and as many rows as resources."""
    res = [
        _mk("data_sources", organism="human"),
        _mk("data_source_metadata", organism="human", data_source="sra"),
    ]
    out = tmp_path / "m.tsv"
    cli._write_tsv(res, out)
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith(
        "resource_type\torganism\tdata_source\tgenomic_unit"
    ), "Stable header not found."
    assert len(lines) == 1 + len(res)
