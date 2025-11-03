"""Integration-ish tests for `recount3 search` across all modes.

Validates:
  * Required filter handling per mode.
  * Output formats: JSONL and TSV.
  * Destinations: explicit `--output` and timestamped `--outdir`.
  * Column ordering for TSV output is stable.

All `recount3.search.*` functions are monkeypatched to return tiny, concrete
`R3Resource` objects; there is no network I/O.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Callable, Iterable

import pytest

from recount3 import search as r3_search  # module namespace for monkeypatching
from .helpers import call_cli_inprocess, make_resource


def _annotations_resource():
    return make_resource(
        resource_type="annotations",
        organism="human",
        genomic_unit="gene",
        annotation_file_extension="G026",
    )


def _gene_exon_resource():
    return make_resource(
        resource_type="count_files_gene_or_exon",
        organism="human",
        data_source="sra",
        genomic_unit="gene",
        project="SRP012345",
        annotation_file_extension="G026",
    )


def _junctions_resource():
    return make_resource(
        resource_type="count_files_junctions",
        organism="human",
        data_source="sra",
        project="SRP012345",
        junction_type="ALL",
        junction_file_extension="MM",
    )


def _metadata_resource():
    return make_resource(
        resource_type="metadata_files",
        organism="human",
        data_source="sra",
        table_name="recount_pred",
        project="SRP096765",
    )


def _bigwig_resource():
    return make_resource(
        resource_type="bigwig_files",
        organism="mouse",
        data_source="sra",
        project="DRP001299",
        sample="DRR014697",
    )


def _sources_resource():
    # Data-sources are catalog-like; one per data_source.
    return make_resource(resource_type="data_sources", organism="human")


def _source_meta_resource():
    return make_resource(
        resource_type="data_source_metadata",
        organism="human",
        data_source="sra",
    )


@pytest.mark.parametrize(
    ("mode", "filters", "patch_target", "factory"),
    [
        ("annotations", ["organism=human", "genomic_unit=gene", "annotation_file_extension=G026"],
         "search_annotations", _annotations_resource),
        ("gene-exon", ["organism=human", "data_source=sra", "genomic_unit=gene", "project=SRP012345"],
         "search_count_files_gene_or_exon", _gene_exon_resource),
        ("junctions", ["organism=human", "data_source=sra", "project=SRP012345"],
         "search_count_files_junctions", _junctions_resource),
        ("metadata", ["organism=human", "data_source=sra", "table_name=recount_pred", "project=SRP096765"],
         "search_metadata_files", _metadata_resource),
        ("bigwig", ["organism=mouse", "data_source=sra", "project=DRP001299", "sample=DRR014697"],
         "search_bigwig_files", _bigwig_resource),
        ("sources", ["organism=human"], "search_data_sources", _sources_resource),
        ("source-meta", ["organism=human", "data_source=sra"],
         "search_data_source_metadata", _source_meta_resource),
    ],
)
def test_search_jsonl_and_tsv(
    monkeypatch,
    tmp_path: Path,
    mode: str,
    filters: list[str],
    patch_target: str,
    factory: Callable[[], object],
) -> None:
    """Each mode should emit records in JSONL and TSV and honor destinations."""
    # Return two fake resources to verify multi-line output.
    monkeypatch.setattr(r3_search, patch_target, lambda **_: [factory(), factory()], raising=True)

    # JSONL to --output
    out_file = tmp_path / "out.jsonl"
    code, out, err = call_cli_inprocess(["search", mode, *filters, "--format=jsonl", "--output", str(out_file)])
    assert code == 0, err
    text = out_file.read_text(encoding="utf-8").strip().splitlines()
    assert len(text) == 2
    for line in text:
        obj = json.loads(line)
        assert obj["resource_type"]
        assert "url" in obj and "arcname" in obj  # convenience keys

    # TSV to --outdir (timestamped filename)
    code, out, err = call_cli_inprocess(["search", mode, *filters, "--format=tsv", "--outdir", str(tmp_path)])
    assert code == 0, err
    files = list(tmp_path.glob(f"{mode}-*.tsv"))
    assert len(files) == 1
    tsv = files[0].read_text(encoding="utf-8").splitlines()
    assert tsv[0].split("\t")[:5] == [
        "resource_type",
        "organism",
        "data_source",
        "genomic_unit",
        "project",
    ]
    assert len(tsv) == 3  # header + 2 rows

    # When neither --output nor --outdir is given, results go to stdout.
    code, out, err = call_cli_inprocess(["search", mode, *filters, "--format=jsonl"])
    assert code == 0 and out.count("\n") == 2
