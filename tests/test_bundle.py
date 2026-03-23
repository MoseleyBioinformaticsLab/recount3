# pylint: disable=redefined-outer-name,protected-access

from __future__ import annotations

import gzip
import io
import logging
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

import recount3._utils as _utils_module
import recount3.bundle as bmod
import recount3.search as search_module
from recount3.bundle import (
    R3ResourceBundle,
    _align_ranges_to_features,
    _choose_alignment_key,
    _coerce_gtf_bp_length,
    _coerce_gtf_phase_column,
    _collapse_rows_by_key,
    _construct_ranged_summarized_experiment,
    _construct_summarized_experiment,
    _count_compat_keys,
    _dedupe_ranges_on_feature_id,
    _ensure_unique_columns,
    _make_unique_names,
    _metadata_origin,
    _namespace_metadata_columns,
    _outer_merge_metadata_frames,
    _parse_gtf_attributes,
    _maybe_relabel_counts_columns_to_external_id,
    _peek_gtf_feature_counts,
    _ranges_from_gtf,
    _read_gtf_dataframe,
    _read_rr_table,
    _select_gtf_resource_for_unit,
    _standardize_metadata_frame,
    _strip_ensembl_version,
    _to_genomic_ranges,
)
from recount3.config import Config
from recount3.errors import CompatibilityError
from recount3.resource import R3Resource

_TESTS_DIR = Path(__file__).parent
_DATA_DIR = _TESTS_DIR / "data"
_MIRROR = _DATA_DIR / "recount3_mirror" / "recount3"

_GENE_GTF_GZ = (
    _MIRROR / "human" / "annotations" / "gene_sums" / "human.gene_sums.G026.gtf.gz"
)
_EXON_GTF_GZ = (
    _MIRROR / "human" / "annotations" / "exon_sums" / "human.exon_sums.G026.gtf.gz"
)



@pytest.fixture()
def local_config(tmp_path: Path) -> Config:
    """Return a Config whose cache_dir is tmp_path and base_url points at the
    local test-data mirror."""
    return Config(
        base_url=f"file://{_MIRROR}/",
        timeout=5,
        insecure_ssl=False,
        max_retries=1,
        user_agent="test",
        cache_dir=tmp_path / "cache",
        cache_disabled=False,
    )


def _mock_resource(
    resource_type: str = "count_files_gene_or_exon",
    *,
    url: str = "http://example.com/file.gz",
    url_path_val: str = "file.gz",
    loaded_data: Any = None,
    is_loaded_val: bool = False,
    **desc_attrs: Any,
) -> MagicMock:
    """Return a MagicMock shaped like an R3Resource."""
    res = MagicMock(spec=R3Resource)
    res.url = url
    desc = MagicMock()
    desc.resource_type = resource_type
    desc.url_path.return_value = url_path_val
    for key, val in desc_attrs.items():
        setattr(desc, key, val)
    res.description = desc
    if loaded_data is not None:
        res.is_loaded.return_value = True
        res.get_loaded.return_value = loaded_data
    else:
        res.is_loaded.return_value = is_loaded_val
        res.get_loaded.return_value = None
    return res


def _gene_df(features: list[str] | None = None, samples: list[str] | None = None) -> pd.DataFrame:
    """Return a tiny gene-counts DataFrame."""
    rows = features or ["ENSG0001", "ENSG0002"]
    cols = samples or ["SRR001", "SRR002"]
    data = np.arange(len(rows) * len(cols), dtype=float).reshape(len(rows), len(cols))
    return pd.DataFrame(data, index=rows, columns=cols)


class TestEnsureUniqueColumns:
    def test_no_duplicates(self) -> None:
        df = pd.DataFrame({"a": [1], "b": [2]})
        out = _ensure_unique_columns(df)
        assert list(out.columns) == ["a", "b"]

    def test_duplicates_suffixed(self) -> None:
        df = pd.DataFrame([[1, 2, 3]], columns=["a", "b", "a"])
        out = _ensure_unique_columns(df)
        assert list(out.columns) == ["a", "b", "a__2"]

    def test_three_duplicates(self) -> None:
        df = pd.DataFrame([[1, 2, 3]], columns=["x", "x", "x"])
        out = _ensure_unique_columns(df)
        assert list(out.columns) == ["x", "x__2", "x__3"]

    def test_none_columns_replaced_by_prefix(self) -> None:
        df = pd.DataFrame([[1, 2]], columns=[None, None])
        out = _ensure_unique_columns(df)
        assert list(out.columns) == ["col", "col__2"]

    def test_none_columns_custom_prefix(self) -> None:
        df = pd.DataFrame([[1]], columns=[None])
        out = _ensure_unique_columns(df, empty_prefix="myprefix")
        assert list(out.columns) == ["myprefix"]

    def test_original_not_mutated(self) -> None:
        df = pd.DataFrame([[1, 2]], columns=["z", "z"])
        _ = _ensure_unique_columns(df)
        assert list(df.columns) == ["z", "z"]

    def test_empty_string_column_replaced(self) -> None:
        df = pd.DataFrame([[1]], columns=[""])
        out = _ensure_unique_columns(df, empty_prefix="empty")
        assert list(out.columns) == ["empty"]


class TestStandardizeMetadataFrame:
    def test_lowercases_columns(self) -> None:
        df = pd.DataFrame({"RAIL_ID": ["1"], "ExTerNal_Id": ["x"]})
        out = _standardize_metadata_frame(df)
        assert "rail_id" in out.columns
        assert "external_id" in out.columns

    def test_strips_whitespace_from_column_names(self) -> None:
        df = pd.DataFrame({" rail_id ": ["1"]})
        out = _standardize_metadata_frame(df)
        assert "rail_id" in out.columns

    def test_adds_missing_merge_keys(self) -> None:
        df = pd.DataFrame({"other": [1]})
        out = _standardize_metadata_frame(df)
        for key in ("rail_id", "external_id", "study"):
            assert key in out.columns

    def test_synonym_study_acc_renamed(self) -> None:
        df = pd.DataFrame({"study_acc": ["SRP001"]})
        out = _standardize_metadata_frame(df)
        assert "study" in out.columns
        assert "study_acc" not in out.columns

    def test_synonym_run_acc_renamed(self) -> None:
        df = pd.DataFrame({"run_acc": ["SRR001"]})
        out = _standardize_metadata_frame(df)
        assert "external_id" in out.columns

    def test_synonym_run_accession_renamed(self) -> None:
        df = pd.DataFrame({"run_accession": ["SRR001"]})
        out = _standardize_metadata_frame(df)
        assert "external_id" in out.columns

    def test_synonym_run_renamed(self) -> None:
        df = pd.DataFrame({"run": ["SRR001"]})
        out = _standardize_metadata_frame(df)
        assert "external_id" in out.columns

    def test_synonym_not_renamed_when_target_already_present(self) -> None:
        df = pd.DataFrame({"run": ["SRR001"], "external_id": ["SRR002"]})
        out = _standardize_metadata_frame(df)
        assert "external_id" in out.columns
        assert "run" in out.columns

    def test_merge_key_columns_cast_to_string(self) -> None:
        df = pd.DataFrame({"rail_id": [1], "external_id": [2], "study": [3]})
        out = _standardize_metadata_frame(df)
        assert str(out["rail_id"].dtype) == "string"

    def test_original_not_mutated(self) -> None:
        df = pd.DataFrame({"RAIL_ID": ["1"]})
        _standardize_metadata_frame(df)
        assert "RAIL_ID" in df.columns


class TestMetadataOrigin:
    def test_uses_table_name(self) -> None:
        res = _mock_resource("metadata_files", table_name="recount_qc")
        assert _metadata_origin(res) == "recount_qc"

    def test_falls_back_to_resource_type(self) -> None:
        res = _mock_resource("metadata_files")
        res.description.table_name = None
        res.description.resource_type = "custom_type"
        assert _metadata_origin(res) == "custom_type"

    def test_returns_metadata_when_no_attrs(self) -> None:
        res = _mock_resource("metadata_files")
        res.description.table_name = None
        res.description.resource_type = None
        assert _metadata_origin(res) == "metadata"

    def test_strips_and_lowercases(self) -> None:
        res = _mock_resource("metadata_files", table_name=" Recount_QC ")
        assert _metadata_origin(res) == "recount_qc"


class TestNamespaceMetadataColumns:
    def test_key_columns_not_renamed(self) -> None:
        df = pd.DataFrame({"rail_id": ["1"], "external_id": ["x"], "study": ["s"]})
        out, prov = _namespace_metadata_columns(df, origin="tbl")
        assert "rail_id" in out.columns
        assert "external_id" in out.columns
        assert "study" in out.columns
        assert not prov

    def test_non_key_columns_renamed(self) -> None:
        df = pd.DataFrame({"rail_id": ["1"], "score": [5.0]})
        out, prov = _namespace_metadata_columns(df, origin="qc")
        assert "qc__score" in out.columns
        assert "qc__score" in prov
        assert prov["qc__score"] == ("qc", "score")

    def test_custom_separator(self) -> None:
        df = pd.DataFrame({"rail_id": ["1"], "val": [0]})
        out, _ = _namespace_metadata_columns(df, origin="tbl", sep="||")
        assert "tbl||val" in out.columns


class TestOuterMergeMetadataFrames:
    def test_empty_list(self) -> None:
        result = _outer_merge_metadata_frames([])
        assert list(result.columns) == ["rail_id", "external_id", "study"]
        assert len(result) == 0

    def test_single_frame(self) -> None:
        df = pd.DataFrame({
            "rail_id": ["1"],
            "external_id": ["SRR001"],
            "study": ["SRP001"],
            "extra": ["v"],
        })
        result = _outer_merge_metadata_frames([df])
        assert len(result) == 1
        assert "extra" in result.columns

    def test_two_frames_merged(self) -> None:
        df1 = pd.DataFrame({
            "rail_id": ["1", "2"],
            "external_id": ["A", "B"],
            "study": ["S1", "S1"],
            "col_a": [10, 20],
        })
        df2 = pd.DataFrame({
            "rail_id": ["1", "3"],
            "external_id": ["A", "C"],
            "study": ["S1", "S1"],
            "col_b": [100, 300],
        })
        result = _outer_merge_metadata_frames([df1, df2])
        assert len(result) == 3  # outer merge: rows 1, 2, 3
        assert "col_a" in result.columns
        assert "col_b" in result.columns


class TestChooseAlignmentKey:
    def test_prefers_external_id_when_better_match(self) -> None:
        merged = pd.DataFrame({
            "external_id": pd.array(["SRR001", "SRR002", "SRR003"], dtype="string"),
            "rail_id": pd.array(["1", "2", "99"], dtype="string"),
        })
        key = _choose_alignment_key(
            sample_ids=["SRR001", "SRR002"],
            merged=merged,
        )
        assert key == "external_id"

    def test_prefers_rail_id_when_better_match(self) -> None:
        merged = pd.DataFrame({
            "external_id": pd.array(["XX", "YY"], dtype="string"),
            "rail_id": pd.array(["1", "2"], dtype="string"),
        })
        key = _choose_alignment_key(
            sample_ids=["1", "2"],
            merged=merged,
        )
        assert key == "rail_id"

    def test_defaults_to_external_id_when_key_absent(self) -> None:
        merged = pd.DataFrame({"other": ["a"]})
        key = _choose_alignment_key(sample_ids=["a"], merged=merged)
        assert key == "external_id"


class TestCollapseRowsByKey:
    def test_key_not_in_columns_returns_df(self) -> None:
        df = pd.DataFrame({"a": [1, 2]})
        result = _collapse_rows_by_key(df, key="missing_key")
        pd.testing.assert_frame_equal(result, df)

    def test_collapse_duplicates(self) -> None:
        df = pd.DataFrame({
            "key": ["A", "A", "B"],
            "val": [1, 2, 3],
        })
        result = _collapse_rows_by_key(df, key="key")
        assert len(result) == 2
        assert set(result["key"]) == {"A", "B"}

    def test_first_non_null_used(self) -> None:
        df = pd.DataFrame({
            "key": ["A", "A"],
            "val": [pd.NA, 42],
        })
        result = _collapse_rows_by_key(df, key="key")
        assert result.loc[result["key"] == "A", "val"].iloc[0] == 42

    def test_all_null_returns_na(self) -> None:
        df = pd.DataFrame({
            "key": ["A", "A"],
            "val": [pd.NA, pd.NA],
        })
        result = _collapse_rows_by_key(df, key="key")
        val = result.loc[result["key"] == "A", "val"].iloc[0]
        assert pd.isna(val)


class TestMaybeRelabelCountsColumnsToExternalId:
    def test_no_external_id_col(self) -> None:
        counts = _gene_df(samples=["c1", "c2"])
        col_df = pd.DataFrame({"other": ["x", "y"]}, index=["c1", "c2"])
        out_counts, out_col = _maybe_relabel_counts_columns_to_external_id(counts, col_df)
        assert list(out_counts.columns) == list(counts.columns)

    def test_length_mismatch(self) -> None:
        counts = _gene_df(samples=["c1", "c2"])
        col_df = pd.DataFrame({"external_id": ["x"]})
        out_counts, out_col = _maybe_relabel_counts_columns_to_external_id(counts, col_df)
        assert list(out_counts.columns) == list(counts.columns)

    def test_missing_external_id_values(self) -> None:
        counts = _gene_df(samples=["c1", "c2"])
        col_df = pd.DataFrame(
            {"external_id": pd.array(["SRR001", pd.NA], dtype="string")},
            index=["c1", "c2"],
        )
        out_counts, _ = _maybe_relabel_counts_columns_to_external_id(counts, col_df)
        assert list(out_counts.columns) == ["c1", "c2"]

    def test_non_unique_external_ids(self) -> None:
        counts = _gene_df(samples=["c1", "c2"])
        col_df = pd.DataFrame(
            {"external_id": pd.array(["SRR001", "SRR001"], dtype="string")},
            index=["c1", "c2"],
        )
        out_counts, _ = _maybe_relabel_counts_columns_to_external_id(counts, col_df)
        assert list(out_counts.columns) == ["c1", "c2"]

    def test_already_matching_no_rename(self) -> None:
        counts = _gene_df(samples=["SRR001", "SRR002"])
        col_df = pd.DataFrame(
            {"external_id": pd.array(["SRR001", "SRR002"], dtype="string")},
            index=["SRR001", "SRR002"],
        )
        out_counts, _ = _maybe_relabel_counts_columns_to_external_id(counts, col_df)
        assert list(out_counts.columns) == ["SRR001", "SRR002"]

    def test_successful_relabel(self) -> None:
        counts = _gene_df(samples=["1", "2"])
        col_df = pd.DataFrame(
            {"external_id": pd.array(["SRR001", "SRR002"], dtype="string")},
            index=["1", "2"],
        )
        out_counts, out_col = _maybe_relabel_counts_columns_to_external_id(counts, col_df)
        assert list(out_counts.columns) == ["SRR001", "SRR002"]
        assert list(out_col.index) == ["SRR001", "SRR002"]

    def test_empty_string_external_id_treated_as_missing(self) -> None:
        counts = _gene_df(samples=["c1", "c2"])
        col_df = pd.DataFrame(
            {"external_id": pd.array(["SRR001", ""], dtype="string")},
            index=["c1", "c2"],
        )
        out_counts, _ = _maybe_relabel_counts_columns_to_external_id(counts, col_df)
        assert list(out_counts.columns) == ["c1", "c2"]


class TestReadRrTable:
    def test_returns_dataframe_when_load_succeeds(self) -> None:
        rr_df = pd.DataFrame({"col1": [1, 2]})
        res = MagicMock(spec=R3Resource)
        res.url = "http://example.com/rr.gz"
        res.load.return_value = rr_df
        result = _read_rr_table(res)
        pd.testing.assert_frame_equal(result, rr_df)

    def test_fallback_to_file_when_load_returns_non_dataframe(
        self, tmp_path: Path
    ) -> None:
        content = "a\tb\n1\t2\n3\t4\n"
        tsv_path = tmp_path / "test.tsv"
        tsv_path.write_text(content)

        res = MagicMock(spec=R3Resource)
        res.url = "http://example.com/test.tsv"
        res.load.return_value = "not-a-dataframe"
        res._cached_path.return_value = tsv_path

        result = _read_rr_table(res)
        assert list(result.columns) == ["a", "b"]
        assert len(result) == 2

    def test_fallback_reads_gz_file(self, tmp_path: Path) -> None:
        content = b"a\tb\n10\t20\n"
        gz_path = tmp_path / "test.tsv.gz"
        with gzip.open(gz_path, "wb") as f:
            f.write(content)

        res = MagicMock(spec=R3Resource)
        res.url = "http://example.com/test.tsv.gz"
        res.load.return_value = None
        res._cached_path.return_value = gz_path

        result = _read_rr_table(res)
        assert len(result) == 1


class TestParseGtfAttributes:
    def test_empty_series_returns_empty_df(self) -> None:
        attrs = pd.Series([], dtype="object")
        result = _parse_gtf_attributes(attrs)
        assert isinstance(result, pd.DataFrame)
        assert result.empty

    def test_no_matches_returns_empty_df(self) -> None:
        attrs = pd.Series(["no kv pairs here"])
        result = _parse_gtf_attributes(attrs)
        assert isinstance(result, pd.DataFrame)

    def test_parses_key_value_pairs(self) -> None:
        attrs = pd.Series(['gene_id "ENSG001"; transcript_id "ENST001"'])
        result = _parse_gtf_attributes(attrs)
        assert "gene_id" in result.columns
        assert result["gene_id"].iloc[0] == "ENSG001"

    def test_multiple_rows(self) -> None:
        attrs = pd.Series([
            'gene_id "ENSG001"; biotype "protein_coding"',
            'gene_id "ENSG002"; biotype "lncRNA"',
        ])
        result = _parse_gtf_attributes(attrs)
        assert len(result) == 2
        assert result["gene_id"].iloc[1] == "ENSG002"


class TestCoerceGtfPhase:
    def test_valid_phases(self) -> None:
        s = pd.Series(["0", "1", "2", ".", pd.NA])
        result = _coerce_gtf_phase_column(s)
        assert result.iloc[0] == 0
        assert result.iloc[1] == 1
        assert result.iloc[2] == 2
        assert pd.isna(result.iloc[3])
        assert pd.isna(result.iloc[4])

    def test_invalid_phase_coerced_to_na_with_warning(self, caplog: pytest.LogCaptureFixture) -> None:
        s = pd.Series(["0", "5", "9"])
        with caplog.at_level(logging.WARNING):
            result = _coerce_gtf_phase_column(s)
        assert pd.isna(result.iloc[1])
        assert pd.isna(result.iloc[2])
        assert "unexpected values" in caplog.text

    def test_empty_string_becomes_na(self) -> None:
        s = pd.Series([""])
        result = _coerce_gtf_phase_column(s)
        assert pd.isna(result.iloc[0])


class TestCoerceGtfBpLength:
    def test_dot_score_returns_width(self) -> None:
        score = pd.Series([".", ".", "."])
        starts = pd.Series([1, 100, 200])
        ends = pd.Series([10, 110, 210])
        result = _coerce_gtf_bp_length(score, starts=starts, ends=ends)
        assert result.iloc[0] == 10  # 10 - 1 + 1
        assert result.iloc[1] == 11  # 110 - 100 + 1

    def test_score_matches_width_uses_score(self) -> None:
        starts = pd.Series([1, 1, 1])
        ends = pd.Series([10, 20, 30])
        widths = ends - starts + 1  # [10, 20, 30]
        score = widths.astype(str)
        result = _coerce_gtf_bp_length(score, starts=starts, ends=ends)
        assert int(result.iloc[0]) == 10

    def test_score_does_not_match_width(self) -> None:
        starts = pd.Series([1, 1])
        ends = pd.Series([10, 20])
        score = pd.Series(["999", "999"])
        result = _coerce_gtf_bp_length(score, starts=starts, ends=ends)
        assert int(result.iloc[0]) == 10  # 10 - 1 + 1
        assert int(result.iloc[1]) == 20  # 20 - 1 + 1

    def test_score_none_comparable(self) -> None:
        starts = pd.Series([5])
        ends = pd.Series([14])
        score = pd.Series(["."])
        result = _coerce_gtf_bp_length(score, starts=starts, ends=ends)
        assert int(result.iloc[0]) == 10


class TestStripEnsemblVersion:
    def test_strips_version_suffix(self) -> None:
        s = pd.Series(["ENSG00000001.12", "ENSG00000002.5"])
        result = _strip_ensembl_version(s)
        assert result.iloc[0] == "ENSG00000001"
        assert result.iloc[1] == "ENSG00000002"

    def test_no_suffix_unchanged(self) -> None:
        s = pd.Series(["ENSG00000001", "GENE_X"])
        result = _strip_ensembl_version(s)
        assert result.iloc[0] == "ENSG00000001"
        assert result.iloc[1] == "GENE_X"


def _minimal_ranges(feature_ids: list[str]) -> pd.DataFrame:
    return pd.DataFrame({
        "feature_id": feature_ids,
        "seqnames": ["chr1"] * len(feature_ids),
        "starts": list(range(1, len(feature_ids) + 1)),
        "ends": list(range(100, 100 + len(feature_ids))),
        "strand": ["+"] * len(feature_ids),
    })


class TestAlignRangesToFeatures:
    def test_missing_required_column_raises_value_error(self) -> None:
        ranges = pd.DataFrame({"feature_id": ["A"], "seqnames": ["chr1"]})
        with pytest.raises(ValueError, match="missing required columns"):
            _align_ranges_to_features(ranges, feature_ids=["A"])

    def test_exact_match_all_features(self) -> None:
        ranges = _minimal_ranges(["A", "B", "C"])
        result = _align_ranges_to_features(ranges, feature_ids=["A", "B", "C"])
        assert list(result.index) == ["A", "B", "C"]
        assert result["seqnames"].iloc[0] == "chr1"

    def test_partial_missing_uses_version_stripped_fallback(self) -> None:
        ranges = _minimal_ranges(["ENSG001.1", "ENSG002.3"])
        result = _align_ranges_to_features(
            ranges, feature_ids=["ENSG001", "ENSG002"]
        )
        assert result["seqnames"].notna().all()

    def test_conflicting_version_stripped_duplicates_raises(self) -> None:
        ranges = pd.DataFrame({
            "feature_id": ["ENSG001.1", "ENSG001.2"],
            "seqnames": ["chr1", "chr2"],
            "starts": [1, 200],
            "ends": [100, 300],
            "strand": ["+", "-"],
        })
        with pytest.raises(ValueError, match="conflicting coordinates"):
            _align_ranges_to_features(ranges, feature_ids=["ENSG001"])

    def test_fallback_still_missing_returns_partial(self) -> None:
        ranges = _minimal_ranges(["ENSG001"])
        result = _align_ranges_to_features(
            ranges, feature_ids=["ENSG001", "NOTFOUND"]
        )
        assert pd.isna(result.loc["NOTFOUND", "seqnames"])

    def test_version_strip_with_non_conflicting_duplicates(self, caplog: pytest.LogCaptureFixture) -> None:
        ranges = pd.DataFrame({
            "feature_id": ["ENSG001.1", "ENSG001.2"],
            "seqnames": ["chr1", "chr1"],
            "starts": [100, 100],
            "ends": [200, 200],
            "strand": ["+", "+"],
        })
        with caplog.at_level(logging.INFO):
            result = _align_ranges_to_features(ranges, feature_ids=["ENSG001"])
        assert result.loc["ENSG001", "seqnames"] == "chr1"


class TestReadGtfDataframe:
    def test_raises_value_error_when_cached_path_fails(self) -> None:
        res = MagicMock(spec=R3Resource)
        res.url = "http://example.com/x.gtf.gz"
        res._cached_path.side_effect = RuntimeError("no cache")
        with pytest.raises(ValueError, match="Cannot resolve cached GTF path"):
            _read_gtf_dataframe(res)

    def test_reads_gene_gtf_gz(self) -> None:
        res = MagicMock(spec=R3Resource)
        res._cached_path.return_value = _GENE_GTF_GZ
        df = _read_gtf_dataframe(res)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2000
        expected_cols = {"seqname", "source", "feature", "start", "end",
                         "score", "strand", "frame", "attributes"}
        assert expected_cols.issubset(set(df.columns))
        assert set(df["feature"].unique()) == {"gene"}

    def test_reads_plain_gtf(self, tmp_path: Path) -> None:
        gtf_content = (
            "chr1\tref\tgene\t1\t100\t.\t+\t.\t"
            'gene_id "ENSG001"; biotype "pc";\n'
        )
        plain_path = tmp_path / "test.gtf"
        plain_path.write_text(gtf_content)

        res = MagicMock(spec=R3Resource)
        res._cached_path.return_value = plain_path
        df = _read_gtf_dataframe(res)
        assert len(df) == 1
        assert df["feature"].iloc[0] == "gene"


class TestRangesFromGtf:
    def _make_gtf_df(self, rows: list[tuple[str, str, str]]) -> pd.DataFrame:
        """Create a minimal GTF DataFrame; rows = (feature, attrs, seqname)."""
        data = []
        for i, (feat, attrs, seqname) in enumerate(rows):
            data.append({
                "seqname": seqname,
                "source": "ref",
                "feature": feat,
                "start": i * 100 + 1,
                "end": i * 100 + 100,
                "score": ".",
                "strand": "+",
                "frame": ".",
                "attributes": attrs,
            })
        return pd.DataFrame(data)

    def test_empty_dataframe_for_feature_kind(self) -> None:
        gtf = self._make_gtf_df([("exon", 'exon_id "E1"', "chr1")])
        result = _ranges_from_gtf(gtf, feature_kind="gene")
        assert result.empty

    def test_gene_with_gene_id(self) -> None:
        gtf = self._make_gtf_df([("gene", 'gene_id "ENSG001"', "chr1")])
        result = _ranges_from_gtf(gtf, feature_kind="gene")
        assert len(result) == 1
        assert str(result["feature_id"].iloc[0]) == "ENSG001"

    def test_exon_with_recount_exon_id(self) -> None:
        gtf = self._make_gtf_df([
            ("exon", 'recount_exon_id "RCE001"', "chr1")
        ])
        result = _ranges_from_gtf(gtf, feature_kind="exon")
        assert len(result) == 1
        assert str(result["feature_id"].iloc[0]) == "RCE001"

    def test_exon_with_exon_id_fallback(self) -> None:
        gtf = self._make_gtf_df([("exon", 'exon_id "EX001"', "chr1")])
        result = _ranges_from_gtf(gtf, feature_kind="exon")
        assert len(result) == 1
        assert str(result["feature_id"].iloc[0]) == "EX001"

    def test_exon_coord_fallback_when_no_id(self) -> None:
        gtf = self._make_gtf_df([("exon", 'biotype "pc"', "chr1")])
        result = _ranges_from_gtf(gtf, feature_kind="exon")
        assert len(result) == 1
        feat_id = str(result["feature_id"].iloc[0])
        assert "chr1" in feat_id

    def test_duplicate_gene_ids_same_coords_dropped(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        gtf = pd.DataFrame([
            {"seqname": "chr1", "source": "ref", "feature": "gene",
             "start": 1, "end": 100, "score": ".", "strand": "+",
             "frame": ".", "attributes": 'gene_id "DUP"'},
            {"seqname": "chr1", "source": "ref", "feature": "gene",
             "start": 1, "end": 100, "score": ".", "strand": "+",
             "frame": ".", "attributes": 'gene_id "DUP"'},
        ])
        with caplog.at_level(logging.INFO):
            result = _ranges_from_gtf(gtf, feature_kind="gene")
        assert len(result) == 1

    def test_duplicate_gene_ids_conflicting_coords_raises(self) -> None:
        gtf = pd.DataFrame([
            {"seqname": "chr1", "source": "ref", "feature": "gene",
             "start": 1, "end": 100, "score": ".", "strand": "+",
             "frame": ".", "attributes": 'gene_id "DUP"'},
            {"seqname": "chr2", "source": "ref", "feature": "gene",
             "start": 999, "end": 1999, "score": ".", "strand": "-",
             "frame": ".", "attributes": 'gene_id "DUP"'},
        ])
        with pytest.raises(ValueError, match="conflicting"):
            _ranges_from_gtf(gtf, feature_kind="gene")

    def test_level_column_coerced_to_int64(self) -> None:
        gtf = pd.DataFrame([{
            "seqname": "chr1",
            "source": "ref",
            "feature": "gene",
            "start": 1,
            "end": 100,
            "score": ".",
            "strand": "+",
            "frame": ".",
            "attributes": 'gene_id "G1"; level "2"',
        }])
        result = _ranges_from_gtf(gtf, feature_kind="gene")
        assert "level" in result.columns
        assert result["level"].dtype == pd.Int64Dtype()

    def test_reads_gene_gtf_gz(self) -> None:
        res = MagicMock(spec=R3Resource)
        res._cached_path.return_value = _GENE_GTF_GZ
        gtf = _read_gtf_dataframe(res)
        result = _ranges_from_gtf(gtf, feature_kind="gene")
        assert not result.empty
        assert len(result) == 2000
        assert "seqnames" in result.columns
        assert "feature_id" in result.columns
        assert "ENSG00000278704.1" in result["feature_id"].values


class TestPeekGtfFeatureCounts:
    def test_reads_gtf_gz(self) -> None:
        res = MagicMock(spec=R3Resource)
        res._cached_path.return_value = _GENE_GTF_GZ
        counts = _peek_gtf_feature_counts(res)
        assert "gene" in counts
        assert counts["gene"] == 2000

    def test_max_lines_limits_scan(self, tmp_path: Path) -> None:
        lines = ["chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id \"G%d\"\n" % i
                 for i in range(20)]
        gz_path = tmp_path / "test.gtf.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.writelines(lines)

        res = MagicMock(spec=R3Resource)
        res._cached_path.return_value = gz_path
        counts = _peek_gtf_feature_counts(res, max_lines=5)
        assert counts["gene"] == 5

    def test_comment_lines_skipped(self, tmp_path: Path) -> None:
        content = "# comment\nchr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id \"G1\"\n"
        gz_path = tmp_path / "test.gtf.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write(content)

        res = MagicMock(spec=R3Resource)
        res._cached_path.return_value = gz_path
        counts = _peek_gtf_feature_counts(res)
        assert counts.get("gene", 0) == 1
        assert counts.get("#", 0) == 0

    def test_downloads_when_cached_path_raises(self, tmp_path: Path) -> None:
        gz_path = tmp_path / "test.gtf.gz"
        content = "chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id \"G1\"\n"
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write(content)

        call_count = {"n": 0}
        def cached_path_side_effect() -> Path:
            call_count["n"] += 1
            if call_count["n"] == 1:
                raise RuntimeError("not cached")
            return gz_path

        res = MagicMock(spec=R3Resource)
        res._cached_path.side_effect = cached_path_side_effect
        counts = _peek_gtf_feature_counts(res)
        res.download.assert_called_once()
        assert counts["gene"] == 1


class TestSelectGtfResourceForUnit:
    def test_no_annotation_resources_returns_none(self) -> None:
        bundle = R3ResourceBundle()
        result = _select_gtf_resource_for_unit(
            bundle, genomic_unit="gene", annotation_extension=None
        )
        assert result is None

    def test_returns_candidate_containing_feature(self, tmp_path: Path) -> None:
        gz_path = tmp_path / "genes.gtf.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write('chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "G1"\n')

        res = MagicMock(spec=R3Resource)
        res.url = "http://example.com/gene_sums/human.gene_sums.G026.gtf.gz"
        desc = MagicMock()
        desc.resource_type = "annotations"
        desc.url_path.return_value = "human/annotations/gene_sums/human.gene_sums.G026.gtf.gz"
        desc.genomic_unit = "gene"
        desc.annotation_extension = "G026"
        res.description = desc
        res._cached_path.return_value = gz_path

        bundle = R3ResourceBundle(resources=[res])
        result = _select_gtf_resource_for_unit(
            bundle, genomic_unit="gene", annotation_extension=None
        )
        assert result is res

    def test_filters_by_annotation_extension(self, tmp_path: Path) -> None:
        gz_path = tmp_path / "genes.gtf.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write('chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "G1"\n')

        res_wrong = MagicMock(spec=R3Resource)
        desc_wrong = MagicMock()
        desc_wrong.resource_type = "annotations"
        desc_wrong.url_path.return_value = "ann/human.gene_sums.G999.gtf.gz"
        desc_wrong.annotation_extension = "G999"
        desc_wrong.genomic_unit = "gene"
        res_wrong.description = desc_wrong
        res_wrong.url = "http://x/G999.gtf.gz"
        res_wrong._cached_path.return_value = gz_path

        res_right = MagicMock(spec=R3Resource)
        desc_right = MagicMock()
        desc_right.resource_type = "annotations"
        desc_right.url_path.return_value = "ann/human.gene_sums.G026.gtf.gz"
        desc_right.annotation_extension = "G026"
        desc_right.genomic_unit = "gene"
        res_right.description = desc_right
        res_right.url = "http://x/G026.gtf.gz"
        res_right._cached_path.return_value = gz_path

        bundle = R3ResourceBundle(resources=[res_wrong, res_right])
        result = _select_gtf_resource_for_unit(
            bundle, genomic_unit="gene", annotation_extension="G026"
        )
        assert result is res_right

    def test_returns_best_when_no_feature_found(self, tmp_path: Path) -> None:
        gz_path = tmp_path / "exon_only.gtf.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write("chr1\tref\texon\t1\t100\t.\t+\t.\texon_id \"E1\"\n")

        res = MagicMock(spec=R3Resource)
        desc = MagicMock()
        desc.resource_type = "annotations"
        desc.url_path.return_value = "ann/exon_only.gtf.gz"
        desc.annotation_extension = "G026"
        desc.genomic_unit = "gene"
        res.description = desc
        res.url = "http://x/exon_only.gtf.gz"
        res._cached_path.return_value = gz_path

        bundle = R3ResourceBundle(resources=[res])
        result = _select_gtf_resource_for_unit(
            bundle, genomic_unit="gene", annotation_extension=None
        )
        assert result is res


class TestToGenomicRanges:
    def test_uses_from_pandas_when_available(self) -> None:
        mock_gr_cls = MagicMock()
        mock_instance = MagicMock()
        mock_gr_cls.from_pandas.return_value = mock_instance

        ranges_df = pd.DataFrame({
            "seqnames": ["chr1"],
            "starts": [1],
            "ends": [100],
            "strand": ["+"],
        })
        with patch.object(_utils_module, "get_genomicranges_class", return_value=mock_gr_cls):
            result = _to_genomic_ranges(ranges_df)
        assert result is mock_instance
        mock_gr_cls.from_pandas.assert_called_once_with(ranges_df)


class TestConstructSummarizedExperiment:
    def test_raises_if_not_2d(self) -> None:
        counts = pd.DataFrame({"a": [1, 2]}).iloc[:, 0]
        with pytest.raises((ValueError, AttributeError)):
            _construct_summarized_experiment(
                counts_df=counts,  # type: ignore[arg-type]
                row_df=pd.DataFrame(),
                col_df=pd.DataFrame(),
                assay_name="counts",
            )

    def test_raises_empty_assay(self) -> None:
        counts = pd.DataFrame()
        with pytest.raises(ValueError, match="Empty assay"):
            _construct_summarized_experiment(
                counts_df=counts,
                row_df=pd.DataFrame(),
                col_df=pd.DataFrame(),
                assay_name="counts",
            )

    def test_raises_non_numeric_values(self) -> None:
        counts = pd.DataFrame({"s1": ["abc", "def"]}, index=["g1", "g2"])
        row_df = pd.DataFrame({"a": [1, 2]})
        col_df = pd.DataFrame({"b": ["x"]})
        with pytest.raises(ValueError, match="non-numeric"):
            _construct_summarized_experiment(
                counts_df=counts,
                row_df=row_df,
                col_df=col_df,
                assay_name="counts",
            )

    def test_valid_construction(self) -> None:
        counts = _gene_df()
        row_df = pd.DataFrame({"gene_name": ["g1", "g2"]})
        col_df = pd.DataFrame({"sample": ["s1", "s2"]})
        se = _construct_summarized_experiment(
            counts_df=counts,
            row_df=row_df,
            col_df=col_df,
            assay_name="raw_counts",
        )
        assert se is not None

    def test_metadata_passed_through(self) -> None:
        counts = _gene_df()
        row_df = pd.DataFrame({"a": [1, 2]})
        col_df = pd.DataFrame({"b": ["x", "y"]})
        se = _construct_summarized_experiment(
            counts_df=counts,
            row_df=row_df,
            col_df=col_df,
            assay_name="raw",
            metadata={"project": "SRP001"},
        )
        assert se is not None


class TestConstructRangedSummarizedExperiment:
    def _ranges(self, n: int) -> pd.DataFrame:
        return pd.DataFrame({
            "seqnames": ["chr1"] * n,
            "starts": list(range(1, n + 1)),
            "ends": list(range(100, 100 + n)),
            "strand": ["+"] * n,
        })

    def test_raises_if_not_2d(self) -> None:
        counts = pd.Series([1, 2])
        with pytest.raises((ValueError, AttributeError)):
            _construct_ranged_summarized_experiment(
                counts_df=counts,  # type: ignore[arg-type]
                row_df=pd.DataFrame(),
                col_df=pd.DataFrame(),
                ranges_df=pd.DataFrame(),
                assay_name="raw",
            )

    def test_raises_empty_assay(self) -> None:
        with pytest.raises(ValueError, match="Empty assay"):
            _construct_ranged_summarized_experiment(
                counts_df=pd.DataFrame(),
                row_df=pd.DataFrame(),
                col_df=pd.DataFrame(),
                ranges_df=pd.DataFrame(),
                assay_name="raw",
            )

    def test_raises_missing_range_columns(self) -> None:
        counts = _gene_df()
        ranges = pd.DataFrame({"seqnames": ["chr1", "chr1"]})
        with pytest.raises(ValueError, match="missing required columns"):
            _construct_ranged_summarized_experiment(
                counts_df=counts,
                row_df=pd.DataFrame({"a": [1, 2]}),
                col_df=pd.DataFrame({"b": ["x", "y"]}),
                ranges_df=ranges,
                assay_name="raw",
            )

    def test_raises_ranges_length_mismatch(self) -> None:
        counts = _gene_df()
        ranges = self._ranges(1)
        with pytest.raises(ValueError, match="ranges_df length"):
            _construct_ranged_summarized_experiment(
                counts_df=counts,
                row_df=pd.DataFrame({"a": [1, 2]}),
                col_df=pd.DataFrame({"b": ["x", "y"]}),
                ranges_df=ranges,
                assay_name="raw",
            )

    def test_raises_missing_coordinate_values(self) -> None:
        counts = _gene_df()
        ranges = pd.DataFrame({
            "seqnames": [pd.NA, "chr1"],
            "starts": [1, 2],
            "ends": [10, 20],
            "strand": ["+", "+"],
        })
        with pytest.raises(ValueError, match="missing values"):
            _construct_ranged_summarized_experiment(
                counts_df=counts,
                row_df=pd.DataFrame({"a": [1, 2]}),
                col_df=pd.DataFrame({"b": ["x", "y"]}),
                ranges_df=ranges,
                assay_name="raw",
            )

    def test_raises_non_numeric_counts(self) -> None:
        counts = pd.DataFrame({"s1": ["abc", "def"]}, index=["g1", "g2"])
        ranges = self._ranges(2)
        with pytest.raises(ValueError, match="non-numeric"):
            _construct_ranged_summarized_experiment(
                counts_df=counts,
                row_df=pd.DataFrame({"a": [1, 2]}),
                col_df=pd.DataFrame({"b": ["x"]}),
                ranges_df=ranges,
                assay_name="raw",
            )

    def test_valid_construction(self) -> None:
        counts = _gene_df()
        ranges = self._ranges(2)
        rse = _construct_ranged_summarized_experiment(
            counts_df=counts,
            row_df=pd.DataFrame({"a": [1, 2]}),
            col_df=pd.DataFrame({"b": ["x", "y"]}),
            ranges_df=ranges,
            assay_name="raw_counts",
        )
        assert rse is not None

    def test_metadata_passed_through(self) -> None:
        counts = _gene_df()
        ranges = self._ranges(2)
        rse = _construct_ranged_summarized_experiment(
            counts_df=counts,
            row_df=pd.DataFrame({"a": [1, 2]}),
            col_df=pd.DataFrame({"b": ["x", "y"]}),
            ranges_df=ranges,
            assay_name="raw",
            metadata={"project": "SRP001"},
        )
        assert rse is not None


class TestCountCompatKeys:
    def test_gene_or_exon_type(self) -> None:
        res = _mock_resource(
            "count_files_gene_or_exon",
            genomic_unit="gene",
        )
        family, feature_key = _count_compat_keys(res)
        assert family == "gene_or_exon"
        assert "gene" in feature_key

    def test_exon_type(self) -> None:
        res = _mock_resource("count_files_gene_or_exon", genomic_unit="exon")
        family, feature_key = _count_compat_keys(res)
        assert family == "gene_or_exon"
        assert "exon" in feature_key

    def test_junctions_type(self) -> None:
        res = _mock_resource(
            "count_files_junctions",
            junction_type="ALL",
            junction_extension="MM",
        )
        family, feature_key = _count_compat_keys(res)
        assert family == "junctions"
        assert "MM" in feature_key

    def test_unknown_type_raises(self) -> None:
        res = _mock_resource("bigwig_files")
        with pytest.raises(ValueError, match="not a recognized count-file type"):
            _count_compat_keys(res)


class TestMakeUniqueNames:
    def test_no_duplicates(self) -> None:
        assert _make_unique_names(["a", "b", "c"]) == ["a", "b", "c"]

    def test_duplicates_suffixed(self) -> None:
        result = _make_unique_names(["a", "b", "a"])
        assert result == ["a", "b", "a__dup2"]

    def test_triple_duplicates(self) -> None:
        result = _make_unique_names(["x", "x", "x"])
        assert result == ["x", "x__dup2", "x__dup3"]

    def test_empty_list(self) -> None:
        assert _make_unique_names([]) == []

    def test_custom_suffix(self) -> None:
        result = _make_unique_names(["a", "a"], suffix="_copy")
        assert result == ["a", "a_copy2"]


class TestDedupeRangesOnFeatureId:
    def _ranges(self, fids: list[str]) -> pd.DataFrame:
        return pd.DataFrame({
            "feature_id": fids,
            "seqnames": ["chr1"] * len(fids),
            "starts": list(range(1, len(fids) + 1)),
            "ends": list(range(100, 100 + len(fids))),
            "strand": ["+"] * len(fids),
        })

    def test_missing_feature_id_raises(self) -> None:
        df = pd.DataFrame({"seqnames": ["chr1"]})
        with pytest.raises(ValueError, match="missing required column"):
            _dedupe_ranges_on_feature_id(df)

    def test_no_duplicates_returns_unchanged(self) -> None:
        df = self._ranges(["A", "B"])
        result = _dedupe_ranges_on_feature_id(df)
        assert len(result) == 2

    def test_missing_coord_columns_raises(self) -> None:
        df = pd.DataFrame({"feature_id": ["A", "A"], "seqnames": ["chr1", "chr2"]})
        with pytest.raises(ValueError, match="missing one or more"):
            _dedupe_ranges_on_feature_id(df)

    def test_consistent_duplicates_deduped_with_warning(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        df = pd.DataFrame({
            "feature_id": ["A", "A"],
            "seqnames": ["chr1", "chr1"],
            "starts": [1, 1],
            "ends": [100, 100],
            "strand": ["+", "+"],
        })
        with caplog.at_level(logging.WARNING):
            result = _dedupe_ranges_on_feature_id(df)
        assert len(result) == 1
        assert "duplicate feature_id" in caplog.text

    def test_inconsistent_duplicates_raises(self) -> None:
        df = pd.DataFrame({
            "feature_id": ["A", "A"],
            "seqnames": ["chr1", "chr2"],
            "starts": [1, 999],
            "ends": [100, 1100],
            "strand": ["+", "-"],
        })
        with pytest.raises(ValueError, match="different coordinates"):
            _dedupe_ranges_on_feature_id(df)


class TestR3ResourceBundleBasic:
    def test_default_empty_bundle(self) -> None:
        b = R3ResourceBundle()
        assert b.resources == []
        assert b.organism is None
        assert b.data_source is None
        assert b.project is None

    def test_add_resource(self) -> None:
        b = R3ResourceBundle()
        res = _mock_resource()
        b.add(res)
        assert len(b.resources) == 1

    def test_extend_resources(self) -> None:
        b = R3ResourceBundle()
        res1 = _mock_resource()
        res2 = _mock_resource()
        b.extend([res1, res2])
        assert len(b.resources) == 2

    def test_extend_with_generator(self) -> None:
        b = R3ResourceBundle()
        b.extend((_mock_resource() for _ in range(3)))
        assert len(b.resources) == 3


class TestR3ResourceBundleDiscover:
    def _fake_search(self, *args: Any, **kwargs: Any) -> list[MagicMock]:
        desc = MagicMock()
        desc.resource_type = "count_files_gene_or_exon"
        desc.url_path.return_value = f"{kwargs.get('project', 'P')}/counts.gz"
        res = MagicMock(spec=R3Resource)
        res.url = f"http://example.com/{kwargs.get('project', 'P')}/counts.gz"
        res.description = desc
        return [res]

    def test_empty_organism_raises(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            R3ResourceBundle.discover(organism=[], data_source="sra", project="P1")

    def test_empty_data_source_raises(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            R3ResourceBundle.discover(organism="human", data_source=[], project="P1")

    def test_empty_project_raises(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            R3ResourceBundle.discover(organism="human", data_source="sra", project=[])

    def test_single_project_sets_identity(self) -> None:
        with patch.object(search_module, "search_project_all", side_effect=self._fake_search):
            b = R3ResourceBundle.discover(
                organism="human",
                data_source="sra",
                project="SRP001",
            )
        assert b.organism == "human"
        assert b.data_source == "sra"
        assert b.project == "SRP001"
        assert len(b.resources) >= 1

    def test_multi_project_clears_identity(self) -> None:
        with patch.object(search_module, "search_project_all", side_effect=self._fake_search):
            b = R3ResourceBundle.discover(
                organism="human",
                data_source="sra",
                project=["SRP001", "SRP002"],
            )
        assert b.organism is None
        assert b.data_source is None
        assert b.project is None

    def test_deduplicate_removes_duplicate_urls(self) -> None:
        def same_url_search(*args: Any, **kwargs: Any) -> list[MagicMock]:
            desc = MagicMock()
            desc.resource_type = "count_files_gene_or_exon"
            desc.url_path.return_value = "same/path.gz"
            res = MagicMock(spec=R3Resource)
            res.url = "http://example.com/same/path.gz"
            res.description = desc
            return [res]

        with patch.object(search_module, "search_project_all", side_effect=same_url_search):
            b = R3ResourceBundle.discover(
                organism="human",
                data_source="sra",
                project=["SRP001", "SRP002"],
                deduplicate=True,
            )
        assert len(b.resources) == 1

    def test_no_resources_discovered_empty_bundle(self) -> None:
        with patch.object(search_module, "search_project_all", return_value=[]):
            b = R3ResourceBundle.discover(
                organism="human",
                data_source="sra",
                project="SRP999",
            )
        assert len(b.resources) == 0


class TestR3ResourceBundleLoad:
    def test_loads_all_resources(self) -> None:
        res1 = _mock_resource()
        res2 = _mock_resource()
        b = R3ResourceBundle(resources=[res1, res2])
        b.load()
        res1.load.assert_called_once()
        res2.load.assert_called_once()

    def test_returns_self_for_chaining(self) -> None:
        b = R3ResourceBundle()
        result = b.load()
        assert result is b

    def test_strict_raises_on_first_error(self) -> None:
        res = _mock_resource()
        res.load.side_effect = RuntimeError("fail")
        b = R3ResourceBundle(resources=[res])
        with pytest.raises(RuntimeError):
            b.load(strict=True)

    def test_non_strict_skips_errors(self) -> None:
        res1 = _mock_resource()
        res1.load.side_effect = RuntimeError("fail")
        res2 = _mock_resource()
        b = R3ResourceBundle(resources=[res1, res2])
        b.load(strict=False)
        res2.load.assert_called_once()

    def test_force_flag_passed_through(self) -> None:
        res = _mock_resource()
        b = R3ResourceBundle(resources=[res])
        b.load(force=True)
        res.load.assert_called_once_with(force=True)


class TestR3ResourceBundleIterLoaded:
    def test_yields_only_loaded_resources(self) -> None:
        df = _gene_df()
        loaded = _mock_resource(loaded_data=df)
        unloaded = _mock_resource()
        b = R3ResourceBundle(resources=[loaded, unloaded])
        results = list(b.iter_loaded())
        assert len(results) == 1
        assert results[0][0] is loaded

    def test_autoload_triggers_load(self) -> None:
        df = _gene_df()
        res = _mock_resource()
        res.is_loaded.return_value = False
        res.load.side_effect = lambda: setattr(res, "_auto_loaded", True)

        b = R3ResourceBundle(resources=[res])
        res.get_loaded.return_value = df
        res.is_loaded.side_effect = [False, True]
        results = list(b.iter_loaded(autoload=True))
        res.load.assert_called_once()

    def test_autoload_skips_on_exception(self) -> None:
        res = _mock_resource()
        res.is_loaded.return_value = False
        res.load.side_effect = RuntimeError("fail")
        b = R3ResourceBundle(resources=[res])
        results = list(b.iter_loaded(autoload=True))
        assert results == []

    def test_resource_type_filter(self) -> None:
        df = _gene_df()
        gene_res = _mock_resource("count_files_gene_or_exon", loaded_data=df)
        meta_res = _mock_resource("metadata_files", loaded_data=pd.DataFrame())
        b = R3ResourceBundle(resources=[gene_res, meta_res])
        results = list(b.iter_loaded(resource_type="count_files_gene_or_exon"))
        assert all(r.description.resource_type == "count_files_gene_or_exon" for r, _ in results)

    def test_skips_none_loaded_data(self) -> None:
        res = MagicMock(spec=R3Resource)
        desc = MagicMock()
        desc.resource_type = "count_files_gene_or_exon"
        res.description = desc
        res.is_loaded.return_value = True
        res.get_loaded.return_value = None
        b = R3ResourceBundle(resources=[res])
        results = list(b.iter_loaded())
        assert results == []

    def test_get_loaded_returns_list(self) -> None:
        df = _gene_df()
        res = _mock_resource(loaded_data=df)
        b = R3ResourceBundle(resources=[res])
        items = b.get_loaded()
        assert len(items) == 1
        assert items[0] is df


class TestR3ResourceBundleIterBigwig:
    def test_yields_bigwig_resources(self) -> None:
        from recount3._bigwig import BigWigFile
        bw = MagicMock(spec=BigWigFile)
        res = _mock_resource("bigwig_files", loaded_data=bw)
        b = R3ResourceBundle(resources=[res])
        results = list(b.iter_bigwig(autoload=False))
        assert len(results) == 1
        assert results[0][1] is bw

    def test_skips_non_bigwig_objects(self) -> None:
        # Loaded data is a DataFrame, not a BigWigFile
        res = _mock_resource("bigwig_files", loaded_data=_gene_df())
        b = R3ResourceBundle(resources=[res])
        results = list(b.iter_bigwig(autoload=False))
        assert results == []


class TestR3ResourceBundleFilter:
    def test_filter_by_resource_type_string(self) -> None:
        gene_res = _mock_resource("count_files_gene_or_exon")
        meta_res = _mock_resource("metadata_files")
        b = R3ResourceBundle(resources=[gene_res, meta_res])
        filtered = b.filter(resource_type="count_files_gene_or_exon")
        assert len(filtered.resources) == 1

    def test_filter_by_resource_type_tuple(self) -> None:
        gene_res = _mock_resource("count_files_gene_or_exon")
        jxn_res = _mock_resource("count_files_junctions")
        meta_res = _mock_resource("metadata_files")
        b = R3ResourceBundle(resources=[gene_res, jxn_res, meta_res])
        filtered = b.filter(
            resource_type=("count_files_gene_or_exon", "count_files_junctions")
        )
        assert len(filtered.resources) == 2

    def test_filter_invert(self) -> None:
        gene_res = _mock_resource("count_files_gene_or_exon")
        meta_res = _mock_resource("metadata_files")
        b = R3ResourceBundle(resources=[gene_res, meta_res])
        filtered = b.filter(resource_type="metadata_files", invert=True)
        assert len(filtered.resources) == 1
        assert filtered.resources[0].description.resource_type == "count_files_gene_or_exon"

    def test_filter_by_predicate(self) -> None:
        res1 = _mock_resource("count_files_gene_or_exon")
        res1.url = "http://example.com/gene.gz"
        res2 = _mock_resource("count_files_gene_or_exon")
        res2.url = "http://example.com/other.gz"
        b = R3ResourceBundle(resources=[res1, res2])
        filtered = b.filter(predicate=lambda r: "gene" in (r.url or ""))
        assert len(filtered.resources) == 1

    def test_predicate_exception_treated_as_false(self) -> None:
        res = _mock_resource()
        b = R3ResourceBundle(resources=[res])
        filtered = b.filter(predicate=lambda r: 1 / 0)  # type: ignore[arg-type]
        assert len(filtered.resources) == 0

    def test_preserves_bundle_identity(self) -> None:
        b = R3ResourceBundle(
            resources=[_mock_resource()],
            organism="human",
            data_source="sra",
            project="SRP001",
        )
        filtered = b.filter(resource_type="count_files_gene_or_exon")
        assert filtered.organism == "human"
        assert filtered.data_source == "sra"
        assert filtered.project == "SRP001"

    def test_no_criteria_returns_all(self) -> None:
        b = R3ResourceBundle(resources=[_mock_resource(), _mock_resource()])
        assert len(b.filter().resources) == 2


class TestR3ResourceBundleConvenienceFilters:
    def _bundle_with_types(self, *types: str) -> R3ResourceBundle:
        return R3ResourceBundle(resources=[_mock_resource(t) for t in types])

    def test_only_counts(self) -> None:
        b = self._bundle_with_types(
            "count_files_gene_or_exon", "count_files_junctions", "metadata_files"
        )
        result = b.only_counts()
        assert len(result.resources) == 2

    def test_only_metadata(self) -> None:
        b = self._bundle_with_types(
            "count_files_gene_or_exon", "metadata_files"
        )
        result = b.only_metadata()
        assert len(result.resources) == 1

    def test_exclude_metadata(self) -> None:
        b = self._bundle_with_types("count_files_gene_or_exon", "metadata_files")
        result = b.exclude_metadata()
        assert len(result.resources) == 1
        assert result.resources[0].description.resource_type != "metadata_files"

    def test_where_delegates_to_filter(self) -> None:
        res = _mock_resource()
        res.url = "http://example.com/target.gz"
        b = R3ResourceBundle(resources=[res])
        result = b.where(lambda r: "target" in (r.url or ""))
        assert len(result.resources) == 1

    def test_counts_alias(self) -> None:
        b = self._bundle_with_types("count_files_gene_or_exon")
        assert len(b.counts().resources) == 1

    def test_metadata_alias(self) -> None:
        b = self._bundle_with_types("metadata_files")
        assert len(b.metadata().resources) == 1

    def test_bigwigs(self) -> None:
        b = self._bundle_with_types("bigwig_files", "count_files_gene_or_exon")
        assert len(b.bigwigs().resources) == 1


class TestResolveProjectIdentity:
    def test_uses_stored_identity(self) -> None:
        b = R3ResourceBundle(organism="human", data_source="sra", project="SRP001")
        org, src, proj = b._resolve_project_identity(None, None, None)
        assert org == "human"
        assert src == "sra"
        assert proj == "SRP001"

    def test_explicit_overrides_stored(self) -> None:
        b = R3ResourceBundle(organism="human", data_source="sra", project="SRP001")
        org, src, proj = b._resolve_project_identity("human", "sra", "SRP001")
        assert proj == "SRP001"

    def test_missing_organism_raises(self) -> None:
        b = R3ResourceBundle(data_source="sra", project="SRP001")
        with pytest.raises(ValueError, match="incomplete"):
            b._resolve_project_identity(None, None, None)

    def test_missing_data_source_raises(self) -> None:
        b = R3ResourceBundle(organism="human", project="SRP001")
        with pytest.raises(ValueError, match="incomplete"):
            b._resolve_project_identity(None, None, None)

    def test_missing_project_raises(self) -> None:
        b = R3ResourceBundle(organism="human", data_source="sra")
        with pytest.raises(ValueError, match="incomplete"):
            b._resolve_project_identity(None, None, None)

    def test_conflicting_organism_raises(self) -> None:
        b = R3ResourceBundle(organism="human", data_source="sra", project="SRP001")
        with pytest.raises(ValueError, match="does not match"):
            b._resolve_project_identity("mouse", None, None)

    def test_conflicting_data_source_raises(self) -> None:
        b = R3ResourceBundle(organism="human", data_source="sra", project="SRP001")
        with pytest.raises(ValueError, match="does not match"):
            b._resolve_project_identity(None, "gtex", None)

    def test_conflicting_project_raises(self) -> None:
        b = R3ResourceBundle(organism="human", data_source="sra", project="SRP001")
        with pytest.raises(ValueError, match="does not match"):
            b._resolve_project_identity(None, None, "SRP999")


class TestR3ResourceBundleSamples:
    def test_calls_samples_for_project(self) -> None:
        b = R3ResourceBundle(organism="human", data_source="sra", project="SRP001")
        with patch.object(
            search_module,
            "samples_for_project",
            return_value=["SRR001", "SRR002"],
        ) as mock_fn:
            result = b.samples()
        mock_fn.assert_called_once_with(
            organism="human", data_source="sra", project="SRP001"
        )
        assert result == ["SRR001", "SRR002"]

    def test_raises_when_identity_missing(self) -> None:
        b = R3ResourceBundle()
        with pytest.raises(ValueError):
            b.samples()


class TestR3ResourceBundleStackCountMatrices:
    def test_raises_no_count_resources(self) -> None:
        b = R3ResourceBundle(resources=[_mock_resource("metadata_files")])
        with pytest.raises(ValueError, match="No count-file resources"):
            b.stack_count_matrices()

    def test_raises_mixed_families_compat_family(self) -> None:
        gene_res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=_gene_df(),
            genomic_unit="gene",
        )
        jxn_res = _mock_resource(
            "count_files_junctions",
            loaded_data=_gene_df(),
            junction_type="ALL",
            junction_extension="MM",
        )
        b = R3ResourceBundle(resources=[gene_res, jxn_res])
        with pytest.raises(CompatibilityError, match="Incompatible count families"):
            b.stack_count_matrices(compat="family")

    def test_raises_mixed_features_compat_feature(self) -> None:
        gene_res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=_gene_df(),
            genomic_unit="gene",
        )
        exon_res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=_gene_df(),
            genomic_unit="exon",
        )
        b = R3ResourceBundle(resources=[gene_res, exon_res])
        with pytest.raises(CompatibilityError, match="Feature-level incompatibility"):
            b.stack_count_matrices(compat="feature")

    def test_raises_unknown_compat(self) -> None:
        res = _mock_resource("count_files_gene_or_exon", loaded_data=_gene_df(), genomic_unit="gene")
        b = R3ResourceBundle(resources=[res])
        with pytest.raises(ValueError, match="Unknown compat mode"):
            b.stack_count_matrices(compat="invalid")  # type: ignore[arg-type]

    def test_raises_no_loaded_frames(self) -> None:
        res = _mock_resource("count_files_gene_or_exon", genomic_unit="gene")
        b = R3ResourceBundle(resources=[res])
        with pytest.raises(ValueError, match="No loaded count matrices"):
            b.stack_count_matrices(autoload=False)

    def test_raises_non_dataframe_loaded(self) -> None:
        res = MagicMock(spec=R3Resource)
        desc = MagicMock()
        desc.resource_type = "count_files_gene_or_exon"
        desc.url_path.return_value = "x.gz"
        desc.genomic_unit = "gene"
        res.description = desc
        res.url = "http://example.com/x.gz"
        res.is_loaded.return_value = True
        res.get_loaded.return_value = "not-a-dataframe"
        b = R3ResourceBundle(resources=[res])
        with pytest.raises(TypeError, match="not a pandas.DataFrame"):
            b.stack_count_matrices(autoload=False)

    def test_stacks_single_resource(self) -> None:
        df = _gene_df()
        res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=df,
            genomic_unit="gene",
        )
        b = R3ResourceBundle(resources=[res])
        result = b.stack_count_matrices(autoload=False)
        assert isinstance(result, pd.DataFrame)
        assert result.shape == df.shape

    def test_stacks_two_gene_resources(self) -> None:
        df1 = _gene_df(samples=["SRR001", "SRR002"])
        df2 = _gene_df(samples=["SRR003", "SRR004"])
        res1 = _mock_resource("count_files_gene_or_exon", loaded_data=df1, genomic_unit="gene")
        res2 = _mock_resource("count_files_gene_or_exon", loaded_data=df2, genomic_unit="gene")
        b = R3ResourceBundle(resources=[res1, res2])
        result = b.stack_count_matrices(autoload=False, axis=1)
        assert result.shape[1] == 4


class TestR3ResourceBundleStackCountsFor:
    def test_stacks_gene_resources(self) -> None:
        df = _gene_df()
        res = _mock_resource("count_files_gene_or_exon", loaded_data=df, genomic_unit="gene")
        b = R3ResourceBundle(resources=[res])
        result = b._stack_counts_for(genomic_unit="gene", autoload=False)
        assert isinstance(result, pd.DataFrame)

    def test_stacks_exon_resources(self) -> None:
        df = _gene_df()
        res = _mock_resource("count_files_gene_or_exon", loaded_data=df, genomic_unit="exon")
        b = R3ResourceBundle(resources=[res])
        result = b._stack_counts_for(genomic_unit="exon", autoload=False)
        assert isinstance(result, pd.DataFrame)

    def test_stacks_junction_resources(self) -> None:
        df = _gene_df()
        res = _mock_resource(
            "count_files_junctions",
            loaded_data=df,
            junction_type="ALL",
            junction_extension="MM",
        )
        b = R3ResourceBundle(resources=[res])
        result = b._stack_counts_for(genomic_unit="junction", autoload=False)
        assert isinstance(result, pd.DataFrame)

    def test_raises_with_load_errors_gene(self) -> None:
        res = _mock_resource("count_files_gene_or_exon", genomic_unit="gene")
        res.load.side_effect = RuntimeError("load failed")
        b = R3ResourceBundle(resources=[res])
        with pytest.raises(ValueError, match="Failed to load any gene/exon"):
            b._stack_counts_for(genomic_unit="gene", autoload=True)

    def test_raises_with_load_errors_junction(self) -> None:
        res = _mock_resource(
            "count_files_junctions",
            junction_type="ALL",
            junction_extension="MM",
        )
        res.load.side_effect = RuntimeError("load failed")
        b = R3ResourceBundle(resources=[res])
        with pytest.raises(ValueError, match="Failed to load any junction"):
            b._stack_counts_for(genomic_unit="junction", autoload=True)


class TestNormalizeSampleMetadata:
    def test_no_metadata_returns_external_id_only(self) -> None:
        b = R3ResourceBundle()
        result = b._normalize_sample_metadata(sample_ids=["SRR001", "SRR002"])
        assert "external_id" in result.columns
        assert list(result["external_id"]) == ["SRR001", "SRR002"]

    def test_with_metadata_frames(self) -> None:
        meta_df = pd.DataFrame({
            "external_id": pd.array(["SRR001", "SRR002"], dtype="string"),
            "rail_id": pd.array(["1", "2"], dtype="string"),
            "study": pd.array(["SRP001", "SRP001"], dtype="string"),
            "score": [99.0, 88.0],
        })
        res = _mock_resource("metadata_files", loaded_data=meta_df, table_name="recount_qc")
        b = R3ResourceBundle(resources=[res])
        result = b._normalize_sample_metadata(sample_ids=["SRR001", "SRR002"])
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2

    def test_fills_missing_external_id(self) -> None:
        meta_df = pd.DataFrame({
            "external_id": pd.array([pd.NA, "SRR002"], dtype="string"),
            "rail_id": pd.array(["1", "2"], dtype="string"),
            "study": pd.array(["SRP001", "SRP001"], dtype="string"),
        })
        res = _mock_resource("metadata_files", loaded_data=meta_df, table_name="recount_qc")
        b = R3ResourceBundle(resources=[res])
        result = b._normalize_sample_metadata(sample_ids=["1", "SRR002"])
        assert result["external_id"].notna().all()

    def test_provenance_stored_in_attrs(self) -> None:
        meta_df = pd.DataFrame({
            "external_id": pd.array(["SRR001"], dtype="string"),
            "rail_id": pd.array(["1"], dtype="string"),
            "study": pd.array(["SRP001"], dtype="string"),
            "score": [1.0],
        })
        res = _mock_resource("metadata_files", loaded_data=meta_df, table_name="qc_table")
        b = R3ResourceBundle(resources=[res])
        result = b._normalize_sample_metadata(sample_ids=["SRR001"])
        assert "recount3_metadata_provenance" in result.attrs


class TestToSummarizedExperiment:
    def _bundle_with_gene_counts(
        self,
        df: pd.DataFrame | None = None,
        annotation_ext: str = "G026",
    ) -> R3ResourceBundle:
        counts_df = df if df is not None else _gene_df()
        res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=counts_df,
            genomic_unit="gene",
            annotation_extension=annotation_ext,
        )
        return R3ResourceBundle(resources=[res])

    def test_basic_construction(self) -> None:
        b = self._bundle_with_gene_counts()
        se = b.to_summarized_experiment(genomic_unit="gene", autoload=False)
        assert se is not None

    def test_with_annotation_extension_filter(self) -> None:
        df = _gene_df()
        res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=df,
            genomic_unit="gene",
            annotation_extension="G026",
        )
        b = R3ResourceBundle(resources=[res])
        se = b.to_summarized_experiment(
            genomic_unit="gene",
            annotation_extension="G026",
            autoload=False,
        )
        assert se is not None

    def test_duplicate_feature_ids_made_unique(self, caplog: pytest.LogCaptureFixture) -> None:
        df = pd.DataFrame(
            [[1.0, 2.0], [3.0, 4.0]],
            index=["ENSG0001", "ENSG0001"],
            columns=["SRR001", "SRR002"],
        )
        res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=df,
            genomic_unit="gene",
        )
        b = R3ResourceBundle(resources=[res])
        with caplog.at_level(logging.WARNING):
            se = b.to_summarized_experiment(genomic_unit="gene", autoload=False)
        assert se is not None
        assert "duplicate feature IDs" in caplog.text


class TestToRangedSummarizedExperiment:
    @pytest.fixture()
    def _synthetic_gtf_gz(self, tmp_path: Path) -> Path:
        """Write a minimal 3-gene GTF.gz and return its path."""
        gz_path = tmp_path / "genes.gtf.gz"
        content = (
            'chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG001.1"\n'
            'chr1\tref\tgene\t200\t300\t.\t-\t.\tgene_id "ENSG002.1"\n'
            'chr2\tref\tgene\t500\t600\t.\t+\t.\tgene_id "ENSG003.1"\n'
        )
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write(content)
        return gz_path

    def _bundle_with_gene_gtf(self, gtf_path: Path) -> R3ResourceBundle:
        """Bundle with one gene count resource + GTF annotation."""
        feature_ids = ["ENSG001.1", "ENSG002.1", "ENSG003.1"]
        counts_df = pd.DataFrame(
            np.ones((3, 2), dtype=float),
            index=feature_ids,
            columns=["SRR001", "SRR002"],
        )
        res_count = MagicMock(spec=R3Resource)
        desc_count = MagicMock()
        desc_count.resource_type = "count_files_gene_or_exon"
        desc_count.url_path.return_value = "human/gene.gz"
        desc_count.genomic_unit = "gene"
        desc_count.annotation_extension = "G026"
        res_count.description = desc_count
        res_count.url = "http://example.com/gene.gz"
        res_count.is_loaded.return_value = True
        res_count.get_loaded.return_value = counts_df

        res_ann = MagicMock(spec=R3Resource)
        desc_ann = MagicMock()
        desc_ann.resource_type = "annotations"
        desc_ann.url_path.return_value = "human/annotations/gene_sums/human.gene_sums.G026.gtf.gz"
        desc_ann.annotation_extension = "G026"
        desc_ann.genomic_unit = "gene"
        res_ann.description = desc_ann
        res_ann.url = "http://example.com/gene_sums/human.gene_sums.G026.gtf.gz"
        res_ann._cached_path.return_value = gtf_path

        return R3ResourceBundle(resources=[res_count, res_ann])

    def test_gene_with_gtf_annotation(self, _synthetic_gtf_gz: Path) -> None:
        b = self._bundle_with_gene_gtf(_synthetic_gtf_gz)
        rse = b.to_ranged_summarized_experiment(genomic_unit="gene", autoload=False)
        assert rse is not None

    def test_fallback_to_se_when_no_ranges(self) -> None:
        df = _gene_df()
        res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=df,
            genomic_unit="gene",
        )
        b = R3ResourceBundle(resources=[res])
        result = b.to_ranged_summarized_experiment(
            genomic_unit="gene",
            autoload=False,
            allow_fallback_to_se=True,
        )
        assert result is not None

    def test_raises_when_no_ranges_and_no_fallback(self) -> None:
        df = _gene_df()
        res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=df,
            genomic_unit="gene",
        )
        b = R3ResourceBundle(resources=[res])
        with pytest.raises(ValueError, match="Could not derive genomic ranges"):
            b.to_ranged_summarized_experiment(
                genomic_unit="gene",
                autoload=False,
                allow_fallback_to_se=False,
            )

    def test_junction_with_rr_coordinates(self, tmp_path: Path) -> None:
        rr_content = "seqnames\tstarts\tends\tstrand\tjunction_id\nchr1\t1\t100\t+\tJX001\nchr1\t200\t300\t+\tJX002\n"
        rr_gz = tmp_path / "jxn.RR.gz"
        with gzip.open(rr_gz, "wt") as f:
            f.write(rr_content)

        mm_df = pd.DataFrame(
            np.ones((2, 2), dtype=float),
            index=["0", "1"],
            columns=["SRR001", "SRR002"],
        )

        res_mm = _mock_resource(
            "count_files_junctions",
            loaded_data=mm_df,
            junction_type="ALL",
            junction_extension="MM",
        )
        res_rr = MagicMock(spec=R3Resource)
        desc_rr = MagicMock()
        desc_rr.resource_type = "count_files_junctions"
        desc_rr.url_path.return_value = "jxn/RR.gz"
        desc_rr.junction_extension = "RR"
        desc_rr.junction_type = "ALL"
        res_rr.description = desc_rr
        res_rr.url = "http://example.com/jxn/RR.gz"
        res_rr.is_loaded.return_value = True
        rr_df = pd.read_csv(io.StringIO(rr_content), sep="\t")
        res_rr.load.return_value = rr_df
        res_rr.get_loaded.return_value = rr_df
        res_rr._cached_path.return_value = rr_gz

        b = R3ResourceBundle(resources=[res_mm, res_rr])
        rse = b.to_ranged_summarized_experiment(
            genomic_unit="junction",
            prefer_rr_junction_coordinates=True,
            autoload=False,
        )
        assert rse is not None

    def test_junction_rr_missing_required_columns_falls_back(self) -> None:
        bad_rr_df = pd.DataFrame({"col_a": [1, 2]})

        mm_df = _gene_df()
        res_mm = _mock_resource(
            "count_files_junctions",
            loaded_data=mm_df,
            junction_type="ALL",
            junction_extension="MM",
        )
        res_rr = MagicMock(spec=R3Resource)
        desc_rr = MagicMock()
        desc_rr.resource_type = "count_files_junctions"
        desc_rr.url_path.return_value = "jxn/RR.gz"
        desc_rr.junction_extension = "RR"
        desc_rr.junction_type = "ALL"
        res_rr.description = desc_rr
        res_rr.url = "http://example.com/jxn/RR.gz"
        res_rr.is_loaded.return_value = True
        res_rr.load.return_value = bad_rr_df
        res_rr.get_loaded.return_value = bad_rr_df

        b = R3ResourceBundle(resources=[res_mm, res_rr])
        with pytest.raises(ValueError, match="Could not derive"):
            b.to_ranged_summarized_experiment(
                genomic_unit="junction",
                prefer_rr_junction_coordinates=True,
                autoload=False,
                allow_fallback_to_se=False,
            )

    def test_junction_no_rr_resource_raises_when_no_fallback(self) -> None:
        mm_df = _gene_df()
        res_mm = _mock_resource(
            "count_files_junctions",
            loaded_data=mm_df,
            junction_type="ALL",
            junction_extension="MM",
        )
        b = R3ResourceBundle(resources=[res_mm])
        with pytest.raises(ValueError, match="Could not derive"):
            b.to_ranged_summarized_experiment(
                genomic_unit="junction",
                prefer_rr_junction_coordinates=True,
                autoload=False,
                allow_fallback_to_se=False,
            )

    def test_no_rr_fallback_to_se_allowed(self) -> None:
        mm_df = _gene_df()
        res_mm = _mock_resource(
            "count_files_junctions",
            loaded_data=mm_df,
            junction_type="ALL",
            junction_extension="MM",
        )
        b = R3ResourceBundle(resources=[res_mm])
        result = b.to_ranged_summarized_experiment(
            genomic_unit="junction",
            prefer_rr_junction_coordinates=True,
            autoload=False,
            allow_fallback_to_se=True,
        )
        assert result is not None

    def test_raises_with_last_ranges_error_chained(self) -> None:
        df = _gene_df(features=["NONEXISTENT_GENE_XYZ"])
        res_count = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=df,
            genomic_unit="gene",
        )
        res_ann = MagicMock(spec=R3Resource)
        desc_ann = MagicMock()
        desc_ann.resource_type = "annotations"
        desc_ann.url_path.return_value = "ann/gene.gtf.gz"
        desc_ann.annotation_extension = None
        desc_ann.genomic_unit = "gene"
        res_ann.description = desc_ann
        res_ann.url = "http://example.com/gene.gtf.gz"
        res_ann._cached_path.return_value = _GENE_GTF_GZ

        b = R3ResourceBundle(resources=[res_count, res_ann])
        with pytest.raises(ValueError):
            b.to_ranged_summarized_experiment(
                genomic_unit="gene",
                autoload=False,
                allow_fallback_to_se=False,
            )

    def test_rr_row_count_mismatch_falls_back(self) -> None:
        rr_content = (
            "seqnames\tstarts\tends\tstrand\tjunction_id\n"
            "chr1\t1\t100\t+\tJX001\n"
            "chr1\t200\t300\t+\tJX002\n"
            "chr1\t400\t500\t+\tJX003\n"
        )
        rr_df = pd.read_csv(io.StringIO(rr_content), sep="\t")
        mm_df = _gene_df(features=["0", "1"], samples=["SRR001"])

        res_mm = _mock_resource(
            "count_files_junctions",
            loaded_data=mm_df,
            junction_type="ALL",
            junction_extension="MM",
        )
        res_rr = MagicMock(spec=R3Resource)
        desc_rr = MagicMock()
        desc_rr.resource_type = "count_files_junctions"
        desc_rr.url_path.return_value = "jxn/RR.gz"
        desc_rr.junction_extension = "RR"
        desc_rr.junction_type = "ALL"
        res_rr.description = desc_rr
        res_rr.url = "http://example.com/jxn/RR.gz"
        res_rr.load.return_value = rr_df
        res_rr.get_loaded.return_value = rr_df
        res_rr.is_loaded.return_value = True

        b = R3ResourceBundle(resources=[res_mm, res_rr])
        with pytest.raises(ValueError, match="Could not derive"):
            b.to_ranged_summarized_experiment(
                genomic_unit="junction",
                prefer_rr_junction_coordinates=True,
                autoload=False,
                allow_fallback_to_se=False,
            )

    def test_rr_duplicate_row_names_made_unique(self, tmp_path: Path) -> None:
        rr_content = (
            "seqnames\tstarts\tends\tstrand\tjunction_id\n"
            "chr1\t1\t100\t+\tJX_DUP\n"
            "chr1\t200\t300\t+\tJX_DUP\n"
        )
        rr_df = pd.read_csv(io.StringIO(rr_content), sep="\t")
        mm_df = pd.DataFrame(
            np.ones((2, 1), dtype=float),
            index=["0", "1"],
            columns=["SRR001"],
        )
        res_mm = _mock_resource(
            "count_files_junctions",
            loaded_data=mm_df,
            junction_type="ALL",
            junction_extension="MM",
        )
        res_rr = MagicMock(spec=R3Resource)
        desc_rr = MagicMock()
        desc_rr.resource_type = "count_files_junctions"
        desc_rr.url_path.return_value = "jxn/RR.gz"
        desc_rr.junction_extension = "RR"
        desc_rr.junction_type = "ALL"
        res_rr.description = desc_rr
        res_rr.url = "http://example.com/jxn/RR.gz"
        res_rr.load.return_value = rr_df
        res_rr.get_loaded.return_value = rr_df
        res_rr.is_loaded.return_value = True

        b = R3ResourceBundle(resources=[res_mm, res_rr])
        rse = b.to_ranged_summarized_experiment(
            genomic_unit="junction",
            prefer_rr_junction_coordinates=True,
            autoload=False,
        )
        assert rse is not None

    def test_rr_without_strand_column_defaults_to_star(self, tmp_path: Path) -> None:
        rr_content = "seqnames\tstarts\tends\nchr1\t1\t100\nchr1\t200\t300\n"
        rr_df = pd.read_csv(io.StringIO(rr_content), sep="\t")
        mm_df = pd.DataFrame(
            np.ones((2, 1), dtype=float),
            index=["0", "1"],
            columns=["SRR001"],
        )
        res_mm = _mock_resource(
            "count_files_junctions",
            loaded_data=mm_df,
            junction_type="ALL",
            junction_extension="MM",
        )
        res_rr = MagicMock(spec=R3Resource)
        desc_rr = MagicMock()
        desc_rr.resource_type = "count_files_junctions"
        desc_rr.url_path.return_value = "jxn/RR.gz"
        desc_rr.junction_extension = "RR"
        res_rr.description = desc_rr
        res_rr.url = "http://example.com/jxn/RR.gz"
        res_rr.load.return_value = rr_df
        res_rr.get_loaded.return_value = rr_df
        res_rr.is_loaded.return_value = True

        b = R3ResourceBundle(resources=[res_mm, res_rr])
        rse = b.to_ranged_summarized_experiment(
            genomic_unit="junction",
            prefer_rr_junction_coordinates=True,
            autoload=False,
        )
        assert rse is not None

    def test_duplicate_feature_ids_in_counts(self, caplog: pytest.LogCaptureFixture) -> None:
        df = pd.DataFrame(
            [[1.0], [2.0]],
            index=["GENE_DUP", "GENE_DUP"],
            columns=["SRR001"],
        )
        res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=df,
            genomic_unit="gene",
        )
        b = R3ResourceBundle(resources=[res])
        with caplog.at_level(logging.WARNING):
            result = b.to_ranged_summarized_experiment(
                genomic_unit="gene",
                autoload=False,
                allow_fallback_to_se=True,
            )
        assert result is not None


class TestR3ResourceBundleDownload:
    def test_calls_download_on_each_resource(self) -> None:
        res1 = _mock_resource()
        res2 = _mock_resource()
        b = R3ResourceBundle(resources=[res1, res2])
        b.download(dest="/some/dir")
        res1.download.assert_called_once_with(
            path="/some/dir", cache_mode="enable", overwrite=False
        )
        res2.download.assert_called_once_with(
            path="/some/dir", cache_mode="enable", overwrite=False
        )

    def test_overwrite_flag_forwarded(self) -> None:
        res = _mock_resource()
        b = R3ResourceBundle(resources=[res])
        b.download(dest=".", overwrite=True)
        res.download.assert_called_once_with(
            path=".", cache_mode="enable", overwrite=True
        )

    def test_cache_mode_forwarded(self) -> None:
        res = _mock_resource()
        b = R3ResourceBundle(resources=[res])
        b.download(dest=".", cache="update")
        res.download.assert_called_once_with(
            path=".", cache_mode="update", overwrite=False
        )

    def test_empty_bundle_does_nothing(self) -> None:
        b = R3ResourceBundle()
        b.download()


class TestModuleConstants:
    def test_gtf_attr_pair_re_matches(self) -> None:
        m = bmod._GTF_ATTR_PAIR_RE.search('gene_id "ENSG001"')
        assert m is not None
        assert m.group(1) == "gene_id"
        assert m.group(2) == "ENSG001"

    def test_exon_id_re_matches(self) -> None:
        m = bmod._EXON_ID_ATTR_RE.search('exon_id "E001"')
        assert m is not None
        assert m.group(1) == "E001"

    def test_recount_exon_id_re_matches(self) -> None:
        m = bmod._RECOUNT_EXON_ID_ATTR_RE.search('recount_exon_id "RCE001"')
        assert m is not None
        assert m.group(1) == "RCE001"

    def test_metadata_namespace_separator(self) -> None:
        assert bmod._METADATA_NAMESPACE_SEPARATOR == "__"

    def test_metadata_merge_keys(self) -> None:
        assert "rail_id" in bmod._METADATA_MERGE_KEYS
        assert "external_id" in bmod._METADATA_MERGE_KEYS
        assert "study" in bmod._METADATA_MERGE_KEYS


class TestParseGtfAttributesExtractedEmpty:
    def test_non_matching_attrs_returns_empty_df(self) -> None:
        attrs = pd.Series(["."])
        result = _parse_gtf_attributes(attrs)
        assert isinstance(result, pd.DataFrame)

    def test_single_word_no_value_returns_empty(self) -> None:
        attrs = pd.Series(["justkey"])
        result = _parse_gtf_attributes(attrs)
        assert isinstance(result, pd.DataFrame)


class TestRangesFromGtfAdditional:
    def test_exon_with_dot_attrs_uses_coord_fallback(self) -> None:
        gtf = pd.DataFrame([{
            "seqname": "chr1",
            "source": "ref",
            "feature": "exon",
            "start": 1,
            "end": 100,
            "score": ".",
            "strand": "+",
            "frame": ".",
            "attributes": ".",
        }])
        result = _ranges_from_gtf(gtf, feature_kind="exon")
        assert len(result) == 1
        assert "chr1" in str(result["feature_id"].iloc[0])

    def test_gene_with_empty_attrs_uses_coord_fallback(self) -> None:
        gtf = pd.DataFrame([{
            "seqname": "chrX",
            "source": "ref",
            "feature": "gene",
            "start": 500,
            "end": 600,
            "score": ".",
            "strand": "-",
            "frame": ".",
            "attributes": ".",
        }])
        result = _ranges_from_gtf(gtf, feature_kind="gene")
        assert len(result) == 1

    def test_extra_attrs_joined_to_output(self) -> None:
        gtf = pd.DataFrame([{
            "seqname": "chr1",
            "source": "ref",
            "feature": "gene",
            "start": 1,
            "end": 100,
            "score": ".",
            "strand": "+",
            "frame": ".",
            "attributes": 'gene_id "G1"; gene_name "MYC"',
        }])
        result = _ranges_from_gtf(gtf, feature_kind="gene")
        assert "gene_name" in result.columns


class TestPeekGtfShortLines:
    def test_line_with_two_fields_skipped(self, tmp_path: Path) -> None:
        gz_path = tmp_path / "short.gtf.gz"
        content = (
            "chr1\tref\n"
            'chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "G1"\n'
        )
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write(content)
        res = MagicMock(spec=R3Resource)
        res._cached_path.return_value = gz_path
        counts = _peek_gtf_feature_counts(res)
        assert counts.get("gene", 0) == 1
        assert counts.get("ref", 0) == 0


class TestSelectGtfResourceForUnitAdditional:
    def test_exon_unit_scoring(self, tmp_path: Path) -> None:
        gz_path = tmp_path / "exons.gtf.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write('chr1\tref\texon\t1\t100\t.\t+\t.\texon_id "E1"\n')

        res = MagicMock(spec=R3Resource)
        desc = MagicMock()
        desc.resource_type = "annotations"
        desc.url_path.return_value = "ann/exon_sums.gtf.gz"
        desc.annotation_extension = "G026"
        desc.genomic_unit = "exon"
        res.description = desc
        res.url = "http://x/exon_sums.gtf.gz"
        res._cached_path.return_value = gz_path

        bundle = R3ResourceBundle(resources=[res])
        result = _select_gtf_resource_for_unit(
            bundle, genomic_unit="exon", annotation_extension=None
        )
        assert result is res

    def test_resource_url_without_gtf_still_selected(self, tmp_path: Path) -> None:
        gz_path = tmp_path / "genes.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write('chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "G1"\n')

        res = MagicMock(spec=R3Resource)
        desc = MagicMock()
        desc.resource_type = "annotations"
        desc.url_path.return_value = "ann/genes.gz"
        desc.annotation_extension = "G026"
        desc.genomic_unit = None
        res.description = desc
        res.url = "http://x/genes.gz"
        res._cached_path.return_value = gz_path

        bundle = R3ResourceBundle(resources=[res])
        result = _select_gtf_resource_for_unit(
            bundle, genomic_unit="gene", annotation_extension=None
        )
        assert result is res

    def test_genomic_unit_matches_adds_score(self, tmp_path: Path) -> None:
        """Cover 766->770 True branch (s += 200)."""
        gz_path = tmp_path / "genes.gtf.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write('chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "G1"\n')

        res = MagicMock(spec=R3Resource)
        desc = MagicMock()
        desc.resource_type = "annotations"
        desc.url_path.return_value = "ann/gene_sums.gtf.gz"
        desc.annotation_extension = "G026"
        desc.genomic_unit = "gene"  # matches genomic_unit
        res.description = desc
        res.url = "http://x/gene_sums.gtf.gz"
        res._cached_path.return_value = gz_path

        bundle = R3ResourceBundle(resources=[res])
        result = _select_gtf_resource_for_unit(
            bundle, genomic_unit="gene", annotation_extension=None
        )
        assert result is res


class TestDiscoverDeduplicate:
    """Cover line 1289->1297: deduplicate=False branch."""

    def test_deduplicate_false_skips_dedup_step(self) -> None:
        call_count = {"n": 0}

        def search_side_effect(*args: Any, **kwargs: Any) -> list[MagicMock]:
            call_count["n"] += 1
            desc = MagicMock()
            desc.resource_type = "count_files_gene_or_exon"
            desc.url_path.return_value = f"proj/counts_{call_count['n']}.gz"
            res = MagicMock(spec=R3Resource)
            res.url = f"http://example.com/proj/counts_{call_count['n']}.gz"
            res.description = desc
            return [res]

        with patch.object(search_module, "search_project_all", side_effect=search_side_effect):
            b = R3ResourceBundle.discover(
                organism="human",
                data_source="sra",
                project="SRP001",
                deduplicate=False,
            )
        assert b.organism == "human"


class TestStackCountMatricesValueErrorContinue:
    """Cover lines 1786-1787: _count_compat_keys raises ValueError."""

    def test_compat_keys_value_error_skipped(self) -> None:
        df = _gene_df()
        good_res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=df,
            genomic_unit="gene",
        )
        bad_res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=df,
            genomic_unit="gene",
        )
        b = R3ResourceBundle(resources=[good_res, bad_res])

        original = bmod._count_compat_keys
        call_count = {"n": 0}

        def patched(res: Any) -> tuple[str, str]:
            call_count["n"] += 1
            if call_count["n"] == 2:
                raise ValueError("forced error")
            return original(res)

        with patch.object(bmod, "_count_compat_keys", side_effect=patched):
            result = b.stack_count_matrices(autoload=False)
        assert isinstance(result, pd.DataFrame)


class TestStackCountMatricesSkipsNonCount:
    """Cover line 1824: continue when rtype not in wanted during iter_loaded."""

    def test_loaded_non_count_resource_skipped(self) -> None:
        gene_df = _gene_df()
        gene_res = _mock_resource("count_files_gene_or_exon", loaded_data=gene_df, genomic_unit="gene")
        meta_res = _mock_resource("metadata_files", loaded_data=pd.DataFrame({"a": [1]}))
        b = R3ResourceBundle(resources=[gene_res, meta_res])
        result = b.stack_count_matrices(autoload=False)
        assert isinstance(result, pd.DataFrame)
        assert result.shape == gene_df.shape


class TestStackCountsForRaise:
    """Cover lines 1902, 1930: re-raise when no load errors."""

    def test_gene_reraises_value_error_no_load_errors(self) -> None:
        b = R3ResourceBundle()
        with pytest.raises(ValueError):
            b._stack_counts_for(genomic_unit="gene", autoload=False)

    def test_junction_reraises_value_error_no_load_errors(self) -> None:
        b = R3ResourceBundle()
        with pytest.raises(ValueError):
            b._stack_counts_for(genomic_unit="junction", autoload=False)


class TestNormalizeSampleMetadataNonDataFrame:
    def test_non_dataframe_metadata_object_skipped(self) -> None:
        res = MagicMock(spec=R3Resource)
        desc = MagicMock()
        desc.resource_type = "metadata_files"
        desc.url_path.return_value = "meta/file.gz"
        desc.table_name = "qc"
        res.description = desc
        res.url = "http://example.com/meta.gz"
        res.is_loaded.return_value = True
        res.get_loaded.return_value = "not-a-dataframe"

        b = R3ResourceBundle(resources=[res])
        result = b._normalize_sample_metadata(sample_ids=["SRR001"])
        assert "external_id" in result.columns


class TestNormalizeSampleMetadataExternalIdFillback:
    def test_aligned_without_external_id_gets_filled(self) -> None:
        meta_df = pd.DataFrame({
            "rail_id": pd.array(["1", "2"], dtype="string"),
            "study": pd.array(["SRP001", "SRP001"], dtype="string"),
            "score": [10.0, 20.0],
        })
        res = _mock_resource("metadata_files", loaded_data=meta_df, table_name="qc")
        b = R3ResourceBundle(resources=[res])
        result = b._normalize_sample_metadata(sample_ids=["1", "2"])
        assert "external_id" in result.columns


class TestToRangedSEAutoload:
    def test_autoload_true_downloads_gtf(self, tmp_path: Path) -> None:
        gz_path = tmp_path / "genes.gtf.gz"
        content = (
            'chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG001"; gene_name "MYC"\n'
            'chr2\tref\tgene\t200\t300\t.\t-\t.\tgene_id "ENSG002"; gene_name "TP53"\n'
        )
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write(content)

        counts_df = pd.DataFrame(
            np.ones((2, 1), dtype=float),
            index=["ENSG001", "ENSG002"],
            columns=["SRR001"],
        )
        res_count = MagicMock(spec=R3Resource)
        desc_count = MagicMock()
        desc_count.resource_type = "count_files_gene_or_exon"
        desc_count.url_path.return_value = "human/gene.gz"
        desc_count.genomic_unit = "gene"
        desc_count.annotation_extension = "G026"
        res_count.description = desc_count
        res_count.url = "http://example.com/gene.gz"
        res_count.is_loaded.return_value = True
        res_count.get_loaded.return_value = counts_df

        res_ann = MagicMock(spec=R3Resource)
        desc_ann = MagicMock()
        desc_ann.resource_type = "annotations"
        desc_ann.url_path.return_value = "human/annotations/gene_sums/human.gene_sums.G026.gtf.gz"
        desc_ann.annotation_extension = "G026"
        desc_ann.genomic_unit = "gene"
        res_ann.description = desc_ann
        res_ann.url = "http://example.com/gene_sums/human.gene_sums.G026.gtf.gz"
        res_ann._cached_path.return_value = gz_path

        b = R3ResourceBundle(resources=[res_count, res_ann])
        rse = b.to_ranged_summarized_experiment(genomic_unit="gene", autoload=True)
        res_ann.download.assert_called()
        assert rse is not None

    def test_enrich_cols_joined_to_row_data(self, tmp_path: Path) -> None:
        gz_path = tmp_path / "genes_with_name.gtf.gz"
        content = (
            'chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG001"; gene_name "MYC"\n'
        )
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write(content)

        counts_df = pd.DataFrame(
            [[5.0]], index=["ENSG001"], columns=["SRR001"]
        )
        res_count = MagicMock(spec=R3Resource)
        desc_count = MagicMock()
        desc_count.resource_type = "count_files_gene_or_exon"
        desc_count.url_path.return_value = "human/gene.gz"
        desc_count.genomic_unit = "gene"
        desc_count.annotation_extension = "G026"
        res_count.description = desc_count
        res_count.url = "http://example.com/gene.gz"
        res_count.is_loaded.return_value = True
        res_count.get_loaded.return_value = counts_df

        res_ann = MagicMock(spec=R3Resource)
        desc_ann = MagicMock()
        desc_ann.resource_type = "annotations"
        desc_ann.url_path.return_value = "human/ann/gene.gtf.gz"
        desc_ann.annotation_extension = "G026"
        desc_ann.genomic_unit = "gene"
        res_ann.description = desc_ann
        res_ann.url = "http://example.com/ann/gene.gtf.gz"
        res_ann._cached_path.return_value = gz_path

        b = R3ResourceBundle(resources=[res_count, res_ann])
        rse = b.to_ranged_summarized_experiment(genomic_unit="gene", autoload=False)
        assert rse is not None

    def test_junction_prefer_rr_false(self) -> None:
        mm_df = _gene_df(features=["jxn1", "jxn2"], samples=["SRR001"])
        res_mm = _mock_resource(
            "count_files_junctions",
            loaded_data=mm_df,
            junction_type="ALL",
            junction_extension="MM",
        )
        b = R3ResourceBundle(resources=[res_mm])
        with pytest.raises(ValueError, match="Could not derive"):
            b.to_ranged_summarized_experiment(
                genomic_unit="junction",
                prefer_rr_junction_coordinates=False,
                autoload=False,
                allow_fallback_to_se=False,
            )

    def test_junction_autoload_downloads_rr(self) -> None:
        rr_content = (
            "seqnames\tstarts\tends\tstrand\tjunction_id\n"
            "chr1\t1\t100\t+\tJX001\n"
        )
        rr_df = pd.read_csv(io.StringIO(rr_content), sep="\t")
        mm_df = pd.DataFrame([[1.0]], index=["0"], columns=["SRR001"])

        res_mm = _mock_resource(
            "count_files_junctions",
            loaded_data=mm_df,
            junction_type="ALL",
            junction_extension="MM",
        )
        res_rr = MagicMock(spec=R3Resource)
        desc_rr = MagicMock()
        desc_rr.resource_type = "count_files_junctions"
        desc_rr.url_path.return_value = "jxn/RR.gz"
        desc_rr.junction_extension = "RR"
        desc_rr.junction_type = "ALL"
        res_rr.description = desc_rr
        res_rr.url = "http://example.com/jxn/RR.gz"
        res_rr.load.return_value = rr_df
        res_rr.get_loaded.return_value = rr_df
        res_rr.is_loaded.return_value = True

        b = R3ResourceBundle(resources=[res_mm, res_rr])
        rse = b.to_ranged_summarized_experiment(
            genomic_unit="junction",
            prefer_rr_junction_coordinates=True,
            autoload=True,
        )
        res_rr.download.assert_called()
        assert rse is not None


# Gene IDs from the first rows of the real trimmed GTF
_REAL_GENE_IDS = [
    "ENSG00000278704.1",
    "ENSG00000277400.1",
    "ENSG00000274847.1",
    "ENSG00000277428.1",
    "ENSG00000276256.1",
]


class TestRealGtfFixtureIntegration:
    """Integration tests using the real (trimmed) GTF fixture files.

    These tests confirm that the trimmed fixture GTFs contain valid data and
    that the full GTF-parsing pipeline produces correct results against real
    GENCODE v26 annotation content.
    """

    def test_gene_gtf_has_correct_row_count(self) -> None:
        res = MagicMock(spec=R3Resource)
        res._cached_path.return_value = _GENE_GTF_GZ
        df = _read_gtf_dataframe(res)
        assert len(df) == 2000
        assert set(df["feature"].unique()) == {"gene"}

    def test_exon_gtf_has_correct_row_count(self) -> None:
        res = MagicMock(spec=R3Resource)
        res._cached_path.return_value = _EXON_GTF_GZ
        df = _read_gtf_dataframe(res)
        assert len(df) == 3000
        assert "exon" in df["feature"].unique()

    def test_gene_gtf_attributes_parsed(self) -> None:
        """Real GTF has gene_id, gene_type, gene_name, level attributes."""
        res = MagicMock(spec=R3Resource)
        res._cached_path.return_value = _GENE_GTF_GZ
        df = _read_gtf_dataframe(res)
        ranges = _ranges_from_gtf(df, feature_kind="gene")
        assert not ranges.empty
        assert len(ranges) == 2000
        assert "gene_name" in ranges.columns
        assert "gene_type" in ranges.columns
        assert "level" in ranges.columns
        assert ranges["level"].dtype == pd.Int64Dtype()

    def test_gene_gtf_versioned_ensembl_ids(self) -> None:
        """Real IDs like ENSG00000278704.1 are preserved (not stripped)."""
        res = MagicMock(spec=R3Resource)
        res._cached_path.return_value = _GENE_GTF_GZ
        df = _read_gtf_dataframe(res)
        ranges = _ranges_from_gtf(df, feature_kind="gene")
        for gene_id in _REAL_GENE_IDS:
            assert gene_id in ranges["feature_id"].values, (
                f"{gene_id} not found in GTF ranges"
            )

    def test_exon_gtf_recount_exon_id_parsed(self) -> None:
        """Exon GTF uses recount_exon_id as the feature_id."""
        res = MagicMock(spec=R3Resource)
        res._cached_path.return_value = _EXON_GTF_GZ
        df = _read_gtf_dataframe(res)
        ranges = _ranges_from_gtf(df, feature_kind="exon")
        assert not ranges.empty
        # recount_exon_id format: "seqname|start|end|strand"
        first_id = str(ranges["feature_id"].iloc[0])
        assert "|" in first_id

    def test_end_to_end_ranged_se_with_real_gene_gtf(self) -> None:
        """Count DataFrame + real GTF -> RangedSummarizedExperiment."""
        counts_df = pd.DataFrame(
            np.ones((5, 3), dtype=float),
            index=_REAL_GENE_IDS,
            columns=["SRR001", "SRR002", "SRR003"],
        )

        res_count = MagicMock(spec=R3Resource)
        desc_count = MagicMock()
        desc_count.resource_type = "count_files_gene_or_exon"
        desc_count.url_path.return_value = "human/gene.gz"
        desc_count.genomic_unit = "gene"
        desc_count.annotation_extension = "G026"
        res_count.description = desc_count
        res_count.url = "http://example.com/gene.gz"
        res_count.is_loaded.return_value = True
        res_count.get_loaded.return_value = counts_df

        res_ann = MagicMock(spec=R3Resource)
        desc_ann = MagicMock()
        desc_ann.resource_type = "annotations"
        desc_ann.url_path.return_value = (
            "human/annotations/gene_sums/human.gene_sums.G026.gtf.gz"
        )
        desc_ann.annotation_extension = "G026"
        desc_ann.genomic_unit = "gene"
        res_ann.description = desc_ann
        res_ann.url = "http://example.com/gene_sums/human.gene_sums.G026.gtf.gz"
        res_ann._cached_path.return_value = _GENE_GTF_GZ

        b = R3ResourceBundle(resources=[res_count, res_ann])
        rse = b.to_ranged_summarized_experiment(genomic_unit="gene", autoload=False)

        assert rse is not None
        # Should have 5 aligned features (all 5 IDs exist in the real GTF)
        assert rse.shape[0] == 5
        assert rse.shape[1] == 3
        assert list(rse.colnames) == ["SRR001", "SRR002", "SRR003"]  # type: ignore


class TestReadRrTableExceptionPaths:
    def test_fallback_to_file_when_load_raises(self, tmp_path: Path) -> None:
        content = "col_a\tcol_b\n10\t20\n"
        tsv_path = tmp_path / "rr.tsv"
        tsv_path.write_text(content)

        res = MagicMock(spec=R3Resource)
        res.url = "http://example.com/rr.tsv"
        res.load.side_effect = RuntimeError("load failed")
        res._cached_path.return_value = tsv_path

        result = _read_rr_table(res)
        assert list(result.columns) == ["col_a", "col_b"]
        assert len(result) == 1

    def test_raises_value_error_when_cached_path_raises(self) -> None:
        res = MagicMock(spec=R3Resource)
        res.url = "http://example.com/rr.gz"
        res.load.side_effect = RuntimeError("load failed")
        res._cached_path.side_effect = OSError("no path")

        with pytest.raises(ValueError, match="Cannot resolve cached RR path"):
            _read_rr_table(res)


class TestSelectGtfResourcePeekException:
    def test_peek_exception_is_logged_and_skipped(self, tmp_path: Path) -> None:
        """When _peek_gtf_feature_counts raises, log a warning and continue."""
        gz_path = tmp_path / "genes.gtf.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write('chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "G1"\n')

        res = MagicMock(spec=R3Resource)
        desc = MagicMock()
        desc.resource_type = "annotations"
        desc.url_path.return_value = "ann/genes.gtf.gz"
        desc.annotation_extension = "G026"
        desc.genomic_unit = "gene"
        res.description = desc
        res.url = "http://x/genes.gtf.gz"
        res._cached_path.return_value = gz_path

        bundle = R3ResourceBundle(resources=[res])

        with patch.object(bmod, "_peek_gtf_feature_counts", side_effect=OSError("peek fail")):
            result = _select_gtf_resource_for_unit(
                bundle, genomic_unit="gene", annotation_extension=None
            )
        assert result is res


class TestConstructSELengthGuards:
    def test_raises_when_row_data_length_changed_by_unique_columns(self) -> None:
        counts = _gene_df()
        row_df = pd.DataFrame({"a": [1, 2]})
        col_df = pd.DataFrame({"b": ["x", "y"]})

        original = bmod._ensure_unique_columns
        call_count: dict[str, int] = {"n": 0}

        def patched(df: pd.DataFrame, *, empty_prefix: str = "col") -> pd.DataFrame:
            call_count["n"] += 1
            if call_count["n"] == 1:
                return df.iloc[:1]
            return original(df, empty_prefix=empty_prefix)

        with patch.object(bmod, "_ensure_unique_columns", side_effect=patched):
            with pytest.raises(ValueError, match="row_data length"):
                _construct_summarized_experiment(
                    counts_df=counts,
                    row_df=row_df,
                    col_df=col_df,
                    assay_name="counts",
                )

    def test_raises_when_col_data_length_changed_by_unique_columns(self) -> None:
        counts = _gene_df()
        row_df = pd.DataFrame({"a": [1, 2]})
        col_df = pd.DataFrame({"b": ["x", "y"]})

        original = bmod._ensure_unique_columns
        call_count: dict[str, int] = {"n": 0}

        def patched(df: pd.DataFrame, *, empty_prefix: str = "col") -> pd.DataFrame:
            call_count["n"] += 1
            if call_count["n"] == 2:
                return df.iloc[:1]
            return original(df, empty_prefix=empty_prefix)

        with patch.object(bmod, "_ensure_unique_columns", side_effect=patched):
            with pytest.raises(ValueError, match="column_data length"):
                _construct_summarized_experiment(
                    counts_df=counts,
                    row_df=row_df,
                    col_df=col_df,
                    assay_name="counts",
                )


class TestConstructRSELengthGuards:
    def _ranges(self) -> pd.DataFrame:
        return pd.DataFrame({
            "seqnames": ["chr1", "chr1"],
            "starts": [1, 200],
            "ends": [100, 300],
            "strand": ["+", "+"],
        })

    def test_raises_when_row_data_length_changed_by_unique_columns(self) -> None:
        counts = _gene_df()
        row_df = pd.DataFrame({"a": [1, 2]})
        col_df = pd.DataFrame({"b": ["x", "y"]})

        original = bmod._ensure_unique_columns
        call_count: dict[str, int] = {"n": 0}

        def patched(df: pd.DataFrame, *, empty_prefix: str = "col") -> pd.DataFrame:
            call_count["n"] += 1
            if call_count["n"] == 1:
                return df.iloc[:1]
            return original(df, empty_prefix=empty_prefix)

        with patch.object(bmod, "_ensure_unique_columns", side_effect=patched):
            with pytest.raises(ValueError, match="row_data length"):
                _construct_ranged_summarized_experiment(
                    counts_df=counts,
                    row_df=row_df,
                    col_df=col_df,
                    ranges_df=self._ranges(),
                    assay_name="raw_counts",
                )

    def test_raises_when_col_data_length_changed_by_unique_columns(self) -> None:
        counts = _gene_df()
        row_df = pd.DataFrame({"a": [1, 2]})
        col_df = pd.DataFrame({"b": ["x", "y"]})

        original = bmod._ensure_unique_columns
        call_count: dict[str, int] = {"n": 0}

        def patched(df: pd.DataFrame, *, empty_prefix: str = "col") -> pd.DataFrame:
            call_count["n"] += 1
            if call_count["n"] == 2:
                return df.iloc[:1]
            return original(df, empty_prefix=empty_prefix)

        with patch.object(bmod, "_ensure_unique_columns", side_effect=patched):
            with pytest.raises(ValueError, match="column_data length"):
                _construct_ranged_summarized_experiment(
                    counts_df=counts,
                    row_df=row_df,
                    col_df=col_df,
                    ranges_df=self._ranges(),
                    assay_name="raw_counts",
                )


class TestNormalizeSampleMetadataMissingAlignKey:
    def test_raises_when_collapse_drops_align_key(self) -> None:
        meta_df = pd.DataFrame({"external_id": ["SRR001", "SRR002"]})
        res = _mock_resource("metadata_files", loaded_data=meta_df, table_name="qc")
        b = R3ResourceBundle(resources=[res])

        def drop_key(df: pd.DataFrame, *, key: str) -> pd.DataFrame:
            return df.drop(columns=[key], errors="ignore")

        with patch.object(bmod, "_collapse_rows_by_key", side_effect=drop_key):
            with pytest.raises(ValueError, match="missing alignment key"):
                b._normalize_sample_metadata(sample_ids=["SRR001", "SRR002"])


class TestNormalizeSampleMetadataExternalIdAbsent:
    def test_external_id_set_when_absent_from_aligned(self) -> None:
        merged_df = pd.DataFrame({
            "rail_id": pd.array(["SRR001", "SRR002"], dtype="string"),
            "study":   pd.array(["SRP001", "SRP001"], dtype="string"),
            "extra":   ["a", "b"],
        })

        meta_df = pd.DataFrame({"rail_id": ["SRR001", "SRR002"]})
        res = _mock_resource("metadata_files", loaded_data=meta_df, table_name="qc")
        b = R3ResourceBundle(resources=[res])

        with patch.object(bmod, "_outer_merge_metadata_frames", return_value=merged_df):
            result = b._normalize_sample_metadata(sample_ids=["SRR001", "SRR002"])

        assert "external_id" in result.columns
        assert list(result["external_id"]) == ["SRR001", "SRR002"]


class TestAddBigwigUrls:
    def test_no_external_id_column_returns_na(self) -> None:
        b = R3ResourceBundle()
        col_df = pd.DataFrame({"other": ["SRR001"]})
        result = b._add_bigwig_urls(col_df)
        assert "BigWigURL" in result.columns
        assert pd.isna(result["BigWigURL"].iloc[0])

    def test_no_count_resources_sets_bigwig_na(self) -> None:
        b = R3ResourceBundle()
        col_df = pd.DataFrame({"external_id": pd.array(["SRR001"], dtype="string")})
        result = b._add_bigwig_urls(col_df)
        assert "BigWigURL" in result.columns
        assert pd.isna(result["BigWigURL"].iloc[0])

    def test_non_count_resource_skipped_then_returns_na(self) -> None:
        meta_res = _mock_resource("metadata_files")
        b = R3ResourceBundle(resources=[meta_res])
        col_df = pd.DataFrame({"external_id": pd.array(["SRR001"], dtype="string")})
        result = b._add_bigwig_urls(col_df)
        assert pd.isna(result["BigWigURL"].iloc[0])

    def test_count_resource_missing_organism_continues_then_na(self) -> None:
        res = _mock_resource(
            "count_files_gene_or_exon",
            organism=None,
            data_source="sra",
            project="SRP001",
        )
        b = R3ResourceBundle(resources=[res])
        col_df = pd.DataFrame({"external_id": pd.array(["SRR001"], dtype="string")})
        result = b._add_bigwig_urls(col_df)
        assert pd.isna(result["BigWigURL"].iloc[0])

    def test_na_external_id_appends_none_url(self) -> None:
        res = _mock_resource(
            "count_files_gene_or_exon",
            organism="human",
            data_source="sra",
            project="SRP001",
        )
        b = R3ResourceBundle(resources=[res])
        col_df = pd.DataFrame({"external_id": pd.array([pd.NA], dtype="string")})
        result = b._add_bigwig_urls(col_df)
        assert "BigWigURL" in result.columns
        assert result["BigWigURL"].iloc[0] is None

    def test_file_source_col_is_found_and_applied(self) -> None:
        res = _mock_resource(
            "count_files_gene_or_exon",
            organism="human",
            data_source="sra",
            project="SRP001",
        )
        b = R3ResourceBundle(resources=[res])
        col_df = pd.DataFrame({
            "external_id": pd.array(["SRR001"], dtype="string"),
            "recount_seq_qc__file_source": ["some/path/sra"],
        })
        mock_bw = MagicMock()
        mock_bw.url_path.return_value = "human/data_sources/sra/base_sums/SRP001/SRR001.bw"
        mock_cfg = MagicMock()
        mock_cfg.base_url = "https://duffel.example.com/recount3"
        with patch("recount3._descriptions.R3BigWig", return_value=mock_bw), \
             patch("recount3.config.default_config", return_value=mock_cfg):
            result = b._add_bigwig_urls(col_df)
        assert "BigWigURL" in result.columns
        assert result["BigWigURL"].iloc[0] is not None

    def test_file_source_col_non_string_value_skips_override(self) -> None:
        res = _mock_resource(
            "count_files_gene_or_exon",
            organism="human",
            data_source="sra",
            project="SRP001",
        )
        b = R3ResourceBundle(resources=[res])
        col_df = pd.DataFrame({
            "external_id": pd.array([pd.NA], dtype="string"),
            "recount_seq_qc__file_source": [42],
        })
        result = b._add_bigwig_urls(col_df)
        assert "BigWigURL" in result.columns

    def test_file_source_col_slash_only_skips_override(self) -> None:
        res = _mock_resource(
            "count_files_gene_or_exon",
            organism="human",
            data_source="sra",
            project="SRP001",
        )
        b = R3ResourceBundle(resources=[res])
        col_df = pd.DataFrame({
            "external_id": pd.array([pd.NA], dtype="string"),
            "recount_seq_qc__file_source": ["/"],
        })
        result = b._add_bigwig_urls(col_df)
        assert "BigWigURL" in result.columns

    def test_successful_bigwig_url_constructed(self) -> None:
        res = _mock_resource(
            "count_files_gene_or_exon",
            organism="human",
            data_source="sra",
            project="SRP001",
        )
        b = R3ResourceBundle(resources=[res])
        col_df = pd.DataFrame({"external_id": pd.array(["SRR001"], dtype="string")})
        mock_bw = MagicMock()
        mock_bw.url_path.return_value = "human/data_sources/sra/base_sums/SRP001/SRR001.bw"
        mock_cfg = MagicMock()
        mock_cfg.base_url = "https://duffel.example.com/recount3"
        with patch("recount3._descriptions.R3BigWig", return_value=mock_bw), \
             patch("recount3.config.default_config", return_value=mock_cfg):
            result = b._add_bigwig_urls(col_df)
        assert "BigWigURL" in result.columns
        url = result["BigWigURL"].iloc[0]
        assert url is not None
        assert url.startswith("https://duffel.example.com/recount3/")
