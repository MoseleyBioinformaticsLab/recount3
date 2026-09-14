# Copyright (c) 2026, Alexander A. Alsalihi, Robert M. Flight,
# Hunter N.B. Moseley. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * All advertising materials mentioning features or use of this software must
#   display the following acknowledgement: This product includes software
#   developed by the copyright holder.
# * Neither the name of the copyright holder nor the names of its contributors
#   may be used to endorse or promote products derived from this software
#   without specific prior written permission.
# * If the source code is used in a published work, then proper citation of the
#   source code must be included with the published work.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# pylint: disable=redefined-outer-name

from __future__ import annotations

import dataclasses
import gzip
import io
import logging
import threading
import urllib.error
import urllib.parse
import zipfile
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
import scipy.sparse

import recount3._utils as _utils_module
import recount3.bundle as bmod
import recount3.resource as resource_module
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
    _classify_ranges_failure,
    _dedupe_ranges_on_feature_id,
    _ensure_unique_columns,
    _load_optional_metadata,
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
from recount3._descriptions import R3ResourceDescription
from recount3.config import Config
from recount3.errors import (
    CompatibilityError,
    DownloadError,
    LoadError,
    MissingRangesError,
    RangesCoverageError,
    RangesError,
)
from recount3.resource import R3Resource

_TESTS_DIR = Path(__file__).parent
_DATA_DIR = _TESTS_DIR / "data"
_MIRROR = _DATA_DIR / "recount3_mirror" / "recount3"

_GENE_GTF_GZ = (
    _MIRROR
    / "human"
    / "annotations"
    / "gene_sums"
    / "human.gene_sums.G026.gtf.gz"
)
_EXON_GTF_GZ = (
    _MIRROR
    / "human"
    / "annotations"
    / "exon_sums"
    / "human.exon_sums.G026.gtf.gz"
)

# Path.as_uri(), not "file://" + str(path): on Windows the latter yields
# file://D:\..., where the drive letter parses as the URL authority rather
# than as part of the path. It is well formed on POSIX only because an
# absolute path already starts with the separator.
_MIRROR_URL = _MIRROR.resolve().as_uri() + "/"


@pytest.fixture()
def local_config(tmp_path: Path) -> Config:
    """Return a Config whose cache_dir is tmp_path and base_url points at the
    local test-data mirror."""
    return Config(
        base_url=_MIRROR_URL,
        timeout=5,
        insecure_ssl=False,
        max_retries=1,
        user_agent="test",
        cache_dir=tmp_path / "cache",
        cache_disabled=False,
    )


def test_mirror_url_is_a_well_formed_file_url() -> None:
    """Guard the fixture that reaches the download code path.

    A "file://" URL carrying a Windows drive letter puts it in the
    authority, where urllib treats it as a remote host and refuses the
    download. This assertion is what fails first if that creeps back.
    """
    parts = urllib.parse.urlsplit(_MIRROR_URL)
    assert parts.scheme == "file"
    assert parts.netloc == ""
    assert parts.path.endswith("/")


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
    res.config = None
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


def _gene_df(
    features: list[str] | None = None, samples: list[str] | None = None
) -> pd.DataFrame:
    """Return a tiny gene-counts DataFrame."""
    rows = features or ["ENSG0001", "ENSG0002"]
    cols = samples or ["SRR001", "SRR002"]
    data = np.arange(len(rows) * len(cols), dtype=float).reshape(
        len(rows), len(cols)
    )
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

    def test_merge_key_text_is_independent_of_parsed_dtype(self) -> None:
        """An int and a float column of the same IDs render identically.

        A single blank cell widens a ``rail_id`` column to ``float64``, so
        without canonical rendering the two tables would carry ``"100"``
        and ``"100.0"`` and could never be joined.
        """
        as_int = _standardize_metadata_frame(pd.DataFrame({"rail_id": [100]}))
        as_float = _standardize_metadata_frame(
            pd.DataFrame({"rail_id": [100.0]})
        )
        assert list(as_int["rail_id"]) == ["100"]
        assert list(as_float["rail_id"]) == list(as_int["rail_id"])

    def test_tables_typed_differently_still_join(self) -> None:
        """The regression this guards: an inner join must not empty out."""
        left = _standardize_metadata_frame(
            pd.DataFrame(
                {
                    "rail_id": [100, 101],
                    "external_id": ["S1", "S2"],
                    "study": ["P", "P"],
                    "a": [1, 2],
                }
            )
        )
        # Same samples, but a blank third row made pandas pick float64.
        right = _standardize_metadata_frame(
            pd.DataFrame(
                {
                    "rail_id": [100.0, 101.0, np.nan],
                    "external_id": ["S1", "S2", "S3"],
                    "study": ["P", "P", "P"],
                    "b": [3, 4, 5],
                }
            )
        )
        merged = pd.merge(
            left, right, on=["rail_id", "external_id", "study"], how="inner"
        )
        assert list(merged["external_id"]) == ["S1", "S2"]

    def test_non_integral_merge_key_kept_verbatim(self) -> None:
        out = _standardize_metadata_frame(pd.DataFrame({"rail_id": [1.5]}))
        assert list(out["rail_id"]) == ["1.5"]

    def test_missing_merge_key_values_stay_missing(self) -> None:
        out = _standardize_metadata_frame(
            pd.DataFrame({"rail_id": [100.0, np.nan]})
        )
        assert list(out["rail_id"])[0] == "100"
        assert pd.isna(list(out["rail_id"])[1])

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
        df = pd.DataFrame(
            {"rail_id": ["1"], "external_id": ["x"], "study": ["s"]}
        )
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
        df = pd.DataFrame(
            {
                "rail_id": ["1"],
                "external_id": ["SRR001"],
                "study": ["SRP001"],
                "extra": ["v"],
            }
        )
        result = _outer_merge_metadata_frames([df])
        assert len(result) == 1
        assert "extra" in result.columns

    def test_two_frames_merged(self) -> None:
        df1 = pd.DataFrame(
            {
                "rail_id": ["1", "2"],
                "external_id": ["A", "B"],
                "study": ["S1", "S1"],
                "col_a": [10, 20],
            }
        )
        df2 = pd.DataFrame(
            {
                "rail_id": ["1", "3"],
                "external_id": ["A", "C"],
                "study": ["S1", "S1"],
                "col_b": [100, 300],
            }
        )
        result = _outer_merge_metadata_frames([df1, df2])
        assert len(result) == 3  # outer merge: rows 1, 2, 3
        assert "col_a" in result.columns
        assert "col_b" in result.columns


class TestChooseAlignmentKey:
    def test_prefers_external_id_when_better_match(self) -> None:
        merged = pd.DataFrame(
            {
                "external_id": pd.array(
                    ["SRR001", "SRR002", "SRR003"], dtype="string"
                ),
                "rail_id": pd.array(["1", "2", "99"], dtype="string"),
            }
        )
        key = _choose_alignment_key(
            sample_ids=["SRR001", "SRR002"],
            merged=merged,
        )
        assert key == "external_id"

    def test_prefers_rail_id_when_better_match(self) -> None:
        merged = pd.DataFrame(
            {
                "external_id": pd.array(["XX", "YY"], dtype="string"),
                "rail_id": pd.array(["1", "2"], dtype="string"),
            }
        )
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
        df = pd.DataFrame(
            {
                "key": ["A", "A", "B"],
                "val": [1, 2, 3],
            }
        )
        result = _collapse_rows_by_key(df, key="key")
        assert len(result) == 2
        assert set(result["key"]) == {"A", "B"}

    def test_first_non_null_used(self) -> None:
        df = pd.DataFrame(
            {
                "key": ["A", "A"],
                "val": [pd.NA, 42],
            }
        )
        result = _collapse_rows_by_key(df, key="key")
        assert result.loc[result["key"] == "A", "val"].iloc[0] == 42

    def test_all_null_returns_na(self) -> None:
        df = pd.DataFrame(
            {
                "key": ["A", "A"],
                "val": [pd.NA, pd.NA],
            }
        )
        result = _collapse_rows_by_key(df, key="key")
        val = result.loc[result["key"] == "A", "val"].iloc[0]
        assert pd.isna(val)


class TestMaybeRelabelCountsColumnsToExternalId:
    def test_no_external_id_col(self) -> None:
        counts = _gene_df(samples=["c1", "c2"])
        col_df = pd.DataFrame({"other": ["x", "y"]}, index=["c1", "c2"])
        out_counts, out_col = _maybe_relabel_counts_columns_to_external_id(
            counts, col_df
        )
        assert list(out_counts.columns) == list(counts.columns)

    def test_length_mismatch(self) -> None:
        counts = _gene_df(samples=["c1", "c2"])
        col_df = pd.DataFrame({"external_id": ["x"]})
        out_counts, out_col = _maybe_relabel_counts_columns_to_external_id(
            counts, col_df
        )
        assert list(out_counts.columns) == list(counts.columns)

    def test_missing_external_id_values(self) -> None:
        counts = _gene_df(samples=["c1", "c2"])
        col_df = pd.DataFrame(
            {"external_id": pd.array(["SRR001", pd.NA], dtype="string")},
            index=["c1", "c2"],
        )
        out_counts, _ = _maybe_relabel_counts_columns_to_external_id(
            counts, col_df
        )
        assert list(out_counts.columns) == ["c1", "c2"]

    def test_non_unique_external_ids(self) -> None:
        counts = _gene_df(samples=["c1", "c2"])
        col_df = pd.DataFrame(
            {"external_id": pd.array(["SRR001", "SRR001"], dtype="string")},
            index=["c1", "c2"],
        )
        out_counts, _ = _maybe_relabel_counts_columns_to_external_id(
            counts, col_df
        )
        assert list(out_counts.columns) == ["c1", "c2"]

    def test_already_matching_no_rename(self) -> None:
        counts = _gene_df(samples=["SRR001", "SRR002"])
        col_df = pd.DataFrame(
            {"external_id": pd.array(["SRR001", "SRR002"], dtype="string")},
            index=["SRR001", "SRR002"],
        )
        out_counts, _ = _maybe_relabel_counts_columns_to_external_id(
            counts, col_df
        )
        assert list(out_counts.columns) == ["SRR001", "SRR002"]

    def test_successful_relabel(self) -> None:
        counts = _gene_df(samples=["1", "2"])
        col_df = pd.DataFrame(
            {"external_id": pd.array(["SRR001", "SRR002"], dtype="string")},
            index=["1", "2"],
        )
        out_counts, out_col = _maybe_relabel_counts_columns_to_external_id(
            counts, col_df
        )
        assert list(out_counts.columns) == ["SRR001", "SRR002"]
        assert list(out_col.index) == ["SRR001", "SRR002"]

    def test_empty_string_external_id_treated_as_missing(self) -> None:
        counts = _gene_df(samples=["c1", "c2"])
        col_df = pd.DataFrame(
            {"external_id": pd.array(["SRR001", ""], dtype="string")},
            index=["c1", "c2"],
        )
        out_counts, _ = _maybe_relabel_counts_columns_to_external_id(
            counts, col_df
        )
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
        attrs = pd.Series(
            [
                'gene_id "ENSG001"; biotype "protein_coding"',
                'gene_id "ENSG002"; biotype "lncRNA"',
            ]
        )
        result = _parse_gtf_attributes(attrs)
        assert len(result) == 2
        assert result["gene_id"].iloc[1] == "ENSG002"

    def test_quoted_value_keeps_embedded_semicolon(self) -> None:
        """A quoted value is taken whole, as rtracklayer does."""
        attrs = pd.Series(['gene_id "A"; gene_name "has;semicolon";'])
        result = _parse_gtf_attributes(attrs)
        assert result["gene_name"].iloc[0] == "has;semicolon"

    def test_quoted_value_with_separators_adds_no_columns(self) -> None:
        """The tail of a quoted value must not become its own attribute.

        Ending the value at a separator inside the quotes would leave
        ``c (d)"`` to be rescanned and parsed as ``c = (d)``, inventing a
        column that the annotation never declared.
        """
        attrs = pd.Series(['gene_id "A"; note "a, b; c (d)"; gene_name "N";'])
        result = _parse_gtf_attributes(attrs)
        assert list(result.columns) == ["gene_id", "note", "gene_name"]
        assert result["note"].iloc[0] == "a, b; c (d)"
        assert result["gene_name"].iloc[0] == "N"

    def test_empty_quoted_value_is_empty_string(self) -> None:
        """``key ""`` is an empty value, not an absent attribute."""
        attrs = pd.Series(['gene_id "A"; gene_name "";'])
        result = _parse_gtf_attributes(attrs)
        assert result["gene_name"].iloc[0] == ""

    def test_unquoted_value_parsed(self) -> None:
        attrs = pd.Series(['gene_id "A"; level 2;'])
        result = _parse_gtf_attributes(attrs)
        assert result["level"].iloc[0] == "2"

    def test_quoted_value_keeps_spaces(self) -> None:
        attrs = pd.Series(['gene_id "A"; gene_name "has space";'])
        result = _parse_gtf_attributes(attrs)
        assert result["gene_name"].iloc[0] == "has space"

    def test_repeated_key_keeps_last_value(self) -> None:
        """Matches rtracklayer, which reports the final occurrence."""
        attrs = pd.Series(['gene_id "A"; tag "x"; tag "y";'])
        result = _parse_gtf_attributes(attrs)
        assert result["tag"].iloc[0] == "y"

    def test_trailing_pair_without_semicolon(self) -> None:
        attrs = pd.Series(['gene_id "A"; gene_name "N"'])
        result = _parse_gtf_attributes(attrs)
        assert result["gene_name"].iloc[0] == "N"


class TestCoerceGtfPhase:
    def test_valid_phases(self) -> None:
        s = pd.Series(["0", "1", "2", ".", pd.NA])
        result = _coerce_gtf_phase_column(s)
        assert result.iloc[0] == 0
        assert result.iloc[1] == 1
        assert result.iloc[2] == 2
        assert pd.isna(result.iloc[3])
        assert pd.isna(result.iloc[4])

    def test_invalid_phase_coerced_to_na_with_warning(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
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
    def test_dot_score_preserves_missing_length(self) -> None:
        score = pd.Series([".", ".", "."])
        starts = pd.Series([1, 100, 200])
        ends = pd.Series([10, 110, 210])
        result = _coerce_gtf_bp_length(score, starts=starts, ends=ends)
        assert pd.isna(result.iloc[0])
        assert pd.isna(result.iloc[1])

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
        assert int(result.iloc[0]) == 999
        assert int(result.iloc[1]) == 999

    def test_score_none_comparable(self) -> None:
        starts = pd.Series([5])
        ends = pd.Series([14])
        score = pd.Series(["."])
        result = _coerce_gtf_bp_length(score, starts=starts, ends=ends)
        assert pd.isna(result.iloc[0])


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
    return pd.DataFrame(
        {
            "feature_id": feature_ids,
            "seqnames": ["chr1"] * len(feature_ids),
            "starts": list(range(1, len(feature_ids) + 1)),
            "ends": list(range(100, 100 + len(feature_ids))),
            "strand": ["+"] * len(feature_ids),
        }
    )


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
        ranges = pd.DataFrame(
            {
                "feature_id": ["ENSG001.1", "ENSG001.2"],
                "seqnames": ["chr1", "chr2"],
                "starts": [1, 200],
                "ends": [100, 300],
                "strand": ["+", "-"],
            }
        )
        with pytest.raises(ValueError, match="conflicting coordinates"):
            _align_ranges_to_features(ranges, feature_ids=["ENSG001"])

    def test_fallback_still_missing_returns_partial(self) -> None:
        ranges = _minimal_ranges(["ENSG001"])
        result = _align_ranges_to_features(
            ranges, feature_ids=["ENSG001", "NOTFOUND"]
        )
        assert pd.isna(result.loc["NOTFOUND", "seqnames"])

    def test_version_strip_with_non_conflicting_duplicates(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        ranges = pd.DataFrame(
            {
                "feature_id": ["ENSG001.1", "ENSG001.2"],
                "seqnames": ["chr1", "chr1"],
                "starts": [100, 100],
                "ends": [200, 200],
                "strand": ["+", "+"],
            }
        )
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
        expected_cols = {
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attributes",
        }
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

    @pytest.mark.parametrize("suffix", [".gtf", ".gtf.gz"])
    def test_normalizes_unstranded_features(
        self, tmp_path: Path, suffix: str
    ) -> None:
        path = tmp_path / f"annotation{suffix}"
        content = "".join(
            f'chr1\tref\tgene\t1\t100\t.\t{strand}\t.\tgene_id "G{i}";\n'
            for i, strand in enumerate([".", "+", "-", "*"])
        )
        opener = gzip.open if suffix.endswith(".gz") else open
        with opener(path, "wt", encoding="utf-8") as handle:
            handle.write(content)
        res = MagicMock(spec=R3Resource)
        res._cached_path.return_value = path

        df = _read_gtf_dataframe(res)

        assert df["strand"].tolist() == ["*", "+", "-", "*"]
        assert df["score"].tolist() == ["."] * 4
        assert df["frame"].tolist() == ["."] * 4
        ranges = _ranges_from_gtf(df, feature_kind="gene")
        bmod._validate_coordinates(ranges)


class TestRangesFromGtf:
    def _make_gtf_df(self, rows: list[tuple[str, str, str]]) -> pd.DataFrame:
        """Create a minimal GTF DataFrame; rows = (feature, attrs, seqname)."""
        data = []
        for i, (feat, attrs, seqname) in enumerate(rows):
            data.append(
                {
                    "seqname": seqname,
                    "source": "ref",
                    "feature": feat,
                    "start": i * 100 + 1,
                    "end": i * 100 + 100,
                    "score": ".",
                    "strand": "+",
                    "frame": ".",
                    "attributes": attrs,
                }
            )
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
        gtf = self._make_gtf_df([("exon", 'recount_exon_id "RCE001"', "chr1")])
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

    def test_duplicate_gene_ids_same_coords_preserved(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        gtf = pd.DataFrame(
            [
                {
                    "seqname": "chr1",
                    "source": "ref",
                    "feature": "gene",
                    "start": 1,
                    "end": 100,
                    "score": ".",
                    "strand": "+",
                    "frame": ".",
                    "attributes": 'gene_id "DUP"',
                },
                {
                    "seqname": "chr1",
                    "source": "ref",
                    "feature": "gene",
                    "start": 1,
                    "end": 100,
                    "score": ".",
                    "strand": "+",
                    "frame": ".",
                    "attributes": 'gene_id "DUP"',
                },
            ]
        )
        with caplog.at_level(logging.INFO):
            result = _ranges_from_gtf(gtf, feature_kind="gene")
        assert len(result) == 2

    def test_duplicate_gene_ids_conflicting_coords_raises(self) -> None:
        gtf = pd.DataFrame(
            [
                {
                    "seqname": "chr1",
                    "source": "ref",
                    "feature": "gene",
                    "start": 1,
                    "end": 100,
                    "score": ".",
                    "strand": "+",
                    "frame": ".",
                    "attributes": 'gene_id "DUP"',
                },
                {
                    "seqname": "chr2",
                    "source": "ref",
                    "feature": "gene",
                    "start": 999,
                    "end": 1999,
                    "score": ".",
                    "strand": "-",
                    "frame": ".",
                    "attributes": 'gene_id "DUP"',
                },
            ]
        )
        with pytest.raises(ValueError, match="conflicting"):
            _ranges_from_gtf(gtf, feature_kind="gene")

    def test_level_column_coerced_to_int64(self) -> None:
        gtf = pd.DataFrame(
            [
                {
                    "seqname": "chr1",
                    "source": "ref",
                    "feature": "gene",
                    "start": 1,
                    "end": 100,
                    "score": ".",
                    "strand": "+",
                    "frame": ".",
                    "attributes": 'gene_id "G1"; level "2"',
                }
            ]
        )
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
        res.ensure_cached.return_value = _GENE_GTF_GZ
        counts = _peek_gtf_feature_counts(res)
        assert "gene" in counts
        assert counts["gene"] == 2000

    def test_autoload_is_passed_to_the_resource(self) -> None:
        """autoload is the bundle's policy; download= is the resource's."""
        res = MagicMock(spec=R3Resource)
        res.ensure_cached.return_value = _GENE_GTF_GZ
        _peek_gtf_feature_counts(res, autoload=False)
        res.ensure_cached.assert_called_once_with(download=False)

    def test_max_lines_limits_scan(self, tmp_path: Path) -> None:
        lines = [
            'chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "G%d"\n' % i
            for i in range(20)
        ]
        gz_path = tmp_path / "test.gtf.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.writelines(lines)

        res = MagicMock(spec=R3Resource)
        res.ensure_cached.return_value = gz_path
        counts = _peek_gtf_feature_counts(res, max_lines=5)
        assert counts["gene"] == 5

    def test_comment_lines_skipped(self, tmp_path: Path) -> None:
        content = '# comment\nchr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "G1"\n'
        gz_path = tmp_path / "test.gtf.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write(content)

        res = MagicMock(spec=R3Resource)
        res.ensure_cached.return_value = gz_path
        counts = _peek_gtf_feature_counts(res)
        assert counts.get("gene", 0) == 1
        assert counts.get("#", 0) == 0


def _warnings(caplog: pytest.LogCaptureFixture) -> list[str]:
    """Return the messages of every captured record at WARNING or above."""
    return [
        rec.getMessage()
        for rec in caplog.records
        if rec.levelno >= logging.WARNING
    ]


class TestAnnotationCachePreparation:
    """Cache preparation for annotations, using real resources and caches.

    These exercise the actual contract of
    :meth:`R3Resource._cached_path`: it computes a path and never consults
    the filesystem, so a cache miss is an absent file, not an exception.
    Mocking ``_cached_path`` to raise cannot reproduce that, which is how
    the first-use ``FileNotFoundError`` warning went unnoticed.
    """

    @staticmethod
    def _gene_annotation(cfg: Config) -> R3Resource:
        """Return the mirror's real human G026 gene annotation resource."""
        return R3Resource(
            description=R3ResourceDescription(
                resource_type="annotations",
                organism="human",
                genomic_unit="gene",
                annotation_extension="G026",
            ),
            config=cfg,
        )

    def test_cache_miss_is_an_absent_file_not_an_exception(
        self, local_config: Config
    ) -> None:
        res = self._gene_annotation(local_config)
        assert not res._cached_path().exists()

    def test_ensure_cached_downloads_on_cold_cache(
        self, local_config: Config
    ) -> None:
        res = self._gene_annotation(local_config)
        path = res.ensure_cached()
        assert path.exists()
        assert path == res._cached_path()

    def test_ensure_cached_reuses_a_warm_cache(
        self, local_config: Config, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        res = self._gene_annotation(local_config)
        res.download(path=None, cache_mode="enable")

        def no_network(*_args: Any, **_kwargs: Any) -> None:
            raise AssertionError("a warm cache must not be re-fetched")

        monkeypatch.setattr(_utils_module, "http_open", no_network)
        assert res.ensure_cached().exists()

    def test_ensure_cached_refuses_to_download_when_told_not_to(
        self, local_config: Config
    ) -> None:
        res = self._gene_annotation(local_config)
        with pytest.raises(FileNotFoundError):
            res.ensure_cached(download=False)
        assert not res._cached_path().exists()

    def test_peek_succeeds_on_a_cold_cache_without_warning(
        self, local_config: Config, caplog: pytest.LogCaptureFixture
    ) -> None:
        """The regression: first use used to warn, then work on the retry."""
        res = self._gene_annotation(local_config)
        with caplog.at_level(logging.DEBUG):
            counts = _peek_gtf_feature_counts(res)
        assert counts["gene"] == 2000
        assert res._cached_path().exists()
        assert _warnings(caplog) == []

    def test_peek_with_autoload_off_raises_and_stays_offline(
        self, local_config: Config
    ) -> None:
        res = self._gene_annotation(local_config)
        with pytest.raises(FileNotFoundError):
            _peek_gtf_feature_counts(res, autoload=False)
        assert not res._cached_path().exists()

    def test_selection_on_a_cold_cache_emits_no_warning(
        self, local_config: Config, caplog: pytest.LogCaptureFixture
    ) -> None:
        res = self._gene_annotation(local_config)
        bundle = R3ResourceBundle(resources=[res])
        with caplog.at_level(logging.DEBUG):
            picked = _select_gtf_resource_for_unit(
                bundle, genomic_unit="gene", annotation_extension="G026"
            )
        assert picked is res
        assert not res._cached_path().exists()  # Selection needs no file I/O.
        assert _warnings(caplog) == []

    def test_selection_on_a_warm_cache_issues_no_request(
        self, local_config: Config, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        res = self._gene_annotation(local_config)
        res.download(path=None, cache_mode="enable")

        def no_network(*_args: Any, **_kwargs: Any) -> None:
            raise AssertionError("a warm cache must not be re-fetched")

        monkeypatch.setattr(_utils_module, "http_open", no_network)
        bundle = R3ResourceBundle(resources=[res])
        picked = _select_gtf_resource_for_unit(
            bundle, genomic_unit="gene", annotation_extension="G026"
        )
        assert picked is res

    def test_selection_with_autoload_off_stays_offline(
        self, local_config: Config, caplog: pytest.LogCaptureFixture
    ) -> None:
        res = self._gene_annotation(local_config)
        bundle = R3ResourceBundle(resources=[res])
        with caplog.at_level(logging.DEBUG):
            picked = _select_gtf_resource_for_unit(
                bundle,
                genomic_unit="gene",
                annotation_extension="G026",
                autoload=False,
            )
        # Ranking alone decides; nothing is fetched and nothing is wrong.
        assert picked is res
        assert not res._cached_path().exists()
        assert _warnings(caplog) == []

    def test_transient_retrieval_failure_is_retried_not_warned(
        self,
        local_config: Config,
        caplog: pytest.LogCaptureFixture,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """The download layer's existing retries cover a flaky first attempt."""
        cfg = dataclasses.replace(local_config, max_retries=3)
        res = self._gene_annotation(cfg)

        real_http_open = _utils_module.http_open
        attempts = {"n": 0}

        def flaky_http_open(url: str, **kwargs: Any) -> Any:
            attempts["n"] += 1
            if attempts["n"] == 1:
                raise urllib.error.URLError("transient")
            return real_http_open(url, **kwargs)

        monkeypatch.setattr(_utils_module, "http_open", flaky_http_open)
        monkeypatch.setattr(_utils_module.time, "sleep", lambda _s: None)

        with caplog.at_level(logging.DEBUG):
            counts = _peek_gtf_feature_counts(res)

        assert attempts["n"] == 2
        assert counts["gene"] == 2000
        assert _warnings(caplog) == []

    def test_persistent_retrieval_failure_is_reported_as_retrieval(
        self,
        local_config: Config,
        caplog: pytest.LogCaptureFixture,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        res = self._gene_annotation(local_config)

        def always_down(*_args: Any, **_kwargs: Any) -> None:
            raise urllib.error.URLError("down")

        monkeypatch.setattr(_utils_module, "http_open", always_down)

        with pytest.raises(DownloadError):
            _peek_gtf_feature_counts(res)

        bundle = R3ResourceBundle(resources=[res])
        with caplog.at_level(logging.DEBUG):
            picked = _select_gtf_resource_for_unit(
                bundle, genomic_unit="gene", annotation_extension="G026"
            )
        # Unusable, but still the best-ranked candidate on offer.
        assert picked is res
        messages = _warnings(caplog)
        assert messages == []  # Selection uses the unambiguous descriptor.

    def test_corrupt_cache_entry_is_reported_as_a_parse_failure(
        self, local_config: Config
    ) -> None:
        res = self._gene_annotation(local_config)
        res.download(path=None, cache_mode="enable")
        res._cached_path().write_bytes(b"this is not gzipped GTF")

        counts_df = pd.DataFrame(
            [[1.0]], index=["ENSG00000278704.1"], columns=["SRR001"]
        )
        res_count = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=counts_df,
            genomic_unit="gene",
            annotation_extension="G026",
        )
        bundle = R3ResourceBundle(resources=[res_count, res])

        with pytest.raises(ValueError, match="could not be parsed"):
            bundle.to_ranged_summarized_experiment(
                genomic_unit="gene", annotation_extension="G026"
            )


class TestClassifyRangesFailure:
    """The failure modes are told apart because the user's fix differs."""

    def test_download_failure_is_retrieval(self) -> None:
        assert (
            _classify_ranges_failure(DownloadError("boom"))
            == "the annotation could not be retrieved"
        )

    def test_absent_file_is_retrieval(self) -> None:
        assert (
            _classify_ranges_failure(FileNotFoundError(2, "nope"))
            == "the annotation could not be retrieved"
        )

    def test_corrupt_archive_is_content_not_retrieval(self) -> None:
        """gzip.BadGzipFile is an OSError, but the bytes did arrive."""
        assert (
            _classify_ranges_failure(gzip.BadGzipFile("not gzipped"))
            == "the annotation could not be parsed"
        )

    def test_parse_failure_is_content(self) -> None:
        assert (
            _classify_ranges_failure(ValueError("bad columns"))
            == "the annotation could not be parsed"
        )

    def test_uncovered_features_is_a_mismatch(self) -> None:
        assert (
            _classify_ranges_failure(RangesCoverageError("missing"))
            == "the annotation does not cover every counted feature"
        )

    def test_nothing_to_read_is_reported_as_absent(self) -> None:
        assert (
            _classify_ranges_failure(MissingRangesError("none"))
            == "no annotation providing genomic ranges was in the bundle"
        )

    def test_no_attempt_at_all_is_reported_as_absent(self) -> None:
        assert (
            _classify_ranges_failure(None)
            == "no annotation providing genomic ranges was in the bundle"
        )

    def test_source_renames_the_subject_of_the_phrase(self) -> None:
        """Junction ranges come from an RR file, not from an annotation."""
        assert (
            _classify_ranges_failure(
                DownloadError("boom"), source="RR coordinate file"
            )
            == "the RR coordinate file could not be retrieved"
        )

    def test_ranges_errors_still_read_as_value_errors(self) -> None:
        """Callers catching ValueError keep working across this split."""
        assert issubclass(RangesCoverageError, ValueError)
        assert issubclass(MissingRangesError, ValueError)
        assert issubclass(RangesError, ValueError)


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
        desc.url_path.return_value = (
            "human/annotations/gene_sums/human.gene_sums.G026.gtf.gz"
        )
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
            f.write('chr1\tref\texon\t1\t100\t.\t+\t.\texon_id "E1"\n')

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

        ranges_df = pd.DataFrame(
            {
                "seqnames": ["chr1"],
                "starts": [1],
                "ends": [100],
                "strand": ["+"],
            }
        )
        with patch.object(
            _utils_module, "get_genomicranges_class", return_value=mock_gr_cls
        ):
            result = _to_genomic_ranges(ranges_df)
        assert result is mock_instance.set_names.return_value
        pd.testing.assert_frame_equal(
            mock_gr_cls.from_pandas.call_args.args[0],
            ranges_df.reset_index(drop=True),
        )


@pytest.mark.requires_biocpy
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

    def test_empty_assay_is_valid(self) -> None:
        result = _construct_summarized_experiment(
            counts_df=pd.DataFrame(),
            row_df=pd.DataFrame(),
            col_df=pd.DataFrame(),
            assay_name="counts",
        )
        assert result.shape == (0, 0)

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


@pytest.mark.requires_biocpy
class TestConstructRangedSummarizedExperiment:
    def _ranges(self, n: int) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "seqnames": ["chr1"] * n,
                "starts": list(range(1, n + 1)),
                "ends": list(range(100, 100 + n)),
                "strand": ["+"] * n,
            }
        )

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

    def test_empty_assay_is_valid_with_empty_ranges(self) -> None:
        result = _construct_ranged_summarized_experiment(
            counts_df=pd.DataFrame(),
            row_df=pd.DataFrame(),
            col_df=pd.DataFrame(),
            ranges_df=pd.DataFrame(
                columns=["seqnames", "starts", "ends", "strand"]
            ),
            assay_name="counts",
        )
        assert result.shape == (0, 0)

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
        ranges = pd.DataFrame(
            {
                "seqnames": [pd.NA, "chr1"],
                "starts": [1, 2],
                "ends": [10, 20],
                "strand": ["+", "+"],
            }
        )
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
        with pytest.raises(
            ValueError, match="not a recognized count-file type"
        ):
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
        return pd.DataFrame(
            {
                "feature_id": fids,
                "seqnames": ["chr1"] * len(fids),
                "starts": list(range(1, len(fids) + 1)),
                "ends": list(range(100, 100 + len(fids))),
                "strand": ["+"] * len(fids),
            }
        )

    def test_missing_feature_id_raises(self) -> None:
        df = pd.DataFrame({"seqnames": ["chr1"]})
        with pytest.raises(ValueError, match="missing required column"):
            _dedupe_ranges_on_feature_id(df)

    def test_no_duplicates_returns_unchanged(self) -> None:
        df = self._ranges(["A", "B"])
        result = _dedupe_ranges_on_feature_id(df)
        assert len(result) == 2

    def test_missing_coord_columns_raises(self) -> None:
        df = pd.DataFrame(
            {"feature_id": ["A", "A"], "seqnames": ["chr1", "chr2"]}
        )
        with pytest.raises(ValueError, match="missing one or more"):
            _dedupe_ranges_on_feature_id(df)

    def test_consistent_duplicates_deduped_with_warning(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        df = pd.DataFrame(
            {
                "feature_id": ["A", "A"],
                "seqnames": ["chr1", "chr1"],
                "starts": [1, 1],
                "ends": [100, 100],
                "strand": ["+", "+"],
            }
        )
        with caplog.at_level(logging.WARNING):
            result = _dedupe_ranges_on_feature_id(df)
        assert len(result) == 1
        assert "duplicate feature_id" in caplog.text

    def test_inconsistent_duplicates_raises(self) -> None:
        df = pd.DataFrame(
            {
                "feature_id": ["A", "A"],
                "seqnames": ["chr1", "chr2"],
                "starts": [1, 999],
                "ends": [100, 1100],
                "strand": ["+", "-"],
            }
        )
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
            R3ResourceBundle.discover(
                organism=[], data_source="sra", project="P1"
            )

    def test_empty_data_source_raises(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            R3ResourceBundle.discover(
                organism="human", data_source=[], project="P1"
            )

    def test_empty_project_raises(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            R3ResourceBundle.discover(
                organism="human", data_source="sra", project=[]
            )

    def test_single_project_sets_identity(self) -> None:
        with patch.object(
            search_module, "search_project_all", side_effect=self._fake_search
        ):
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
        with patch.object(
            search_module, "search_project_all", side_effect=self._fake_search
        ):
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

        with patch.object(
            search_module, "search_project_all", side_effect=same_url_search
        ):
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
        assert all(
            r.description.resource_type == "count_files_gene_or_exon"
            for r, _ in results
        )

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
        assert (
            filtered.resources[0].description.resource_type
            == "count_files_gene_or_exon"
        )

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
            "count_files_gene_or_exon",
            "count_files_junctions",
            "metadata_files",
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
        b = self._bundle_with_types(
            "count_files_gene_or_exon", "metadata_files"
        )
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
        b = R3ResourceBundle(
            organism="human", data_source="sra", project="SRP001"
        )
        org, src, proj = b._resolve_project_identity(None, None, None)
        assert org == "human"
        assert src == "sra"
        assert proj == "SRP001"

    def test_explicit_overrides_stored(self) -> None:
        b = R3ResourceBundle(
            organism="human", data_source="sra", project="SRP001"
        )
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
        b = R3ResourceBundle(
            organism="human", data_source="sra", project="SRP001"
        )
        with pytest.raises(ValueError, match="does not match"):
            b._resolve_project_identity("mouse", None, None)

    def test_conflicting_data_source_raises(self) -> None:
        b = R3ResourceBundle(
            organism="human", data_source="sra", project="SRP001"
        )
        with pytest.raises(ValueError, match="does not match"):
            b._resolve_project_identity(None, "gtex", None)

    def test_conflicting_project_raises(self) -> None:
        b = R3ResourceBundle(
            organism="human", data_source="sra", project="SRP001"
        )
        with pytest.raises(ValueError, match="does not match"):
            b._resolve_project_identity(None, None, "SRP999")


class TestR3ResourceBundleSamples:
    def test_calls_samples_for_project(self) -> None:
        b = R3ResourceBundle(
            organism="human", data_source="sra", project="SRP001"
        )
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
        with pytest.raises(
            CompatibilityError, match="Incompatible count families"
        ):
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
        with pytest.raises(
            CompatibilityError, match="Feature-level incompatibility"
        ):
            b.stack_count_matrices(compat="feature")

    def test_raises_unknown_compat(self) -> None:
        res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=_gene_df(),
            genomic_unit="gene",
        )
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
        res1 = _mock_resource(
            "count_files_gene_or_exon", loaded_data=df1, genomic_unit="gene"
        )
        res2 = _mock_resource(
            "count_files_gene_or_exon", loaded_data=df2, genomic_unit="gene"
        )
        b = R3ResourceBundle(resources=[res1, res2])
        result = b.stack_count_matrices(autoload=False, axis=1)
        assert result.shape[1] == 4


class TestR3ResourceBundleStackCountsFor:
    def test_stacks_gene_resources(self) -> None:
        df = _gene_df()
        res = _mock_resource(
            "count_files_gene_or_exon", loaded_data=df, genomic_unit="gene"
        )
        b = R3ResourceBundle(resources=[res])
        result = b._stack_counts_for(genomic_unit="gene", autoload=False)
        assert isinstance(result, pd.DataFrame)

    def test_stacks_exon_resources(self) -> None:
        df = _gene_df()
        res = _mock_resource(
            "count_files_gene_or_exon", loaded_data=df, genomic_unit="exon"
        )
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
        with pytest.raises(
            ValueError, match="Failed to load requested count matrix"
        ):
            b._stack_counts_for(genomic_unit="gene", autoload=True)

    def test_raises_with_load_errors_junction(self) -> None:
        res = _mock_resource(
            "count_files_junctions",
            junction_type="ALL",
            junction_extension="MM",
        )
        res.load.side_effect = RuntimeError("load failed")
        b = R3ResourceBundle(resources=[res])
        with pytest.raises(
            ValueError, match="Failed to load requested count matrix"
        ):
            b._stack_counts_for(genomic_unit="junction", autoload=True)


class TestNormalizeSampleMetadata:
    def test_no_metadata_returns_external_id_only(self) -> None:
        b = R3ResourceBundle()
        result = b._normalize_sample_metadata(sample_ids=["SRR001", "SRR002"])
        assert "external_id" in result.columns
        assert list(result["external_id"]) == ["SRR001", "SRR002"]

    def test_with_metadata_frames(self) -> None:
        meta_df = pd.DataFrame(
            {
                "external_id": pd.array(["SRR001", "SRR002"], dtype="string"),
                "rail_id": pd.array(["1", "2"], dtype="string"),
                "study": pd.array(["SRP001", "SRP001"], dtype="string"),
                "score": [99.0, 88.0],
            }
        )
        res = _mock_resource(
            "metadata_files", loaded_data=meta_df, table_name="recount_qc"
        )
        b = R3ResourceBundle(resources=[res])
        result = b._normalize_sample_metadata(sample_ids=["SRR001", "SRR002"])
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2

    def test_fills_missing_external_id(self) -> None:
        meta_df = pd.DataFrame(
            {
                "external_id": pd.array([pd.NA, "SRR002"], dtype="string"),
                "rail_id": pd.array(["1", "2"], dtype="string"),
                "study": pd.array(["SRP001", "SRP001"], dtype="string"),
            }
        )
        res = _mock_resource(
            "metadata_files", loaded_data=meta_df, table_name="recount_qc"
        )
        b = R3ResourceBundle(resources=[res])
        result = b._normalize_sample_metadata(
            sample_ids=["1", "SRR002"], metadata_join="outer"
        )
        assert result["external_id"].notna().all()

    def test_provenance_stored_in_attrs(self) -> None:
        meta_df = pd.DataFrame(
            {
                "external_id": pd.array(["SRR001"], dtype="string"),
                "rail_id": pd.array(["1"], dtype="string"),
                "study": pd.array(["SRP001"], dtype="string"),
                "score": [1.0],
            }
        )
        res = _mock_resource(
            "metadata_files", loaded_data=meta_df, table_name="qc_table"
        )
        b = R3ResourceBundle(resources=[res])
        result = b._normalize_sample_metadata(sample_ids=["SRR001"])
        assert "recount3_metadata_provenance" in result.attrs

    @staticmethod
    def _meta_resource(
        external_ids: list[str],
        rail_ids: list[str],
        *,
        table_name: str,
        column: str,
    ) -> MagicMock:
        """Return a loaded metadata resource for the given samples."""
        frame = pd.DataFrame(
            {
                "external_id": pd.array(external_ids, dtype="string"),
                "rail_id": pd.array(rail_ids, dtype="string"),
                "study": pd.array(
                    ["SRP001"] * len(external_ids), dtype="string"
                ),
                column: list(range(len(external_ids))),
            }
        )
        return _mock_resource(
            "metadata_files", loaded_data=frame, table_name=table_name
        )

    def test_inner_join_matching_nothing_raises(self) -> None:
        """An inner join that matches no sample must not be silent.

        The count samples would otherwise all be dropped and the caller
        handed an experiment with no columns at all.
        """
        b = R3ResourceBundle(
            resources=[
                self._meta_resource(
                    ["SRR001", "SRR002"],
                    ["1", "2"],
                    table_name="sra",
                    column="a",
                ),
                self._meta_resource(
                    ["SRR901", "SRR902"],
                    ["901", "902"],
                    table_name="recount_qc",
                    column="b",
                ),
            ]
        )
        with pytest.raises(
            ValueError, match="No count sample matched the merged sample"
        ):
            b._normalize_sample_metadata(sample_ids=["SRR001", "SRR002"])

    def test_outer_join_matching_nothing_keeps_count_samples(self) -> None:
        """The error message's suggested remedy actually works."""
        b = R3ResourceBundle(
            resources=[
                self._meta_resource(
                    ["SRR001", "SRR002"],
                    ["1", "2"],
                    table_name="sra",
                    column="a",
                ),
                self._meta_resource(
                    ["SRR901", "SRR902"],
                    ["901", "902"],
                    table_name="recount_qc",
                    column="b",
                ),
            ]
        )
        result = b._normalize_sample_metadata(
            sample_ids=["SRR001", "SRR002"], metadata_join="outer"
        )
        assert list(result["external_id"]) == ["SRR001", "SRR002"]

    def test_no_sample_ids_does_not_raise(self) -> None:
        """The guard is about lost samples, not about an empty request."""
        b = R3ResourceBundle()
        result = b._normalize_sample_metadata(sample_ids=[])
        assert len(result) == 0


@pytest.mark.requires_biocpy
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

    def test_duplicate_feature_ids_made_unique(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
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

    def test_sample_alignment_losing_every_sample_raises(self) -> None:
        """No code path may hand back an experiment with zero samples.

        `_normalize_sample_metadata` already refuses an inner join that
        matches nothing, so this backstop in `_prepare_experiment` is
        reached by forcing that method to return empty alignment.
        """
        b = self._bundle_with_gene_counts()
        empty = pd.DataFrame({"external_id": pd.array([], dtype="string")})
        empty.attrs["recount3_metadata_provenance"] = {}
        with patch.object(
            R3ResourceBundle,
            "_normalize_sample_metadata",
            return_value=empty,
        ):
            with pytest.raises(
                ValueError, match="Sample alignment produced no samples"
            ):
                b.to_summarized_experiment(genomic_unit="gene", autoload=False)


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
        desc_ann.url_path.return_value = (
            "human/annotations/gene_sums/human.gene_sums.G026.gtf.gz"
        )
        desc_ann.annotation_extension = "G026"
        desc_ann.genomic_unit = "gene"
        res_ann.description = desc_ann
        res_ann.url = "http://example.com/gene_sums/human.gene_sums.G026.gtf.gz"
        res_ann._cached_path.return_value = gtf_path

        return R3ResourceBundle(resources=[res_count, res_ann])

    @pytest.mark.requires_biocpy
    def test_gene_with_gtf_annotation(self, _synthetic_gtf_gz: Path) -> None:
        b = self._bundle_with_gene_gtf(_synthetic_gtf_gz)
        rse = b.to_ranged_summarized_experiment(
            genomic_unit="gene", autoload=False
        )
        assert rse is not None

    @pytest.mark.requires_biocpy
    def test_gene_with_unstranded_gtf_annotation(
        self, _synthetic_gtf_gz: Path
    ) -> None:
        with gzip.open(_synthetic_gtf_gz, "rt", encoding="utf-8") as handle:
            content = handle.read()
        content = content.replace("\t+\t", "\t.\t").replace("\t-\t", "\t.\t")
        with gzip.open(_synthetic_gtf_gz, "wt", encoding="utf-8") as handle:
            handle.write(content)
        b = self._bundle_with_gene_gtf(_synthetic_gtf_gz)
        rse = b.to_ranged_summarized_experiment(
            genomic_unit="gene", autoload=False
        )
        assert rse.shape == (3, 2)
        assert list(rse.row_ranges.strand) == [0, 0, 0]

    @pytest.mark.requires_biocpy
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

    @pytest.mark.requires_biocpy
    def test_junction_with_rr_coordinates(self, tmp_path: Path) -> None:
        rr_content = (
            "seqnames\tstarts\tends\tstrand\tjunction_id\n"
            "chr1\t1\t100\t+\tJX001\n"
            "chr1\t200\t300\t+\tJX002\n"
        )
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

    @pytest.mark.requires_biocpy
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

    @pytest.mark.requires_biocpy
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

    @pytest.mark.requires_biocpy
    def test_rr_without_strand_column_defaults_to_star(
        self, tmp_path: Path
    ) -> None:
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

    @pytest.mark.requires_biocpy
    def test_duplicate_feature_ids_in_counts(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
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

    def test_parallel_calls_download_on_each_resource(self) -> None:
        resources = [_mock_resource() for _ in range(5)]
        b = R3ResourceBundle(resources=resources)
        b.download(dest="/some/dir", max_workers=4)
        for res in resources:
            res.download.assert_called_once_with(
                path="/some/dir", cache_mode="enable", overwrite=False
            )

    def test_max_workers_one_is_sequential(self) -> None:
        res1 = _mock_resource()
        res2 = _mock_resource()
        b = R3ResourceBundle(resources=[res1, res2])
        b.download(dest="/some/dir", max_workers=1)
        res1.download.assert_called_once_with(
            path="/some/dir", cache_mode="enable", overwrite=False
        )
        res2.download.assert_called_once_with(
            path="/some/dir", cache_mode="enable", overwrite=False
        )

    def test_download_error_propagates(self) -> None:
        good = _mock_resource()
        bad = _mock_resource()
        bad.download.side_effect = RuntimeError("boom")
        b = R3ResourceBundle(resources=[good, bad])
        with pytest.raises(RuntimeError, match="boom"):
            b.download(max_workers=4)

    def test_empty_bundle_parallel_noop(self) -> None:
        b = R3ResourceBundle()
        b.download(max_workers=4)

    def test_parallel_cached_transfers_overlap_and_fill_one_zip(
        self, local_config: Config, tmp_path: Path
    ) -> None:
        """Distinct cached transfers overlap and produce a complete ZIP.

        The barrier only trips once every transfer is in flight, so a bundle
        that serialized its downloads would hang instead of finishing. The
        archive is then written under a single writer lock, which the ZIP's
        own integrity check confirms.
        """
        barrier = threading.Barrier(4, timeout=5)
        resources = [
            R3Resource(
                R3ResourceDescription(
                    resource_type="count_files_gene_or_exon",
                    organism="human",
                    data_source="sra",
                    project=f"SRP00{index}",
                    genomic_unit="gene",
                    annotation_extension="G026",
                ),
                config=local_config,
            )
            for index in range(4)
        ]

        def transfer(url: str, path: Path, **_: Any) -> None:
            barrier.wait()
            path.write_bytes(url.encode())

        target = tmp_path / "bundle.zip"
        with patch.object(resource_module, "download_to_file", transfer):
            R3ResourceBundle(resources=resources).download(
                dest=str(target), max_workers=4
            )

        with zipfile.ZipFile(target) as archive:
            assert archive.testzip() is None
            assert set(archive.namelist()) == {res.arcname for res in resources}
            for res in resources:
                assert archive.read(res.arcname) == res.url.encode()


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
        gtf = pd.DataFrame(
            [
                {
                    "seqname": "chr1",
                    "source": "ref",
                    "feature": "exon",
                    "start": 1,
                    "end": 100,
                    "score": ".",
                    "strand": "+",
                    "frame": ".",
                    "attributes": ".",
                }
            ]
        )
        result = _ranges_from_gtf(gtf, feature_kind="exon")
        assert len(result) == 1
        assert "chr1" in str(result["feature_id"].iloc[0])

    def test_gene_with_empty_attrs_uses_coord_fallback(self) -> None:
        gtf = pd.DataFrame(
            [
                {
                    "seqname": "chrX",
                    "source": "ref",
                    "feature": "gene",
                    "start": 500,
                    "end": 600,
                    "score": ".",
                    "strand": "-",
                    "frame": ".",
                    "attributes": ".",
                }
            ]
        )
        result = _ranges_from_gtf(gtf, feature_kind="gene")
        assert len(result) == 1

    def test_extra_attrs_joined_to_output(self) -> None:
        gtf = pd.DataFrame(
            [
                {
                    "seqname": "chr1",
                    "source": "ref",
                    "feature": "gene",
                    "start": 1,
                    "end": 100,
                    "score": ".",
                    "strand": "+",
                    "frame": ".",
                    "attributes": 'gene_id "G1"; gene_name "MYC"',
                }
            ]
        )
        result = _ranges_from_gtf(gtf, feature_kind="gene")
        assert "gene_name" in result.columns


class TestPeekGtfShortLines:
    def test_line_with_two_fields_skipped(self, tmp_path: Path) -> None:
        gz_path = tmp_path / "short.gtf.gz"
        content = (
            "chr1\tref\n" 'chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "G1"\n'
        )
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write(content)
        res = MagicMock(spec=R3Resource)
        res.ensure_cached.return_value = gz_path
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

    def test_resource_url_without_gtf_still_selected(
        self, tmp_path: Path
    ) -> None:
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

        with patch.object(
            search_module, "search_project_all", side_effect=search_side_effect
        ):
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
        gene_res = _mock_resource(
            "count_files_gene_or_exon", loaded_data=gene_df, genomic_unit="gene"
        )
        meta_res = _mock_resource(
            "metadata_files", loaded_data=pd.DataFrame({"a": [1]})
        )
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
    def test_non_dataframe_metadata_is_an_error(self) -> None:
        res = _mock_resource(
            "metadata_files", loaded_data="bad table", table_name="sra"
        )
        with pytest.raises(TypeError, match="not a DataFrame"):
            R3ResourceBundle([res])._normalize_sample_metadata(
                sample_ids=["SRR001"]
            )


class TestNormalizeSampleMetadataExternalIdFillback:
    def test_aligned_without_external_id_gets_filled(self) -> None:
        meta_df = pd.DataFrame(
            {
                "rail_id": pd.array(["1", "2"], dtype="string"),
                "study": pd.array(["SRP001", "SRP001"], dtype="string"),
                "score": [10.0, 20.0],
            }
        )
        res = _mock_resource(
            "metadata_files", loaded_data=meta_df, table_name="qc"
        )
        b = R3ResourceBundle(resources=[res])
        result = b._normalize_sample_metadata(
            sample_ids=["1", "2"], metadata_join="outer"
        )
        assert "external_id" in result.columns


class TestToRangedSEAutoload:
    @pytest.mark.requires_biocpy
    def test_autoload_true_reuses_cached_gtf(self, tmp_path: Path) -> None:
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
        desc_ann.url_path.return_value = (
            "human/annotations/gene_sums/human.gene_sums.G026.gtf.gz"
        )
        desc_ann.annotation_extension = "G026"
        desc_ann.genomic_unit = "gene"
        res_ann.description = desc_ann
        res_ann.url = "http://example.com/gene_sums/human.gene_sums.G026.gtf.gz"
        res_ann._cached_path.return_value = gz_path

        b = R3ResourceBundle(resources=[res_count, res_ann])
        rse = b.to_ranged_summarized_experiment(
            genomic_unit="gene", autoload=True
        )
        # The annotation is already in the cache, so no download is issued.
        res_ann.download.assert_not_called()
        assert rse is not None

    @pytest.mark.requires_biocpy
    def test_autoload_becomes_the_resources_download_policy(
        self, tmp_path: Path
    ) -> None:
        """autoload reaches the annotation, which is what the bug missed.

        Whether an absent file is then fetched is R3Resource.ensure_cached's
        job; see TestEnsureCached in test_resource.py.
        """
        gz_path = tmp_path / "cache" / "genes.gtf.gz"
        content = (
            'chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG001"\n'
            'chr2\tref\tgene\t200\t300\t.\t-\t.\tgene_id "ENSG002"\n'
        )

        counts_df = pd.DataFrame(
            np.ones((2, 1), dtype=float),
            index=["ENSG001", "ENSG002"],
            columns=["SRR001"],
        )
        res_count = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=counts_df,
            genomic_unit="gene",
            annotation_extension="G026",
        )

        res_ann = MagicMock(spec=R3Resource)
        desc_ann = MagicMock()
        desc_ann.resource_type = "annotations"
        desc_ann.url_path.return_value = "human/ann/gene.gtf.gz"
        desc_ann.annotation_extension = "G026"
        desc_ann.genomic_unit = "gene"
        res_ann.description = desc_ann
        res_ann.url = "http://example.com/ann/gene.gtf.gz"
        res_ann._cached_path.return_value = gz_path

        def fake_ensure_cached(*, download: bool) -> Path:
            assert download is True
            gz_path.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(gz_path, "wt", encoding="utf-8") as fh:
                fh.write(content)
            return gz_path

        res_ann.ensure_cached.side_effect = fake_ensure_cached

        b = R3ResourceBundle(resources=[res_count, res_ann])
        rse = b.to_ranged_summarized_experiment(
            genomic_unit="gene", autoload=True
        )
        res_ann.ensure_cached.assert_called_with(download=True)
        assert rse is not None

    @pytest.mark.requires_biocpy
    def test_enrich_cols_joined_to_row_data(self, tmp_path: Path) -> None:
        gz_path = tmp_path / "genes_with_name.gtf.gz"
        content = 'chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG001"; gene_name "MYC"\n'
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write(content)

        counts_df = pd.DataFrame([[5.0]], index=["ENSG001"], columns=["SRR001"])
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
        rse = b.to_ranged_summarized_experiment(
            genomic_unit="gene", autoload=False
        )
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

    @pytest.mark.requires_biocpy
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
        res_rr.download.assert_not_called()  # Already-loaded RR is reused.
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
            assert (
                gene_id in ranges["feature_id"].values
            ), f"{gene_id} not found in GTF ranges"

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

    @pytest.mark.requires_biocpy
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
        rse = b.to_ranged_summarized_experiment(
            genomic_unit="gene", autoload=False
        )

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

        with patch.object(
            bmod, "_peek_gtf_feature_counts", side_effect=OSError("peek fail")
        ):
            result = _select_gtf_resource_for_unit(
                bundle, genomic_unit="gene", annotation_extension=None
            )
        assert result is res


class TestAddBigwigUrls:
    def test_no_external_id_column_returns_na(self) -> None:
        b = R3ResourceBundle()
        col_df = pd.DataFrame({"other": ["SRR001"]})
        result = b._add_bigwig_urls(col_df)
        assert "BigWigURL" in result.columns
        assert pd.isna(result["BigWigURL"].iloc[0])

    def test_no_count_resources_sets_bigwig_na(self) -> None:
        b = R3ResourceBundle()
        col_df = pd.DataFrame(
            {"external_id": pd.array(["SRR001"], dtype="string")}
        )
        result = b._add_bigwig_urls(col_df)
        assert "BigWigURL" in result.columns
        assert pd.isna(result["BigWigURL"].iloc[0])

    def test_non_count_resource_skipped_then_returns_na(self) -> None:
        meta_res = _mock_resource("metadata_files")
        b = R3ResourceBundle(resources=[meta_res])
        col_df = pd.DataFrame(
            {"external_id": pd.array(["SRR001"], dtype="string")}
        )
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
        col_df = pd.DataFrame(
            {"external_id": pd.array(["SRR001"], dtype="string")}
        )
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
        col_df = pd.DataFrame(
            {"external_id": pd.array([pd.NA], dtype="string")}
        )
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
        col_df = pd.DataFrame(
            {
                "external_id": pd.array(["SRR001"], dtype="string"),
                "recount_seq_qc__file_source": ["some/path/sra"],
            }
        )
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
        col_df = pd.DataFrame(
            {
                "external_id": pd.array([pd.NA], dtype="string"),
                "recount_seq_qc__file_source": [42],
            }
        )
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
        col_df = pd.DataFrame(
            {
                "external_id": pd.array([pd.NA], dtype="string"),
                "recount_seq_qc__file_source": ["/"],
            }
        )
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
        col_df = pd.DataFrame(
            {"external_id": pd.array(["SRR001"], dtype="string")}
        )
        mock_cfg = MagicMock()
        mock_cfg.base_url = "https://duffel.example.com/recount3/"
        with patch("recount3.resource.default_config", return_value=mock_cfg):
            result = b._add_bigwig_urls(col_df)
        assert "BigWigURL" in result.columns
        url = result["BigWigURL"].iloc[0]
        assert url is not None
        assert url.startswith("https://duffel.example.com/recount3/")


def _annotation_resource(
    *,
    cached_path: Path,
    url_path: str,
    url: str,
    genomic_unit: Any = None,
    annotation_extension: str = "G026",
) -> MagicMock:
    """Return an annotation resource whose cache lookups hit a real file.

    ``_select_gtf_resource_for_unit`` inspects candidates by reading them, so
    ``ensure_cached`` has to yield a path that actually exists; a bare
    ``MagicMock`` return value silently turns into a file descriptor.
    """
    res = MagicMock(spec=R3Resource)
    desc = MagicMock()
    desc.resource_type = "annotations"
    desc.url_path.return_value = url_path
    desc.annotation_extension = annotation_extension
    desc.genomic_unit = genomic_unit
    res.description = desc
    res.url = url
    res.ensure_cached.return_value = cached_path
    res._cached_path.return_value = cached_path
    return res


def _write_gtf(path: Path, *, feature: str = "gene") -> Path:
    """Write a one-row uncompressed GTF naming ``feature`` and return its path."""
    path.write_text(
        f'chr1\tref\t{feature}\t1\t100\t.\t+\t.\t{feature}_id "F1"\n',
        encoding="utf-8",
    )
    return path


def _junction_resources(
    *,
    project: str,
    samples: list[str],
    rr_frame: pd.DataFrame,
    features: list[str] | None = None,
) -> tuple[MagicMock, MagicMock]:
    """Return a loaded MM/RR junction resource pair for one project."""
    counts = pd.DataFrame(
        np.ones((len(rr_frame), len(samples)), dtype=float),
        index=features or [str(i) for i in range(len(rr_frame))],
        columns=samples,
    )
    shared = {
        "junction_type": "ALL",
        "organism": "human",
        "data_source": "sra",
        "project": project,
    }
    mm = _mock_resource(
        "count_files_junctions",
        url=f"http://example.com/{project}.MM.gz",
        loaded_data=counts,
        junction_extension="MM",
        **shared,
    )
    rr = _mock_resource(
        "count_files_junctions",
        url=f"http://example.com/{project}.RR.gz",
        loaded_data=rr_frame,
        junction_extension="RR",
        **shared,
    )
    return mm, rr


def _rr_frame() -> pd.DataFrame:
    """Return a two-row RR sidecar using R's column spellings."""
    return pd.DataFrame(
        {
            "chromosome": ["chr1", "chr1"],
            "start": [1, 200],
            "end": [100, 300],
            "strand": ["+", "-"],
        }
    )


class TestEnsureUniqueColumnsSuffixCollision:
    def test_generated_suffix_skips_a_name_already_in_use(self) -> None:
        """A ``__2`` suffix must not collide with a literal ``__2`` column."""
        df = pd.DataFrame([[1, 2, 3]], columns=["a", "a", "a__2"])
        out = _ensure_unique_columns(df)
        assert list(out.columns) == ["a", "a__3", "a__2"]


class TestAlignRangesToFeaturesWithDuplicates:
    @staticmethod
    def _ranges(feature_ids: list[str]) -> pd.DataFrame:
        n = len(feature_ids)
        return pd.DataFrame(
            {
                "feature_id": feature_ids,
                "seqnames": [f"chr{i + 1}" for i in range(n)],
                "starts": [(i + 1) * 10 for i in range(n)],
                "ends": [(i + 1) * 10 + 5 for i in range(n)],
                "strand": ["+"] * n,
            }
        )

    def test_repeated_ids_are_matched_by_occurrence_order(self) -> None:
        ranges = self._ranges(["G1", "G1", "G2"])
        out = _align_ranges_to_features(ranges, feature_ids=["G2", "G1", "G1"])
        assert list(out.index) == ["G2", "G1", "G1"]
        # The two G1 rows keep their original order after reordering.
        assert list(out["starts"]) == [30, 10, 20]

    def test_differing_duplicate_counts_are_ambiguous(self) -> None:
        ranges = self._ranges(["G1", "G1"])
        with pytest.raises(
            RangesCoverageError, match="Ambiguous duplicate feature occurrences"
        ):
            _align_ranges_to_features(ranges, feature_ids=["G1"])

    def test_unannotated_feature_alongside_duplicates_is_rejected(self) -> None:
        ranges = self._ranges(["G1", "G1"])
        with pytest.raises(
            RangesCoverageError, match="Cannot match duplicate feature"
        ):
            _align_ranges_to_features(ranges, feature_ids=["G1", "G1", "G2"])


class TestSelectGtfResourceForUnitRanking:
    def test_two_exact_unit_matches_are_ambiguous(self, tmp_path: Path) -> None:
        path = _write_gtf(tmp_path / "genes.gtf")
        resources = [
            _annotation_resource(
                cached_path=path,
                url_path=f"ann/{name}.gtf.gz",
                url=f"http://example.com/{name}.gtf.gz",
                genomic_unit="gene",
            )
            for name in ("first", "second")
        ]
        bundle = R3ResourceBundle(resources=resources)
        with pytest.raises(
            CompatibilityError, match="Multiple matching annotations"
        ):
            _select_gtf_resource_for_unit(
                bundle, genomic_unit="gene", annotation_extension=None
            )

    def test_case_insensitive_unit_outranks_an_unlabelled_candidate(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """``"Gene"`` misses the exact match but still scores highest."""
        labelled = _annotation_resource(
            cached_path=_write_gtf(tmp_path / "a.gtf"),
            url_path="ann/a.txt",
            url="http://example.com/a.txt",
            genomic_unit="Gene",
        )
        unlabelled = _annotation_resource(
            cached_path=_write_gtf(tmp_path / "b.gtf"),
            url_path="ann/b.txt",
            url="http://example.com/b.txt",
            genomic_unit=None,
        )
        bundle = R3ResourceBundle(resources=[unlabelled, labelled])
        with caplog.at_level(logging.INFO):
            result = _select_gtf_resource_for_unit(
                bundle, genomic_unit="gene", annotation_extension=None
            )
        assert result is labelled
        assert "Selected annotation resource for gene" in caplog.text

    def test_exon_named_gtf_scores_on_name_and_extension(
        self, tmp_path: Path
    ) -> None:
        res = _annotation_resource(
            cached_path=_write_gtf(tmp_path / "e.gtf", feature="exon"),
            url_path="ann/human.exon_sums.G026.gtf.gz",
            url="http://example.com/human.exon_sums.G026.gtf.gz",
            genomic_unit=None,
        )
        bundle = R3ResourceBundle(resources=[res])
        result = _select_gtf_resource_for_unit(
            bundle, genomic_unit="exon", annotation_extension=None
        )
        assert result is res

    def test_uncached_candidate_is_skipped_when_autoload_is_off(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        res = _annotation_resource(
            cached_path=tmp_path / "absent.gz",
            url_path="ann/x.gz",
            url="http://example.com/x.gz",
            genomic_unit=None,
        )
        res.ensure_cached.side_effect = FileNotFoundError("not cached")
        bundle = R3ResourceBundle(resources=[res])
        with caplog.at_level(logging.DEBUG):
            result = _select_gtf_resource_for_unit(
                bundle,
                genomic_unit="gene",
                annotation_extension=None,
                autoload=False,
            )
        # Nothing could be inspected, so the ranking alone decides.
        assert result is res
        assert "autoload is disabled" in caplog.text

    def test_unreadable_candidate_is_logged_and_skipped(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        corrupt = tmp_path / "bad.gtf.gz"
        corrupt.write_text("this is not gzip", encoding="utf-8")
        higher_scoring = _annotation_resource(
            cached_path=corrupt,
            url_path="ann/gene_sums.gtf.gz",
            url="http://example.com/gene_sums.gtf.gz",
            genomic_unit=None,
        )
        readable = _annotation_resource(
            cached_path=_write_gtf(tmp_path / "good.gtf"),
            url_path="ann/g.txt",
            url="http://example.com/g.txt",
            genomic_unit=None,
        )
        bundle = R3ResourceBundle(resources=[higher_scoring, readable])
        with caplog.at_level(logging.WARNING):
            result = _select_gtf_resource_for_unit(
                bundle, genomic_unit="gene", annotation_extension=None
            )
        assert result is readable
        assert "Skipping annotation candidate" in caplog.text


class TestNumericAssayValidation:
    def test_sparse_counts_need_a_zero_fill_value(self) -> None:
        counts = pd.DataFrame(
            {"s1": pd.arrays.SparseArray([1.0, 2.0], fill_value=1.0)},
            index=["g1", "g2"],
        )
        with pytest.raises(ValueError, match="zero fill values"):
            bmod._numeric_assay(counts)

    def test_masked_columns_with_missing_counts_are_rejected(self) -> None:
        """Mixing masked and plain integer columns yields an object matrix."""
        counts = pd.DataFrame(
            {
                "s1": pd.array([1, pd.NA], dtype="Int64"),
                "s2": np.array([3, 4]),
            },
            index=["g1", "g2"],
        )
        with pytest.raises(ValueError, match="missing or non-finite"):
            bmod._numeric_assay(counts)

    def test_negative_counts_are_rejected(self) -> None:
        counts = pd.DataFrame({"s1": [1.0, -2.0]}, index=["g1", "g2"])
        with pytest.raises(ValueError, match="negative counts"):
            bmod._numeric_assay(counts)


class TestValidateCountFrame:
    def test_sparse_column_needs_a_zero_fill_value(self) -> None:
        frame = pd.DataFrame(
            {"s1": pd.arrays.SparseArray([1.0, 2.0], fill_value=3.0)}
        )
        with pytest.raises(ValueError, match="zero fill values"):
            bmod._validate_count_frame(frame)

    def test_masked_boolean_column_with_missing_values_is_rejected(
        self,
    ) -> None:
        frame = pd.DataFrame({"s1": pd.array([True, pd.NA], dtype="boolean")})
        with pytest.raises(
            ValueError, match="finite, non-missing and nonnegative"
        ):
            bmod._validate_count_frame(frame)


@pytest.mark.requires_biocpy
class TestExperimentFrameAlignment:
    def test_rows_are_reordered_to_match_the_assay(self) -> None:
        frame = pd.DataFrame({"score": [2.0, 1.0]}, index=["s2", "s1"])
        out = bmod._experiment_frame(frame, ["s1", "s2"], "column_data")
        assert list(out.get_column("score")) == [1.0, 2.0]
        assert list(out.get_row_names()) == ["s1", "s2"]


@pytest.mark.requires_biocpy
class TestConstructRangedSummarizedExperimentRangeAlignment:
    @staticmethod
    def _ranges(index: list[str]) -> pd.DataFrame:
        n = len(index)
        return pd.DataFrame(
            {
                "seqnames": ["chr1"] * n,
                "starts": [(i + 1) * 10 for i in range(n)],
                "ends": [(i + 1) * 10 + 5 for i in range(n)],
                "strand": ["+"] * n,
            },
            index=index,
        )

    def test_ranges_are_reindexed_to_the_assay_order(self) -> None:
        counts = _gene_df()
        ranges = self._ranges(["ENSG0002", "ENSG0001"])
        rse = _construct_ranged_summarized_experiment(
            counts_df=counts,
            row_df=pd.DataFrame({"a": [1, 2]}),
            col_df=pd.DataFrame({"b": ["x", "y"]}),
            ranges_df=ranges,
            assay_name="raw",
        )
        assert list(rse.get_row_ranges().get_start()) == [20, 10]

    def test_ranges_identifiers_must_match_the_assay(self) -> None:
        counts = _gene_df()
        ranges = self._ranges(["other1", "other2"])
        with pytest.raises(
            ValueError, match="ranges_df identifiers do not match"
        ):
            _construct_ranged_summarized_experiment(
                counts_df=counts,
                row_df=pd.DataFrame({"a": [1, 2]}),
                col_df=pd.DataFrame({"b": ["x", "y"]}),
                ranges_df=ranges,
                assay_name="raw",
            )


class TestValidateCoordinates:
    @staticmethod
    def _frame(**overrides: Any) -> pd.DataFrame:
        values: dict[str, Any] = {
            "seqnames": ["chr1"],
            "starts": [10],
            "ends": [20],
            "strand": ["+"],
        }
        values.update(overrides)
        return pd.DataFrame(values)

    def test_fractional_coordinates_are_rejected(self) -> None:
        with pytest.raises(ValueError, match="finite integers"):
            bmod._validate_coordinates(self._frame(starts=[10.5]))

    def test_starts_below_one_are_rejected(self) -> None:
        with pytest.raises(ValueError, match="positive inclusive intervals"):
            bmod._validate_coordinates(self._frame(starts=[0]))

    def test_end_before_start_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="positive inclusive intervals"):
            bmod._validate_coordinates(self._frame(starts=[20], ends=[10]))

    def test_unknown_strand_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="Invalid genomic strand"):
            bmod._validate_coordinates(self._frame(strand=["?"]))


class TestMergeCountFrames:
    def test_rejects_an_unknown_join_policy(self) -> None:
        with pytest.raises(ValueError, match="join_policy must be"):
            bmod._merge_count_frames([_gene_df()], "left")

    def test_rejects_duplicate_feature_ids_across_projects(self) -> None:
        first = _gene_df(features=["G1", "G1"], samples=["S1"])
        second = _gene_df(features=["G2"], samples=["S2"])
        with pytest.raises(
            CompatibilityError, match="Cannot align duplicate feature IDs"
        ):
            bmod._merge_count_frames([first, second], "inner")


class TestStackCountsForValidation:
    def test_rejects_an_unknown_join_policy(self) -> None:
        bundle = R3ResourceBundle()
        with pytest.raises(ValueError, match="join_policy must be"):
            bundle._stack_counts_for(genomic_unit="gene", join_policy="left")

    def test_rejects_mixed_organisms(self) -> None:
        resources = [
            _mock_resource(
                "count_files_gene_or_exon",
                loaded_data=_gene_df(),
                genomic_unit="gene",
                organism=organism,
                annotation_extension="G026",
            )
            for organism in ("human", "mouse")
        ]
        bundle = R3ResourceBundle(resources=resources)
        with pytest.raises(
            CompatibilityError, match="incompatible organisms, annotations"
        ):
            bundle._stack_counts_for(genomic_unit="gene", autoload=False)

    def test_unloaded_resource_without_autoload_is_an_error(self) -> None:
        res = _mock_resource("count_files_gene_or_exon", genomic_unit="gene")
        bundle = R3ResourceBundle(resources=[res])
        with pytest.raises(ValueError, match="Count resource is not loaded"):
            bundle._stack_counts_for(genomic_unit="gene", autoload=False)

    def test_non_dataframe_counts_are_an_error(self) -> None:
        res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data="not a frame",
            genomic_unit="gene",
        )
        bundle = R3ResourceBundle(resources=[res])
        with pytest.raises(
            TypeError, match="Loaded counts are not a DataFrame"
        ):
            bundle._stack_counts_for(genomic_unit="gene", autoload=False)

    def test_duplicate_sample_identifiers_are_an_error(self) -> None:
        counts = pd.DataFrame(
            np.ones((2, 2)), index=["G1", "G2"], columns=["S1", "S1"]
        )
        res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=counts,
            genomic_unit="gene",
        )
        bundle = R3ResourceBundle(resources=[res])
        with pytest.raises(ValueError, match="Duplicate sample identifiers in"):
            bundle._stack_counts_for(genomic_unit="gene", autoload=False)

    def test_blank_sample_identifiers_are_an_error(self) -> None:
        counts = pd.DataFrame(
            np.ones((2, 2)), index=["G1", "G2"], columns=["S1", "   "]
        )
        res = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=counts,
            genomic_unit="gene",
        )
        bundle = R3ResourceBundle(resources=[res])
        with pytest.raises(
            ValueError, match="Missing feature or sample identifiers"
        ):
            bundle._stack_counts_for(genomic_unit="gene", autoload=False)

    def test_multiple_junction_projects_are_keyed_by_rr_coordinates(
        self,
    ) -> None:
        """Stacking two MM matrices relabels both to their RR coordinates."""
        mm1, rr1 = _junction_resources(
            project="SRP001", samples=["S1", "S2"], rr_frame=_rr_frame()
        )
        mm2, rr2 = _junction_resources(
            project="SRP002", samples=["S3", "S4"], rr_frame=_rr_frame()
        )
        bundle = R3ResourceBundle(resources=[mm1, rr1, mm2, rr2])
        out = bundle._stack_counts_for(genomic_unit="junction", autoload=False)
        assert list(out.index) == ["chr1:1-100:+", "chr1:200-300:-"]
        assert list(out.columns) == ["S1", "S2", "S3", "S4"]


class TestJunctionRanges:
    def test_requires_exactly_one_mm_matrix(self) -> None:
        mm1, rr1 = _junction_resources(
            project="SRP001", samples=["S1"], rr_frame=_rr_frame()
        )
        mm2, rr2 = _junction_resources(
            project="SRP002", samples=["S2"], rr_frame=_rr_frame()
        )
        bundle = R3ResourceBundle(resources=[mm1, rr1, mm2, rr2])
        with pytest.raises(
            CompatibilityError, match="needs its own project and RR sidecar"
        ):
            bundle._junction_ranges(mm1.get_loaded(), autoload=False)

    def test_uncached_rr_sidecar_is_read_from_disk(
        self, tmp_path: Path
    ) -> None:
        rr_path = tmp_path / "jxn.RR.tsv"
        rr_path.write_text(
            "chromosome\tstart\tend\tstrand\n"
            "chr1\t1\t100\t+\n"
            "chr1\t200\t300\t-\n",
            encoding="utf-8",
        )
        mm, rr = _junction_resources(
            project="SRP001", samples=["S1"], rr_frame=_rr_frame()
        )
        rr.is_loaded.return_value = False
        rr.get_loaded.return_value = None
        rr._cached_path.return_value = rr_path
        bundle = R3ResourceBundle(resources=[mm, rr])

        ranges = bundle._junction_ranges(mm.get_loaded(), autoload=True)

        rr.ensure_cached.assert_called_once_with(download=True)
        assert list(ranges.index) == ["chr1:1-100:+", "chr1:200-300:-"]
        # The parsed table is cached on the resource for later reuse.
        pd.testing.assert_frame_equal(rr._cached_data, pd.read_table(rr_path))

    def test_duplicate_rr_coordinates_are_rejected(self) -> None:
        duplicated = pd.DataFrame(
            {
                "chromosome": ["chr1", "chr1"],
                "start": [1, 1],
                "end": [100, 100],
                "strand": ["+", "+"],
            }
        )
        mm, rr = _junction_resources(
            project="SRP001", samples=["S1"], rr_frame=duplicated
        )
        bundle = R3ResourceBundle(resources=[mm, rr])
        with pytest.raises(
            RangesCoverageError, match="duplicate junction coordinates"
        ):
            bundle._junction_ranges(mm.get_loaded(), autoload=False)


def _metadata_resource(
    frame: pd.DataFrame, *, table_name: str = "recount_qc", **desc: Any
) -> MagicMock:
    """Return a loaded metadata resource wrapping ``frame``."""
    return _mock_resource(
        "metadata_files", loaded_data=frame, table_name=table_name, **desc
    )


class TestUnpublishedMetadataTables:
    """recount3 does not publish every metadata table for every project.

    R treats an unreachable table as expected: ``file_retrieve()`` warns and
    yields ``NA``, ``read_metadata()`` drops it, and the RSE is still built
    from the tables that do exist. These tests pin the same rule here.
    """

    _PROJECT = {
        "organism": "human",
        "data_source": "sra",
        "project": "SRP001",
    }

    @staticmethod
    def _unretrievable(table_name: str, **desc: Any) -> MagicMock:
        res = _mock_resource(
            "metadata_files",
            url=f"http://example.com/{table_name}.MD.gz",
            table_name=table_name,
            **desc,
        )
        res.load.side_effect = DownloadError(f"Failed to download {table_name}")
        return res

    @staticmethod
    def _usable(table_name: str, **desc: Any) -> MagicMock:
        frame = pd.DataFrame(
            {
                "external_id": pd.array(["SRR001", "SRR002"], dtype="string"),
                "rail_id": pd.array(["1", "2"], dtype="string"),
                "study": pd.array(["SRP001", "SRP001"], dtype="string"),
                "score": [1.0, 2.0],
            }
        )
        return _metadata_resource(frame, table_name=table_name, **desc)

    def test_helper_drops_an_unretrievable_table_with_a_warning(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        res = self._unretrievable("recount_pred")
        with caplog.at_level(logging.WARNING):
            assert _load_optional_metadata(res, autoload=True) is None
        assert "could not be retrieved" in caplog.text
        assert "recount_pred" in caplog.text

    def test_helper_propagates_a_parse_failure(self) -> None:
        """A retrieved-but-damaged table is corruption, not absence."""
        res = _mock_resource("metadata_files", table_name="recount_qc")
        res.load.side_effect = LoadError("not a TSV")
        with pytest.raises(LoadError):
            _load_optional_metadata(res, autoload=True)

    def test_helper_requires_loading_when_autoload_is_off(self) -> None:
        res = self._unretrievable("recount_qc")
        with pytest.raises(ValueError, match="Metadata resource is not loaded"):
            _load_optional_metadata(res, autoload=False)
        res.load.assert_not_called()

    def test_an_absent_table_still_yields_the_remaining_metadata(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        bundle = R3ResourceBundle(
            resources=[
                self._usable("recount_qc", **self._PROJECT),
                self._unretrievable("recount_pred", **self._PROJECT),
            ]
        )
        with caplog.at_level(logging.WARNING):
            col = bundle._normalize_sample_metadata(
                sample_ids=["SRR001", "SRR002"]
            )
        assert list(col.index) == ["SRR001", "SRR002"]
        assert list(col["recount_qc__score"]) == [1.0, 2.0]
        assert "could not be retrieved" in caplog.text

    def test_an_absent_table_still_builds_the_experiment(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        counts = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=_gene_df(),
            genomic_unit="gene",
            annotation_extension="G026",
            **self._PROJECT,
        )
        bundle = R3ResourceBundle(
            resources=[
                counts,
                self._usable("recount_qc", **self._PROJECT),
                self._unretrievable("recount_pred", **self._PROJECT),
            ]
        )
        with caplog.at_level(logging.WARNING):
            _, _, col, _, _ = bundle._prepare_experiment(
                genomic_unit="gene",
                annotation_extension=None,
                join_policy="inner",
                metadata_join="inner",
                autoload=True,
            )
        assert list(col.index) == ["SRR001", "SRR002"]
        assert "recount_pred" in caplog.text

    def test_every_table_absent_is_reported_rather_than_ignored(self) -> None:
        """R stops here too, rather than returning bare sample IDs."""
        bundle = R3ResourceBundle(
            resources=[
                self._unretrievable("recount_qc", **self._PROJECT),
                self._unretrievable("recount_pred", **self._PROJECT),
            ]
        )
        with pytest.raises(ValueError, match="empty or could not be retrieved"):
            bundle._normalize_sample_metadata(sample_ids=["SRR001"])


class TestNormalizeSampleMetadataValidation:
    def test_rejects_an_unknown_metadata_join(self) -> None:
        bundle = R3ResourceBundle()
        with pytest.raises(ValueError, match="metadata_join must be"):
            bundle._normalize_sample_metadata(
                sample_ids=["S1"], metadata_join="left"
            )

    def test_unloaded_metadata_without_autoload_is_an_error(self) -> None:
        res = _mock_resource("metadata_files", table_name="recount_qc")
        bundle = R3ResourceBundle(resources=[res])
        with pytest.raises(ValueError, match="Metadata resource is not loaded"):
            bundle._normalize_sample_metadata(sample_ids=["S1"], autoload=False)

    def test_unloaded_metadata_is_loaded_when_autoload_is_on(self) -> None:
        frame = pd.DataFrame(
            {
                "external_id": pd.array(["S1"], dtype="string"),
                "rail_id": pd.array(["1"], dtype="string"),
                "study": pd.array(["SRP001"], dtype="string"),
                "score": [1.0],
            }
        )
        res = _metadata_resource(frame)
        res.is_loaded.return_value = False
        bundle = R3ResourceBundle(resources=[res])

        result = bundle._normalize_sample_metadata(
            sample_ids=["S1"], autoload=True
        )

        res.load.assert_called_once_with()
        assert list(result.index) == ["S1"]

    def test_empty_metadata_tables_are_dropped_then_reported(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        empty = pd.DataFrame(
            {"external_id": pd.array([], dtype="string"), "score": []}
        )
        bundle = R3ResourceBundle(resources=[_metadata_resource(empty)])
        with caplog.at_level(logging.WARNING):
            with pytest.raises(
                ValueError,
                match="empty or could not be retrieved",
            ):
                bundle._normalize_sample_metadata(sample_ids=["S1"])
        assert "Dropping empty metadata table" in caplog.text

    def test_duplicate_metadata_keys_are_rejected(self) -> None:
        frame = pd.DataFrame(
            {
                "external_id": pd.array(["S1", "S1"], dtype="string"),
                "rail_id": pd.array(["1", "1"], dtype="string"),
                "study": pd.array(["SRP001", "SRP001"], dtype="string"),
                "score": [1.0, 2.0],
            }
        )
        bundle = R3ResourceBundle(resources=[_metadata_resource(frame)])
        with pytest.raises(ValueError, match="Duplicate sample metadata keys"):
            bundle._normalize_sample_metadata(sample_ids=["S1"])

    def test_repeated_columns_across_tables_are_rejected(self) -> None:
        """Two tables sharing an origin namespace collide on their columns."""
        frame = pd.DataFrame(
            {
                "external_id": pd.array(["S1"], dtype="string"),
                "rail_id": pd.array(["1"], dtype="string"),
                "study": pd.array(["SRP001"], dtype="string"),
                "score": [1.0],
            }
        )
        resources = [
            _metadata_resource(frame.copy(), table_name="recount_qc"),
            _metadata_resource(frame.copy(), table_name="recount_qc"),
        ]
        bundle = R3ResourceBundle(resources=resources)
        with pytest.raises(ValueError, match="Repeated metadata table columns"):
            bundle._normalize_sample_metadata(sample_ids=["S1"])

    def test_identifier_reused_across_rows_is_ambiguous(self) -> None:
        frame = pd.DataFrame(
            {
                "external_id": pd.array(["A", "B"], dtype="string"),
                "rail_id": pd.array([pd.NA, "A"], dtype="string"),
                "study": pd.array(["SRP001", "SRP001"], dtype="string"),
                "score": [1.0, 2.0],
            }
        )
        bundle = R3ResourceBundle(resources=[_metadata_resource(frame)])
        with pytest.raises(ValueError, match="Ambiguous sample identifier 'A'"):
            bundle._normalize_sample_metadata(sample_ids=["A", "B"])

    def test_inner_join_rejects_metadata_rows_without_counts(self) -> None:
        frame = pd.DataFrame(
            {
                "external_id": pd.array(["S1", "S2"], dtype="string"),
                "rail_id": pd.array(["1", "2"], dtype="string"),
                "study": pd.array(["SRP001", "SRP001"], dtype="string"),
                "score": [1.0, 2.0],
            }
        )
        bundle = R3ResourceBundle(resources=[_metadata_resource(frame)])
        with pytest.raises(
            ValueError, match="samples missing from the counts matrix"
        ):
            bundle._normalize_sample_metadata(
                sample_ids=["S1"], metadata_join="inner"
            )

    def test_inner_join_rejects_matched_rows_without_external_id(
        self,
    ) -> None:
        frame = pd.DataFrame(
            {
                "rail_id": pd.array(["1"], dtype="string"),
                "study": pd.array(["SRP001"], dtype="string"),
                "score": [1.0],
            }
        )
        bundle = R3ResourceBundle(resources=[_metadata_resource(frame)])
        with pytest.raises(ValueError, match="missing external_id values"):
            bundle._normalize_sample_metadata(
                sample_ids=["1"], metadata_join="inner"
            )


class TestPrepareExperimentValidation:
    def test_rejects_an_unknown_metadata_join(self) -> None:
        bundle = R3ResourceBundle()
        with pytest.raises(ValueError, match="metadata_join must be"):
            bundle._prepare_experiment(
                genomic_unit="gene",
                annotation_extension=None,
                join_policy="inner",
                metadata_join="left",
                autoload=False,
            )

    def test_requires_count_resources(self) -> None:
        bundle = R3ResourceBundle(resources=[_mock_resource("metadata_files")])
        with pytest.raises(
            ValueError, match="No count-file resources available"
        ):
            bundle._prepare_experiment(
                genomic_unit="gene",
                annotation_extension=None,
                join_policy="inner",
                metadata_join="inner",
                autoload=False,
            )

    def test_unloaded_metadata_is_reported_before_loading_counts(self) -> None:
        project = {
            "organism": "human",
            "data_source": "sra",
            "project": "SRP001",
        }
        counts = _mock_resource(
            "count_files_gene_or_exon",
            loaded_data=_gene_df(),
            genomic_unit="gene",
            annotation_extension="G026",
            **project,
        )
        metadata = _mock_resource(
            "metadata_files", table_name="recount_qc", **project
        )
        bundle = R3ResourceBundle(resources=[counts, metadata])
        with pytest.raises(ValueError, match="Metadata resource is not loaded"):
            bundle._prepare_experiment(
                genomic_unit="gene",
                annotation_extension=None,
                join_policy="inner",
                metadata_join="inner",
                autoload=False,
            )

    def test_samples_shared_between_projects_are_rejected(self) -> None:
        resources = [
            _mock_resource(
                "count_files_gene_or_exon",
                url=f"http://example.com/{project}.gz",
                loaded_data=_gene_df(),
                genomic_unit="gene",
                annotation_extension="G026",
                organism="human",
                data_source="sra",
                project=project,
            )
            for project in ("SRP001", "SRP002")
        ]
        bundle = R3ResourceBundle(resources=resources)
        with pytest.raises(
            ValueError,
            match="Duplicate sample identifiers across selected count",
        ):
            bundle._prepare_experiment(
                genomic_unit="gene",
                annotation_extension=None,
                join_policy="inner",
                metadata_join="inner",
                autoload=False,
            )


class TestStackCountMatricesFeatureCompat:
    def test_a_single_feature_key_satisfies_feature_compat(self) -> None:
        resources = [
            _mock_resource(
                "count_files_gene_or_exon",
                loaded_data=_gene_df(samples=samples),
                genomic_unit="gene",
            )
            for samples in (["SRR001", "SRR002"], ["SRR003", "SRR004"])
        ]
        bundle = R3ResourceBundle(resources=resources)
        out = bundle.stack_count_matrices(
            compat="feature", axis=1, autoload=False
        )
        assert list(out.columns) == [
            "SRR001",
            "SRR002",
            "SRR003",
            "SRR004",
        ]


@pytest.mark.requires_biocpy
class TestToRangedSummarizedExperimentResourceUrls:
    def test_annotation_url_is_recorded_only_once(self, tmp_path: Path) -> None:
        """An annotation already listed in the metadata is not re-appended."""
        gz_path = tmp_path / "genes.gtf.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as handle:
            handle.write(
                'chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG001.1"\n'
                'chr1\tref\tgene\t200\t300\t.\t-\t.\tgene_id "ENSG002.1"\n'
            )
        shared_url = "http://example.com/human.gene_sums.G026.gz"
        counts = _gene_df(features=["ENSG001.1", "ENSG002.1"])

        res_count = _mock_resource(
            "count_files_gene_or_exon",
            url=shared_url,
            loaded_data=counts,
            genomic_unit="gene",
            annotation_extension="G026",
        )
        res_ann = _annotation_resource(
            cached_path=gz_path,
            url_path="human/annotations/gene_sums/human.gene_sums.G026.gtf.gz",
            url=shared_url,
            genomic_unit="gene",
        )
        bundle = R3ResourceBundle(resources=[res_count, res_ann])

        rse = bundle.to_ranged_summarized_experiment(
            genomic_unit="gene", autoload=False
        )
        assert rse.get_metadata()["resource_urls"] == [shared_url]


class TestNumericAssaySparseStorage:
    def test_all_sparse_counts_stay_a_sparse_matrix(self) -> None:
        counts = pd.DataFrame(
            {
                "S1": pd.arrays.SparseArray([0.0, 2.0], fill_value=0.0),
                "S2": pd.arrays.SparseArray([3.0, 0.0], fill_value=0.0),
            },
            index=["G1", "G2"],
        )
        matrix = bmod._numeric_assay(counts)
        assert scipy.sparse.issparse(matrix)
        np.testing.assert_array_equal(
            matrix.toarray(), [[0.0, 3.0], [2.0, 0.0]]
        )


class TestValidateCountFrameColumnKinds:
    def test_sparse_column_with_zero_fill_is_accepted(self) -> None:
        frame = pd.DataFrame(
            {"S1": pd.arrays.SparseArray([0.0, 3.0], fill_value=0.0)}
        )
        assert bmod._validate_count_frame(frame) is None

    def test_non_numeric_column_is_rejected(self) -> None:
        frame = pd.DataFrame({"S1": ["abc", "def"]})
        with pytest.raises(ValueError, match="non-numeric values"):
            bmod._validate_count_frame(frame)


class TestExperimentFrameMismatch:
    def test_length_mismatch_names_both_dimensions(self) -> None:
        frame = pd.DataFrame({"a": [1]})
        with pytest.raises(
            ValueError, match=r"column_data length 1 != assay dimension 2"
        ):
            bmod._experiment_frame(frame, ["s1", "s2"], "column_data")

    def test_unrelated_identifiers_are_rejected(self) -> None:
        frame = pd.DataFrame({"a": [1, 2]}, index=["x", "y"])
        with pytest.raises(
            ValueError, match="row_data identifiers do not match"
        ):
            bmod._experiment_frame(frame, ["s1", "s2"], "row_data")


class TestMergeCountFramesAcrossFeatureSpaces:
    @staticmethod
    def _frames() -> list[pd.DataFrame]:
        return [
            _gene_df(features=["G1", "G2"], samples=["S1"]),
            _gene_df(features=["G2", "G3"], samples=["S2"]),
        ]

    def test_inner_join_keeps_only_shared_features(self) -> None:
        out = bmod._merge_count_frames(self._frames(), "inner")
        assert list(out.index) == ["G2"]
        assert list(out.columns) == ["S1", "S2"]

    def test_outer_join_unions_features_and_fills_new_rows(self) -> None:
        out = bmod._merge_count_frames(self._frames(), "outer")
        assert list(out.index) == ["G1", "G2", "G3"]
        # fill_value applies only to rows a frame never had.
        assert out.loc["G1", "S2"] == 0
        assert out.loc["G3", "S1"] == 0


class TestMakeUniqueNamesSuffixCollision:
    def test_generated_suffix_skips_a_real_feature_name(self) -> None:
        assert _make_unique_names(["a", "a", "a__dup2"]) == [
            "a",
            "a__dup3",
            "a__dup2",
        ]


class TestNormalizeSampleMetadataIdentityConflicts:
    def test_one_rail_id_mapping_to_two_external_ids_is_rejected(self) -> None:
        frame = pd.DataFrame(
            {
                "rail_id": pd.array(["1", "1"], dtype="string"),
                "external_id": pd.array(["A", "B"], dtype="string"),
                "study": pd.array(["SRP001", "SRP001"], dtype="string"),
                "score": [1.0, 2.0],
            }
        )
        bundle = R3ResourceBundle(resources=[_metadata_resource(frame)])
        with pytest.raises(
            ValueError, match="Conflicting rail_id/external_id mappings"
        ):
            bundle._normalize_sample_metadata(sample_ids=["A", "B"])


class TestPrepareExperimentAssembly:
    _PROJECT = {
        "organism": "human",
        "data_source": "sra",
        "project": "SRP001",
    }

    def test_mixed_annotation_builds_are_rejected(self) -> None:
        resources = [
            _mock_resource(
                "count_files_gene_or_exon",
                url=f"http://example.com/{ext}.gz",
                loaded_data=_gene_df(),
                genomic_unit="gene",
                annotation_extension=ext,
                **self._PROJECT,
            )
            for ext in ("G026", "G029")
        ]
        bundle = R3ResourceBundle(resources=resources)
        with pytest.raises(
            CompatibilityError, match="Ambiguous annotation, organism"
        ):
            bundle._prepare_experiment(
                genomic_unit="gene",
                annotation_extension=None,
                join_policy="inner",
                metadata_join="inner",
                autoload=False,
            )

    def test_already_loaded_metadata_is_merged_without_reloading(self) -> None:
        counts = _mock_resource(
            "count_files_gene_or_exon",
            url="http://example.com/gene.gz",
            loaded_data=_gene_df(),
            genomic_unit="gene",
            annotation_extension="G026",
            **self._PROJECT,
        )
        meta_frame = pd.DataFrame(
            {
                "external_id": pd.array(["SRR001", "SRR002"], dtype="string"),
                "rail_id": pd.array(["1", "2"], dtype="string"),
                "study": pd.array(["SRP001", "SRP001"], dtype="string"),
                "score": [1.0, 2.0],
            }
        )
        metadata_res = _metadata_resource(meta_frame, **self._PROJECT)
        bundle = R3ResourceBundle(resources=[counts, metadata_res])

        _, _, col, ranges, metadata = bundle._prepare_experiment(
            genomic_unit="gene",
            annotation_extension=None,
            join_policy="inner",
            metadata_join="inner",
            autoload=False,
        )

        metadata_res.load.assert_not_called()
        assert list(col.index) == ["SRR001", "SRR002"]
        assert list(col["recount_qc__score"]) == [1.0, 2.0]
        assert col["BigWigURL"].str.contains("SRP001").all()
        assert ranges is None
        assert metadata["annotation"] == "gencode_v26"
        assert metadata["annotation_extension"] == "G026"

    def test_junction_projects_are_aligned_on_rr_coordinates(self) -> None:
        mm1, rr1 = _junction_resources(
            project="SRP001", samples=["S1", "S2"], rr_frame=_rr_frame()
        )
        mm2, rr2 = _junction_resources(
            project="SRP002", samples=["S3", "S4"], rr_frame=_rr_frame()
        )
        bundle = R3ResourceBundle(resources=[mm1, rr1, mm2, rr2])

        counts, _, _, ranges, metadata = bundle._prepare_experiment(
            genomic_unit="junction",
            annotation_extension=None,
            join_policy="inner",
            metadata_join="inner",
            autoload=False,
        )

        coordinates = ["chr1:1-100:+", "chr1:200-300:-"]
        assert list(counts.index) == coordinates
        assert list(counts.columns) == ["S1", "S2", "S3", "S4"]
        # The per-project range frames are deduplicated back to one row each.
        assert list(ranges.index) == coordinates
        assert metadata["jxn_format"] == "ALL"


@pytest.mark.requires_biocpy
class TestToRangedSummarizedExperimentRangeReuse:
    def test_a_repeated_build_reuses_the_cached_alignment(
        self, tmp_path: Path
    ) -> None:
        gz_path = tmp_path / "genes.gtf.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as handle:
            handle.write(
                'chr1\tref\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG0001"\n'
                'chr1\tref\tgene\t200\t300\t.\t-\t.\tgene_id "ENSG0002"\n'
            )
        res_count = _mock_resource(
            "count_files_gene_or_exon",
            url="http://example.com/gene.gz",
            loaded_data=_gene_df(),
            genomic_unit="gene",
            annotation_extension="G026",
        )
        res_ann = _annotation_resource(
            cached_path=gz_path,
            url_path="human/annotations/gene_sums/human.gene_sums.G026.gtf.gz",
            url="http://example.com/human.gene_sums.G026.gtf.gz",
            genomic_unit="gene",
        )
        bundle = R3ResourceBundle(resources=[res_count, res_ann])

        first = bundle.to_ranged_summarized_experiment(
            genomic_unit="gene", autoload=False
        )
        with patch.object(bmod, "_read_gtf_dataframe") as reader:
            second = bundle.to_ranged_summarized_experiment(
                genomic_unit="gene", autoload=False
            )

        reader.assert_not_called()
        assert list(second.get_row_names()) == list(first.get_row_names())
        assert list(second.get_row_ranges().get_start()) == [1, 200]

    def test_multi_project_junctions_keep_the_prepared_ranges(self) -> None:
        mm1, rr1 = _junction_resources(
            project="SRP001", samples=["S1", "S2"], rr_frame=_rr_frame()
        )
        mm2, rr2 = _junction_resources(
            project="SRP002", samples=["S3", "S4"], rr_frame=_rr_frame()
        )
        bundle = R3ResourceBundle(resources=[mm1, rr1, mm2, rr2])

        rse = bundle.to_ranged_summarized_experiment(
            genomic_unit="junction", autoload=False
        )

        assert list(rse.get_row_names()) == [
            "chr1:1-100:+",
            "chr1:200-300:-",
        ]
        assert list(rse.get_row_ranges().get_start()) == [1, 200]
