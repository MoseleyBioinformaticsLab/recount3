# Copyright (c) 2026, Alexander Alsalihi, Robert M. Flight, Hunter N.B. Moseley.
# All rights reserved.
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


from __future__ import annotations

import dataclasses
import gzip
import io
import logging
import pathlib
import threading
from typing import Any
from unittest.mock import MagicMock
from unittest.mock import patch
import urllib.error
import urllib.parse
import zipfile

import numpy as np
import pandas as pd
import pytest
import scipy.sparse

from recount3._descriptions import R3ResourceDescription
import recount3._utils as _utils_module
import recount3.bundle as bmod
from recount3.bundle import _align_ranges_to_features
from recount3.bundle import _choose_alignment_key
from recount3.bundle import _classify_ranges_failure
from recount3.bundle import _coerce_gtf_bp_length
from recount3.bundle import _coerce_gtf_phase_column
from recount3.bundle import _collapse_rows_by_key
from recount3.bundle import _construct_ranged_summarized_experiment
from recount3.bundle import _construct_summarized_experiment
from recount3.bundle import _count_compat_keys
from recount3.bundle import _dedupe_ranges_on_feature_id
from recount3.bundle import _ensure_unique_columns
from recount3.bundle import _load_optional_metadata
from recount3.bundle import _make_unique_names
from recount3.bundle import _maybe_relabel_counts_columns_to_external_id
from recount3.bundle import _metadata_origin
from recount3.bundle import _namespace_metadata_columns
from recount3.bundle import _outer_merge_metadata_frames
from recount3.bundle import _parse_gtf_attributes
from recount3.bundle import _peek_gtf_feature_counts
from recount3.bundle import _ranges_from_gtf
from recount3.bundle import _read_gtf_dataframe
from recount3.bundle import _read_rr_table
from recount3.bundle import _select_gtf_resource_for_unit
from recount3.bundle import _standardize_metadata_frame
from recount3.bundle import _strip_ensembl_version
from recount3.bundle import _to_genomic_ranges
from recount3.bundle import R3ResourceBundle
from recount3.config import Config
from recount3.errors import CompatibilityError
from recount3.errors import DownloadError
from recount3.errors import LoadError
from recount3.errors import MissingRangesError
from recount3.errors import RangesCoverageError
from recount3.errors import RangesError
from recount3.resource import R3Resource
import recount3.search as search_module

_TESTS_DIR = pathlib.Path(__file__).parent
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
def local_config(tmp_path: pathlib.Path) -> Config:
    """Return a Config using tmp_path and the local test-data mirror."""
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
        self, tmp_path: pathlib.Path
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

    def test_fallback_reads_gz_file(self, tmp_path: pathlib.Path) -> None:
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

    def test_reads_plain_gtf(self, tmp_path: pathlib.Path) -> None:
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
        self, tmp_path: pathlib.Path, suffix: str
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

    def test_max_lines_limits_scan(self, tmp_path: pathlib.Path) -> None:
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

    def test_comment_lines_skipped(self, tmp_path: pathlib.Path) -> None:
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

    def test_returns_candidate_containing_feature(
        self, tmp_path: pathlib.Path
    ) -> None:
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

    def test_filters_by_annotation_extension(
        self, tmp_path: pathlib.Path
    ) -> None:
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

    def test_returns_best_when_no_feature_found(
        self, tmp_path: pathlib.Path
    ) -> None:
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
    def _synthetic_gtf_gz(self, tmp_path: pathlib.Path) -> pathlib.Path:
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

    def _bundle_with_gene_gtf(self, gtf_path: pathlib.Path) -> R3ResourceBundle:
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
    def test_gene_with_gtf_annotation(
        self, _synthetic_gtf_gz: pathlib.Path
    ) -> None:
        b = self._bundle_with_gene_gtf(_synthetic_gtf_gz)
        rse = b.to_ranged_summarized_experiment(
            genomic_unit="gene", autoload=False
        )
        assert rse is not None

    @pytest.mark.requires_biocpy
    def test_gene_with_unstranded_gtf_annotation(
        self, _synthetic_gtf_gz: pathlib.Path
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
    def test_junction_with_rr_coordinates(self, tmp_path: pathlib.Path) -> None:
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
    def test_rr_duplicate_row_names_made_unique(
        self, tmp_path: pathlib.Path
    ) -> None:
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
        self, tmp_path: pathlib.Path
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

    def test_pathlike_dest_forwarded(self) -> None:
        res = _mock_resource()
        b = R3ResourceBundle(resources=[res])
        dest = pathlib.Path("/some/dir")
        b.download(dest=dest)
        res.download.assert_called_once_with(
            path=dest, cache_mode="enable", overwrite=False
        )

    def test_pathlike_dest_forwarded_in_parallel(self) -> None:
        resources = [_mock_resource() for _ in range(3)]
        b = R3ResourceBundle(resources=resources)
        dest = pathlib.Path("/some/dir")
        b.download(dest=dest, max_workers=4)
        for res in resources:
            res.download.assert_called_once_with(
                path=dest, cache_mode="enable", overwrite=False
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
        self, local_config: Config, tmp_path: pathlib.Path
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

        def transfer(url: str, path: pathlib.Path, **_: Any) -> None:
            barrier.wait()
            path.write_bytes(url.encode())

        target = tmp_path / "bundle.zip"
        with patch.object(_utils_module, "download_to_file", transfer):
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

    def test_line_with_two_fields_skipped(self, tmp_path: pathlib.Path) -> None:
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

    def test_exon_unit_scoring(self, tmp_path: pathlib.Path) -> None:
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
        self, tmp_path: pathlib.Path
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

    def test_genomic_unit_matches_adds_score(
        self, tmp_path: pathlib.Path
    ) -> None:
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
    def test_autoload_true_reuses_cached_gtf(
        self, tmp_path: pathlib.Path
    ) -> None:
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
        self, tmp_path: pathlib.Path
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

        def fake_ensure_cached(*, download: bool) -> pathlib.Path:
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
    def test_enrich_cols_joined_to_row_data(
        self, tmp_path: pathlib.Path
    ) -> None:
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

    def test_fallback_to_file_when_load_raises(
        self, tmp_path: pathlib.Path
    ) -> None:
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

    def test_peek_exception_is_logged_and_skipped(
        self, tmp_path: pathlib.Path
    ) -> None:
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
        with patch("recount3.config.default_config", return_value=mock_cfg):
            result = b._add_bigwig_urls(col_df)
        assert "BigWigURL" in result.columns
        url = result["BigWigURL"].iloc[0]
        assert url is not None
        assert url.startswith("https://duffel.example.com/recount3/")


def _annotation_resource(
    *,
    cached_path: pathlib.Path,
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


def _write_gtf(path: pathlib.Path, *, feature: str = "gene") -> pathlib.Path:
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

    def test_two_exact_unit_matches_are_ambiguous(
        self, tmp_path: pathlib.Path
    ) -> None:
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
        self, tmp_path: pathlib.Path, caplog: pytest.LogCaptureFixture
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
        self, tmp_path: pathlib.Path
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
        self, tmp_path: pathlib.Path, caplog: pytest.LogCaptureFixture
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
        self, tmp_path: pathlib.Path, caplog: pytest.LogCaptureFixture
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
        self, tmp_path: pathlib.Path
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

    def test_annotation_url_is_recorded_only_once(
        self, tmp_path: pathlib.Path
    ) -> None:
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
        self, tmp_path: pathlib.Path
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


# ---------------------------------------------------------------------------
# Characterization of the count merge, project grouping, BigWig URL, and
# sample-lookup paths. The reference implementations below are the original
# code, kept verbatim as oracles, so faster rewrites must reproduce them
# exactly (values, dtypes, sparse storage, labels, names and attrs).
# ---------------------------------------------------------------------------


def _reference_merge_count_frames(
    frames: list[pd.DataFrame], join_policy: str
) -> pd.DataFrame:
    """Return the original ``_merge_count_frames`` result."""
    if join_policy not in ("inner", "outer"):
        raise ValueError("join_policy must be 'inner' or 'outer'.")
    if len(frames) == 1:
        return frames[0].copy(deep=False)
    if all(frame.index.equals(frames[0].index) for frame in frames[1:]):
        return pd.concat(frames, axis=1)
    if any(not frame.index.is_unique for frame in frames):
        raise CompatibilityError(
            "Cannot align duplicate feature IDs across projects."
        )
    index = frames[0].index
    for frame in frames[1:]:
        index = (
            index.intersection(frame.index, sort=False)
            if join_policy == "inner"
            else index.union(frame.index, sort=False)
        )
    return pd.concat(
        [frame.reindex(index, fill_value=0) for frame in frames], axis=1
    )


def _assert_identical_frames(
    result: pd.DataFrame, expected: pd.DataFrame
) -> None:
    """Assert equality down to sparse storage, dtypes, names and attrs."""
    pd.testing.assert_frame_equal(result, expected, check_exact=True)
    pd.testing.assert_index_equal(result.index, expected.index, exact=True)
    pd.testing.assert_index_equal(result.columns, expected.columns, exact=True)
    assert result.index.name == expected.index.name
    assert result.columns.name == expected.columns.name
    assert result.attrs == expected.attrs
    assert [str(d) for d in result.dtypes] == [str(d) for d in expected.dtypes]
    for position in range(expected.shape[1]):
        got = result.iloc[:, position].array
        want = expected.iloc[:, position].array
        assert type(got) is type(want)
        if isinstance(want, pd.arrays.SparseArray):
            assert got.kind == want.kind
            assert type(got.fill_value) is type(want.fill_value)
            assert got.sp_values.dtype == want.sp_values.dtype
            np.testing.assert_array_equal(
                got.sp_index.to_int_index().indices,
                want.sp_index.to_int_index().indices,
            )
            np.testing.assert_array_equal(got.sp_values, want.sp_values)


def _loader_frame(
    dense: np.ndarray,
    labels: list[str],
    columns: list[str],
    *,
    stored: np.ndarray | None = None,
    dtype: Any = np.int64,
    name: str | None = None,
) -> pd.DataFrame:
    """Return a sparse frame built the way the junction loader builds one.

    ``stored`` marks entries to store even when they are zero, so the frame
    can carry explicit zeros the way a MatrixMarket file can.
    """
    mask = dense != 0 if stored is None else stored
    rows, cols = np.nonzero(mask)
    matrix = scipy.sparse.coo_matrix(
        (dense[rows, cols].astype(dtype), (rows, cols)), shape=dense.shape
    ).tocsc()
    frame = pd.DataFrame.sparse.from_spmatrix(matrix)
    frame.index = pd.Index(labels, name=name)
    frame.columns = columns
    return frame


_MERGE_FLAVORS = (
    "loader_int64",
    "loader_float64",
    "loader_float64_nan",
    "loader_mixed_subtypes",
    "dense_int64",
    "dense_float64",
    "sparse_int64",
    "sparse_int32",
    "sparse_float32",
    "sparse_uint64",
    "sparse_bool",
    "sparse_block",
    "sparse_float64_numpy_fill",
    "sparse_int64_numpy_fill",
    "mixed_dense_sparse",
)


def _merge_frame(
    rng: np.random.Generator,
    flavor: str,
    labels: list[str],
    columns: list[str],
    name: str | None,
) -> pd.DataFrame:
    """Return one random count frame of the requested storage flavor."""
    shape = (len(labels), len(columns))
    dense = rng.integers(0, 4, size=shape) * (rng.random(shape) < 0.5)
    stored = (dense != 0) | (rng.random(shape) < 0.15)
    if flavor.startswith("loader_"):
        values = dense.astype(float)
        if flavor == "loader_float64_nan":
            values[rng.random(shape) < 0.1] = np.nan
            stored |= np.isnan(values)
        if flavor == "loader_mixed_subtypes":
            frame = _loader_frame(
                values, labels, columns, stored=stored, name=name
            )
            for position in range(0, len(columns), 2):
                column = _loader_frame(
                    values[:, [position]],
                    labels,
                    [columns[position]],
                    stored=stored[:, [position]],
                    dtype=np.float64,
                    name=name,
                )
                frame.isetitem(position, column.iloc[:, 0])
            return frame
        dtype = np.int64 if flavor == "loader_int64" else np.float64
        return _loader_frame(
            values, labels, columns, stored=stored, dtype=dtype, name=name
        )
    index = pd.Index(labels, name=name)
    if flavor.startswith("dense_"):
        dtype = np.int64 if flavor == "dense_int64" else np.float64
        return pd.DataFrame(dense.astype(dtype), index=index, columns=columns)
    specs = {
        "sparse_int64": (np.int64, 0, "integer"),
        "sparse_int32": (np.int32, 0, "integer"),
        "sparse_float32": (np.float32, 0, "integer"),
        "sparse_uint64": (np.uint64, 0, "integer"),
        "sparse_bool": (bool, False, "integer"),
        "sparse_block": (np.int64, 0, "block"),
        "sparse_float64_numpy_fill": (np.float64, np.float64(0.0), "integer"),
        "sparse_int64_numpy_fill": (np.int64, np.int64(0), "integer"),
        "mixed_dense_sparse": (np.int64, 0, "integer"),
    }
    dtype, fill_value, kind = specs[flavor]
    data = {}
    for position, column in enumerate(columns):
        values = dense[:, position].astype(dtype)
        if flavor == "mixed_dense_sparse" and position % 2:
            data[column] = values
        else:
            data[column] = pd.arrays.SparseArray(
                values, fill_value=fill_value, kind=kind
            )
    return pd.DataFrame(data, index=index, columns=columns)


def _random_merge_case(seed: int) -> tuple[list[pd.DataFrame], str]:
    """Return random frames with overlapping features and a join policy."""
    rng = np.random.default_rng(seed)
    pool = [f"F{i}" for i in range(int(rng.integers(1, 25)))]
    name = [None, "gene_id"][int(rng.integers(0, 2))]
    attrs = [{}, {"origin": "test"}][int(rng.integers(0, 2))]
    frames = []
    for project in range(int(rng.integers(2, 6))):
        size = int(rng.integers(0, len(pool) + 1))
        labels = [str(x) for x in rng.choice(pool, size=size, replace=False)]
        columns = [f"P{project}S{j}" for j in range(int(rng.integers(0, 4)))]
        flavor = _MERGE_FLAVORS[int(rng.integers(0, len(_MERGE_FLAVORS)))]
        frame = _merge_frame(rng, flavor, labels, columns, name)
        frame.attrs = dict(attrs)
        frames.append(frame)
    if rng.random() < 0.2:
        frames.append(frames[0].iloc[::-1].copy())
    return frames, ["inner", "outer"][int(rng.integers(0, 2))]


class TestMergeCountFramesCharacterization:
    """Pin every observable property of the cross-project count merge."""

    @staticmethod
    def _three_loader_frames(dtype: Any = np.int64) -> list[pd.DataFrame]:
        first = _loader_frame(
            np.array([[1, 0], [0, 2], [3, 0]]),
            ["a", "b", "c"],
            ["S1", "S2"],
            dtype=dtype,
        )
        second = _loader_frame(
            np.array([[4], [5], [0]]), ["d", "b", "a"], ["S3"], dtype=dtype
        )
        third = _loader_frame(
            np.array([[6, 7], [0, 8]]), ["e", "a"], ["S4", "S5"], dtype=dtype
        )
        return [first, second, third]

    def test_outer_join_lists_features_in_first_seen_order(self) -> None:
        out = bmod._merge_count_frames(self._three_loader_frames(), "outer")
        assert list(out.index) == ["a", "b", "c", "d", "e"]
        assert list(out.columns) == ["S1", "S2", "S3", "S4", "S5"]
        expected = np.array(
            [
                [1, 0, 0, 0, 8],
                [0, 2, 5, 0, 0],
                [3, 0, 0, 0, 0],
                [0, 0, 4, 0, 0],
                [0, 0, 0, 6, 7],
            ]
        )
        np.testing.assert_array_equal(out.to_numpy(), expected)
        assert all(str(d) == "Sparse[int64, 0]" for d in out.dtypes)

    def test_inner_join_keeps_the_first_frames_order(self) -> None:
        frames = self._three_loader_frames()
        frames[2] = _loader_frame(
            np.array([[6], [7], [9]]), ["b", "c", "a"], ["S4"]
        )
        frames[1] = _loader_frame(
            np.array([[4], [5], [0]]), ["c", "a", "b"], ["S3"]
        )
        out = bmod._merge_count_frames(frames, "inner")
        assert list(out.index) == ["a", "b", "c"]
        np.testing.assert_array_equal(
            out.to_numpy(),
            np.array([[1, 0, 5, 9], [0, 2, 0, 6], [3, 0, 4, 7]]),
        )

    def test_reindexed_columns_keep_rows_sorted(self) -> None:
        out = bmod._merge_count_frames(self._three_loader_frames(), "outer")
        for position in range(out.shape[1]):
            rows = out.iloc[:, position].array.sp_index.to_int_index().indices
            assert np.all(np.diff(rows) > 0)

    def test_reindexing_drops_stored_zeros(self) -> None:
        dense = np.array([[0, 2], [5, 0]])
        stored = np.ones_like(dense, dtype=bool)
        moved = _loader_frame(dense, ["y", "x"], ["S2", "S3"], stored=stored)
        anchor = _loader_frame(np.array([[1], [1]]), ["x", "z"], ["S1"])
        out = bmod._merge_count_frames([anchor, moved], "outer")
        assert list(out.index) == ["x", "z", "y"]
        s2 = out["S2"].array
        assert list(s2.sp_index.to_int_index().indices) == [0]
        assert list(s2.sp_values) == [5]

    def test_a_frame_already_on_the_joined_index_keeps_stored_zeros(
        self,
    ) -> None:
        stored = np.ones((2, 1), dtype=bool)
        full = _loader_frame(
            np.array([[0], [3]]), ["x", "y"], ["S1"], stored=stored
        )
        subset = _loader_frame(np.array([[1]]), ["y"], ["S2"])
        out = bmod._merge_count_frames([full, subset], "outer")
        assert list(out["S1"].array.sp_values) == [0, 3]
        assert list(out["S2"].array.sp_index.to_int_index().indices) == [1]

    def test_stored_missing_counts_survive_reindexing(self) -> None:
        dense = np.array([[np.nan], [2.0]])
        moved = _loader_frame(dense, ["y", "x"], ["S2"], dtype=np.float64)
        anchor = _loader_frame(
            np.array([[1.0]]), ["x"], ["S1"], dtype=np.float64
        )
        frames = [anchor, moved]
        out = bmod._merge_count_frames(frames, "outer")
        assert list(out.index) == ["x", "y"]
        assert out["S2"].loc["x"] == 2.0
        assert pd.isna(out["S2"].loc["y"])
        _assert_identical_frames(
            out, _reference_merge_count_frames(frames, "outer")
        )

    def test_an_empty_inner_join_has_no_rows_and_keeps_dtypes(self) -> None:
        first = _loader_frame(np.array([[1]]), ["a"], ["S1"])
        second = _loader_frame(
            np.array([[2.0]]), ["b"], ["S2"], dtype=np.float64
        )
        out = bmod._merge_count_frames([first, second], "inner")
        assert out.shape == (0, 2)
        _assert_identical_frames(
            out, _reference_merge_count_frames([first, second], "inner")
        )

    @pytest.mark.parametrize(
        ("names", "expected"),
        [
            (["gene_id", "gene_id", "gene_id"], "gene_id"),
            (["gene_id", None, "gene_id"], None),
            (["gene_id", "gene_id", "exon_id"], None),
        ],
    )
    def test_index_names_reconcile_like_set_operations(
        self, names: list[str | None], expected: str | None
    ) -> None:
        frames = self._three_loader_frames()
        for frame, name in zip(frames, names):
            frame.index.name = name
        for policy in ("inner", "outer"):
            out = bmod._merge_count_frames(frames, policy)
            assert out.index.name == expected

    def test_shared_attrs_are_kept(self) -> None:
        frames = self._three_loader_frames()
        for frame in frames:
            frame.attrs = {"origin": "recount3"}
        out = bmod._merge_count_frames(frames, "outer")
        assert out.attrs == {"origin": "recount3"}

    def test_inputs_are_not_modified(self) -> None:
        frames = self._three_loader_frames()
        before = [frame.copy() for frame in frames]
        bmod._merge_count_frames(frames, "outer")
        for frame, copy in zip(frames, before):
            _assert_identical_frames(frame, copy)

    def test_repeated_sample_labels_across_frames_are_kept(self) -> None:
        first = _loader_frame(np.array([[1], [2]]), ["a", "b"], ["S1"])
        second = _loader_frame(np.array([[3]]), ["b"], ["S1"])
        out = bmod._merge_count_frames([first, second], "outer")
        assert list(out.columns) == ["S1", "S1"]
        np.testing.assert_array_equal(out.to_numpy(), [[1, 0], [2, 3]])

    def test_mixed_index_dtypes_join_like_pandas(self) -> None:
        first = _gene_df(features=["G1", "G2"], samples=["S1"])
        second = _gene_df(features=["G2", "G3"], samples=["S2"])
        second.index = pd.Index(["G2", "G3"], dtype=object)
        for policy in ("inner", "outer"):
            _assert_identical_frames(
                bmod._merge_count_frames([first, second], policy),
                _reference_merge_count_frames([first, second], policy),
            )

    def test_an_empty_frame_joins_like_pandas(self) -> None:
        empty = _loader_frame(np.zeros((0, 1)), [], ["S0"])
        first = _loader_frame(np.array([[1]]), ["a"], ["S1"])
        for frames in ([empty, first], [first, empty]):
            for policy in ("inner", "outer"):
                _assert_identical_frames(
                    bmod._merge_count_frames(frames, policy),
                    _reference_merge_count_frames(frames, policy),
                )

    def test_frames_without_columns_join_like_pandas(self) -> None:
        bare = _loader_frame(np.zeros((2, 0)), ["b", "z"], [])
        first = _loader_frame(np.array([[1]]), ["a"], ["S1"])
        for policy in ("inner", "outer"):
            _assert_identical_frames(
                bmod._merge_count_frames([first, bare], policy),
                _reference_merge_count_frames([first, bare], policy),
            )

    @pytest.mark.parametrize("seed", range(400))
    def test_matches_the_reference_on_random_frames(self, seed: int) -> None:
        frames, policy = _random_merge_case(seed)
        try:
            expected = _reference_merge_count_frames(frames, policy)
        except (CompatibilityError, ValueError) as exc:
            with pytest.raises(type(exc)) as raised:
                bmod._merge_count_frames(frames, policy)
            assert str(raised.value) == str(exc)
            return
        _assert_identical_frames(
            bmod._merge_count_frames(frames, policy), expected
        )


class TestMergeCountFramesContract:
    def test_a_single_frame_is_returned_without_copying_its_data(
        self,
    ) -> None:
        frame = _gene_df()
        out = bmod._merge_count_frames([frame], "outer")
        assert out is not frame
        assert np.shares_memory(out.to_numpy(), frame.to_numpy())

    def test_identical_indexes_may_repeat_feature_ids(self) -> None:
        first = _gene_df(features=["G1", "G1"], samples=["S1"])
        second = _gene_df(features=["G1", "G1"], samples=["S2"])
        out = bmod._merge_count_frames([first, second], "inner")
        assert list(out.index) == ["G1", "G1"]
        assert list(out.columns) == ["S1", "S2"]

    def test_sparse_reindex_at_the_int32_row_boundary(self) -> None:
        last_row = np.iinfo(np.int32).max
        target = pd.RangeIndex(last_row + 1)
        source = pd.DataFrame.sparse.from_spmatrix(
            scipy.sparse.csc_matrix([[7], [13], [11]]),
            index=pd.Index([last_row, -1, 3], dtype=np.int64),
            columns=["S1"],
        )

        result = bmod._reindex_sparse_frame(source, target)

        pd.testing.assert_index_equal(result.index, target)
        column = result["S1"].array
        np.testing.assert_array_equal(column.sp_index.indices, [3, last_row])
        np.testing.assert_array_equal(column.sp_values, [11, 7])

    def test_error_messages(self) -> None:
        with pytest.raises(
            ValueError, match=r"^join_policy must be 'inner' or 'outer'\.$"
        ):
            bmod._merge_count_frames([_gene_df()], "left")
        repeated = _gene_df(features=["G1", "G1"], samples=["S1"])
        with pytest.raises(
            CompatibilityError,
            match=r"^Cannot align duplicate feature IDs across projects\.$",
        ):
            bmod._merge_count_frames([repeated, _gene_df()], "outer")

    def test_object_labels_keep_their_dtype_in_an_outer_join(self) -> None:
        indexes = [
            pd.Index([1, 2], dtype=object),
            pd.Index([3, 1], dtype=object),
        ]
        joined = bmod._join_feature_index(indexes, "outer")
        pd.testing.assert_index_equal(
            joined, indexes[0].union(indexes[1], sort=False), exact=True
        )
        frames = [
            _loader_frame(np.array([[1], [2]]), [1, 2], ["S1"]),
            _loader_frame(np.array([[3], [4]]), [3, 1], ["S2"]),
        ]
        for frame, index in zip(frames, indexes):
            frame.index = index
        _assert_identical_frames(
            bmod._merge_count_frames(frames, "outer"),
            _reference_merge_count_frames(frames, "outer"),
        )

    @pytest.mark.parametrize(
        ("first_labels", "second_labels"),
        [
            ([2**53 + 1], [float(2**53)]),
            ([5, np.nan], [1, None]),
        ],
        ids=["integers-above-float-precision", "distinct-missing-labels"],
    )
    def test_object_labels_keep_distinct_values(
        self, first_labels: list[Any], second_labels: list[Any]
    ) -> None:
        indexes = [
            pd.Index(labels, dtype=object)
            for labels in (first_labels, second_labels)
        ]
        joined = bmod._join_feature_index(indexes, "outer")
        pd.testing.assert_index_equal(
            joined, indexes[0].union(indexes[1], sort=False), exact=True
        )
        assert len(joined) == sum(map(len, indexes))

        frames = [
            _loader_frame(
                np.arange(1, len(index) + 1)[:, None],
                list(index),
                [f"S{position}"],
            ).set_axis(index)
            for position, index in enumerate(indexes)
        ]
        _assert_identical_frames(
            bmod._merge_count_frames(frames, "outer"),
            _reference_merge_count_frames(frames, "outer"),
        )


def _project_resources(
    project: str,
    samples: list[str],
    *,
    tables: tuple[str, ...] = ("recount_qc",),
    features: list[str] | None = None,
) -> list[MagicMock]:
    """Return a loaded gene-count resource and metadata tables for a project."""
    identity = {"organism": "human", "data_source": "sra", "project": project}
    counts = _mock_resource(
        "count_files_gene_or_exon",
        url=f"http://example.com/{project}.gene.gz",
        loaded_data=_gene_df(features=features, samples=samples),
        genomic_unit="gene",
        annotation_extension="G026",
        **identity,
    )
    tables_out = []
    for number, table in enumerate(tables):
        frame = pd.DataFrame(
            {
                "external_id": pd.array(samples, dtype="string"),
                "rail_id": pd.array(
                    [f"{project}-{s}" for s in samples], dtype="string"
                ),
                "study": pd.array([project] * len(samples), dtype="string"),
                f"value{number}": [float(i) for i in range(len(samples))],
            }
        )
        res = _metadata_resource(frame, table_name=table, **identity)
        res.url = f"http://example.com/{project}.{table}.MD.gz"
        tables_out.append(res)
    return [counts, *tables_out]


class TestPrepareExperimentProjectGrouping:
    """Projects are grouped by identity in first-seen order."""

    @staticmethod
    def _interleaved_bundle() -> tuple[R3ResourceBundle, list[str]]:
        a_counts, a_qc, a_sra = _project_resources(
            "SRP00A", ["A1", "A2"], tables=("recount_qc", "sra")
        )
        b_counts, b_qc = _project_resources("SRP00B", ["B1"])
        c_counts, c_qc = _project_resources("SRP00C", ["C1", "C2", "C3"])
        unrelated = _mock_resource(
            "annotations",
            url="http://example.com/annotation.gtf.gz",
            organism="human",
            data_source="sra",
            project="SRP00A",
        )
        resources = [
            b_qc,
            a_sra,
            b_counts,
            unrelated,
            c_qc,
            a_counts,
            c_counts,
            a_qc,
        ]
        expected_urls = [
            b_counts.url,
            b_qc.url,
            a_counts.url,
            a_sra.url,
            a_qc.url,
            c_counts.url,
            c_qc.url,
        ]
        return R3ResourceBundle(resources=resources), expected_urls

    def test_interleaved_projects_keep_first_seen_order(self) -> None:
        bundle, expected_urls = self._interleaved_bundle()
        counts, rows, col, ranges, metadata = bundle._prepare_experiment(
            genomic_unit="gene",
            annotation_extension=None,
            join_policy="inner",
            metadata_join="inner",
            autoload=False,
        )
        samples = ["B1", "A1", "A2", "C1", "C2", "C3"]
        assert list(counts.columns) == samples
        assert list(col.index) == samples
        assert (
            list(col["study"]) == ["SRP00B"] + ["SRP00A"] * 2 + ["SRP00C"] * 3
        )
        assert metadata["resource_urls"] == expected_urls
        assert metadata["project"] == ["SRP00B", "SRP00A", "SRP00C"]
        assert list(rows["feature_id"]) == ["ENSG0001", "ENSG0002"]
        assert ranges is None

    def test_each_project_receives_only_its_own_metadata(self) -> None:
        bundle, _ = self._interleaved_bundle()
        _, _, col, _, metadata = bundle._prepare_experiment(
            genomic_unit="gene",
            annotation_extension=None,
            join_policy="inner",
            metadata_join="inner",
            autoload=False,
        )
        assert list(col["recount_qc__value0"]) == [0.0, 0.0, 1.0, 0, 1, 2]
        sra = col["sra__value1"]
        assert sra.isna().tolist() == [True, False, False, True, True, True]
        assert list(sra.iloc[1:3]) == [0.0, 1.0]
        assert set(metadata["metadata_columns"]) >= {
            "recount_qc__value0",
            "sra__value1",
        }
        assert col["BigWigURL"].iloc[0].endswith("SRP00B_B1.ALL.bw")
        assert col["BigWigURL"].iloc[1].endswith("SRP00A_A1.ALL.bw")

    def test_junction_groups_carry_their_own_coordinate_sidecars(
        self,
    ) -> None:
        mm1, rr1 = _junction_resources(
            project="SRP001", samples=["S1"], rr_frame=_rr_frame()
        )
        other = pd.DataFrame(
            {
                "chromosome": ["chr2", "chr1"],
                "start": [5, 1],
                "end": [50, 100],
                "strand": ["+", "+"],
            }
        )
        mm2, rr2 = _junction_resources(
            project="SRP002", samples=["S2"], rr_frame=other
        )
        bundle = R3ResourceBundle(resources=[rr2, mm1, mm2, rr1])
        counts, _, _, ranges, metadata = bundle._prepare_experiment(
            genomic_unit="junction",
            annotation_extension=None,
            join_policy="outer",
            metadata_join="inner",
            autoload=False,
        )
        assert list(counts.index) == [
            "chr1:1-100:+",
            "chr1:200-300:-",
            "chr2:5-50:+",
        ]
        assert metadata["resource_urls"] == [mm1.url, rr1.url, mm2.url, rr2.url]
        assert list(ranges["seqnames"]) == ["chr1", "chr1", "chr2"]
        np.testing.assert_array_equal(
            counts.to_numpy(), [[1.0, 1.0], [1.0, 0.0], [0.0, 1.0]]
        )


class TestPrepareExperimentResourceSelection:
    _IDENTITY = {"organism": "human", "data_source": "sra", "project": "SRP001"}

    def test_a_project_split_across_count_files_is_built_once(self) -> None:
        first = _mock_resource(
            "count_files_gene_or_exon",
            url="http://example.com/part1.gz",
            loaded_data=_gene_df(samples=["S1", "S2"]),
            genomic_unit="gene",
            annotation_extension="G026",
            **self._IDENTITY,
        )
        second = _mock_resource(
            "count_files_gene_or_exon",
            url="http://example.com/part2.gz",
            loaded_data=_gene_df(samples=["S3"]),
            genomic_unit="gene",
            annotation_extension="G026",
            **self._IDENTITY,
        )
        bundle = R3ResourceBundle(resources=[first, second])
        counts, _, col, _, metadata = bundle._prepare_experiment(
            genomic_unit="gene",
            annotation_extension=None,
            join_policy="inner",
            metadata_join="inner",
            autoload=False,
        )
        assert list(counts.columns) == ["S1", "S2", "S3"]
        assert list(col.index) == ["S1", "S2", "S3"]
        assert metadata["resource_urls"] == [first.url, second.url]

    def test_gene_builds_record_sidecars_but_not_junction_counts(self) -> None:
        counts_res, metadata_res = _project_resources("SRP001", ["S1", "S2"])
        mm, rr = _junction_resources(
            project="SRP001", samples=["S1", "S2"], rr_frame=_rr_frame()
        )
        bundle = R3ResourceBundle(resources=[counts_res, mm, metadata_res, rr])
        _, _, _, _, metadata = bundle._prepare_experiment(
            genomic_unit="gene",
            annotation_extension=None,
            join_policy="inner",
            metadata_join="inner",
            autoload=False,
        )
        assert metadata["resource_urls"] == [
            counts_res.url,
            metadata_res.url,
            rr.url,
        ]


class TestStackCountsForJunctionGrouping:
    def test_interleaved_projects_are_keyed_by_their_own_sidecars(
        self,
    ) -> None:
        frames = {
            "SRP001": _rr_frame(),
            "SRP002": pd.DataFrame(
                {
                    "chromosome": ["chr3", "chr1"],
                    "start": [7, 200],
                    "end": [70, 300],
                    "strand": ["-", "-"],
                }
            ),
            "SRP003": pd.DataFrame(
                {
                    "chromosome": ["chr1"],
                    "start": [1],
                    "end": [100],
                    "strand": ["+"],
                }
            ),
        }
        pairs = {
            project: _junction_resources(
                project=project, samples=[f"{project}_S"], rr_frame=rr
            )
            for project, rr in frames.items()
        }
        (mm1, rr1), (mm2, rr2), (mm3, rr3) = pairs.values()
        bundle = R3ResourceBundle(resources=[mm2, rr3, mm1, rr2, mm3, rr1])
        out = bundle._stack_counts_for(
            genomic_unit="junction", join_policy="outer", autoload=False
        )
        assert list(out.index) == [
            "chr3:7-70:-",
            "chr1:200-300:-",
            "chr1:1-100:+",
        ]
        assert list(out.columns) == ["SRP002_S", "SRP001_S", "SRP003_S"]
        np.testing.assert_array_equal(
            out.to_numpy(),
            [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 1.0]],
        )


class TestStackCountsForGeneProjects:
    def test_gene_tables_with_different_features_are_merged(self) -> None:
        identity = {"organism": "human", "annotation_extension": "G026"}
        first = _mock_resource(
            "count_files_gene_or_exon",
            url="http://example.com/a.gz",
            loaded_data=_gene_df(features=["G1", "G2"], samples=["S1"]),
            genomic_unit="gene",
            project="SRP001",
            **identity,
        )
        second = _mock_resource(
            "count_files_gene_or_exon",
            url="http://example.com/b.gz",
            loaded_data=_gene_df(features=["G3", "G2"], samples=["S2"]),
            genomic_unit="gene",
            project="SRP002",
            **identity,
        )
        out = R3ResourceBundle(resources=[first, second])._stack_counts_for(
            genomic_unit="gene", join_policy="outer", autoload=False
        )
        assert list(out.index) == ["G1", "G2", "G3"]
        np.testing.assert_array_equal(
            out.to_numpy(), [[0.0, 0.0], [1.0, 1.0], [0.0, 0.0]]
        )


class TestStackCountsForTables:
    """Validation and labels of the tables ``_stack_counts_for`` combines."""

    @staticmethod
    def _bundle(frame: Any, *, url: str = "http://example.com/g.gz") -> Any:
        res = _mock_resource(
            "count_files_gene_or_exon",
            url=url,
            loaded_data=frame,
            genomic_unit="gene",
            annotation_extension="G026",
            organism="human",
        )
        return R3ResourceBundle(resources=[res])

    @pytest.mark.parametrize("bad", [-1.0, np.nan, np.inf])
    def test_invalid_counts_are_rejected(self, bad: float) -> None:
        frame = _gene_df()
        frame.iloc[0, 0] = bad
        with pytest.raises(
            ValueError,
            match=r"^Counts must be finite, non-missing and nonnegative\.$",
        ):
            self._bundle(frame)._stack_counts_for(
                genomic_unit="gene", autoload=False
            )

    def test_labels_become_text_without_touching_the_loaded_table(
        self,
    ) -> None:
        frame = pd.DataFrame(
            np.arange(4.0).reshape(2, 2), index=[7, 8], columns=[1, 2]
        )
        out = self._bundle(frame)._stack_counts_for(
            genomic_unit="gene", autoload=False
        )
        assert list(out.index) == ["7", "8"]
        assert list(out.columns) == ["1", "2"]
        assert list(frame.index) == [7, 8]
        assert list(frame.columns) == [1, 2]
        assert np.shares_memory(out.to_numpy(), frame.to_numpy())

    def test_missing_feature_identifiers_are_rejected(self) -> None:
        frame = _gene_df()
        frame.index = pd.Index(["G1", np.nan], dtype=object)
        with pytest.raises(ValueError) as raised:
            self._bundle(frame)._stack_counts_for(
                genomic_unit="gene", autoload=False
            )
        assert str(raised.value) == (
            "Missing feature or sample identifiers in http://example.com/g.gz"
        )

    def test_error_messages_name_the_resource(self) -> None:
        url = "http://example.com/g.gz"
        repeated = _gene_df(samples=["S1", "S1"])
        with pytest.raises(ValueError) as raised:
            self._bundle(repeated)._stack_counts_for(
                genomic_unit="gene", autoload=False
            )
        assert str(raised.value) == f"Duplicate sample identifiers in {url}"
        with pytest.raises(TypeError) as raised:
            self._bundle(np.zeros((2, 2)))._stack_counts_for(
                genomic_unit="gene", autoload=False
            )
        assert str(raised.value) == f"Loaded counts are not a DataFrame: {url}"
        unloaded = _mock_resource(
            "count_files_gene_or_exon",
            url=url,
            genomic_unit="gene",
            annotation_extension="G026",
            organism="human",
        )
        bundle = R3ResourceBundle(resources=[unloaded])
        with pytest.raises(ValueError) as raised:
            bundle._stack_counts_for(genomic_unit="gene", autoload=False)
        assert str(raised.value) == (
            f"Count resource is not loaded: {url}. "
            "Load it first or use autoload=True."
        )
        unloaded.load.side_effect = OSError("disk")
        with pytest.raises(ValueError) as raised:
            bundle._stack_counts_for(genomic_unit="gene", autoload=True)
        assert str(raised.value) == (
            f"Failed to load requested count matrix {url}: disk"
        )


class TestRowValues:
    """Per-sample loops read columns as ``iterrows()`` presented them."""

    def test_row_values_of_a_single_row_use_the_shared_dtype(self) -> None:
        frame = pd.DataFrame({"external_id": [12], "auc": [1.5]})
        assert bmod._row_values(frame, "external_id").tolist() == [12.0]

    @pytest.mark.parametrize(
        "frame",
        [
            pd.DataFrame({"x": pd.array([1, None], dtype="Int64")}),
            pd.DataFrame({"x": pd.array([1, None], dtype="Int64"), "y": 1.5}),
            pd.DataFrame({"x": pd.array(["a", None], dtype="string")}),
            pd.DataFrame({"x": [True, False], "y": [1, 2]}),
            pd.DataFrame({"x": pd.Categorical([12, None]), "y": ["a", "b"]}),
        ],
        ids=["masked", "masked_and_float", "text", "bool_and_int", "category"],
    )
    def test_row_values_match_iterrows(self, frame: pd.DataFrame) -> None:
        expected = [row["x"] for _, row in frame.iterrows()]
        values = list(bmod._row_values(frame, "x"))
        assert [type(v) for v in values] == [type(v) for v in expected]
        assert [str(v) for v in values] == [str(v) for v in expected]


def _reference_add_bigwig_urls(
    bundle: R3ResourceBundle, col_df: pd.DataFrame
) -> pd.DataFrame:
    """Return the original row-by-row ``_add_bigwig_urls`` result."""
    out = col_df.copy()
    count_types = {"count_files_gene_or_exon", "count_files_junctions"}
    candidates = [
        r
        for r in bundle.resources
        if getattr(r.description, "resource_type", None) in count_types
    ]
    urls = []
    for _, row in out.iterrows():
        study = row.get("study")
        matching = [
            r
            for r in candidates
            if pd.notna(study)
            and bmod._description_value(r, "project") == str(study)
        ]
        if (
            not matching
            and len({bmod._project_key(r) for r in candidates}) == 1
        ):
            matching = candidates[:1]
        res = matching[0] if matching else None
        external = row.get("external_id")
        if res is None or pd.isna(external) or not str(external).strip():
            urls.append(None)
            continue
        org, source, project = bmod._project_key(res)
        source_fields = [
            c for c in out.columns if c.lower().endswith("file_source")
        ]
        if source_fields and isinstance(row[source_fields[0]], str):
            source = row[source_fields[0]].rstrip("/").rsplit("/", 1)[-1]
        if not all((org, source, project)):
            urls.append(None)
            continue
        urls.append(
            bmod.resource.build_url(
                "bigwig_files",
                organism=org,
                data_source=source,
                project=project,
                sample=str(external),
                config=res.config,
            )
        )
    out["BigWigURL"] = urls
    return out


def _mirror_config(base_url: str) -> Config:
    return Config(
        base_url=base_url,
        timeout=5,
        insecure_ssl=False,
        max_retries=1,
        user_agent="test",
        cache_dir="unused-cache",
        cache_disabled=True,
    )


def _count_resource(
    project: str,
    *,
    organism: str = "human",
    data_source: str = "sra",
    genomic_unit: str = "gene",
    base_url: str = "https://mirror-a.example/recount3/",
) -> R3Resource:
    return R3Resource(
        R3ResourceDescription(
            resource_type="count_files_gene_or_exon",
            organism=organism,
            data_source=data_source,
            genomic_unit=genomic_unit,
            project=project,
            annotation_extension="G026",
        ),
        config=_mirror_config(base_url),
    )


class TestAddBigwigUrlsCharacterization:
    """Pin the BigWig URL column exactly, row by row."""

    @staticmethod
    def _multi_project_bundle() -> R3ResourceBundle:
        mouse = _mock_resource(
            "count_files_gene_or_exon",
            organism=None,
            data_source="sra",
            project="SRP004",
        )

        return R3ResourceBundle(
            resources=[
                _count_resource("SRP001"),
                _count_resource(
                    "SRP001",
                    genomic_unit="exon",
                    base_url="https://mirror-b.example/r3/",
                ),
                _mock_resource("metadata_files", project="SRP003"),
                _count_resource(
                    "SRP002",
                    organism="mouse",
                    data_source="gtex",
                    base_url="https://mirror-c.example/",
                ),
                mouse,
            ]
        )

    def test_missing_category_does_not_change_an_integer_id(self) -> None:
        bundle = R3ResourceBundle([_count_resource("SRP001")])
        frame = pd.DataFrame(
            {"external_id": pd.Categorical([12, None]), "study": "SRP001"}
        )
        urls = bundle._add_bigwig_urls(frame)
        _assert_identical_frames(
            urls, _reference_add_bigwig_urls(bundle, frame)
        )
        assert (
            urls["BigWigURL"]
            .iloc[0]
            .endswith("/12/sra.base_sums.SRP001_12.ALL.bw")
        )

    @staticmethod
    def _metadata() -> pd.DataFrame:
        return pd.DataFrame(
            {
                "external_id": pd.array(
                    ["SRR1", "SRR2", "SRR3", "SRR4", " ", "", pd.NA, "SRR8"],
                    dtype="string",
                ),
                "study": pd.array(
                    [
                        "SRP001",
                        "SRP002",
                        "SRP003",
                        pd.NA,
                        "SRP001",
                        "SRP001",
                        "SRP001",
                        "SRP004",
                    ],
                    dtype="string",
                ),
                "score": np.arange(8.0),
            },
            index=[f"row{i}" for i in range(8)],
        )

    def test_urls_follow_each_samples_study(self) -> None:
        out = self._multi_project_bundle()._add_bigwig_urls(self._metadata())
        urls = out["BigWigURL"].tolist()
        assert urls[:2] == [
            "https://mirror-a.example/recount3/human/data_sources/sra/"
            "base_sums/01/SRP001/R1/sra.base_sums.SRP001_SRR1.ALL.bw",
            "https://mirror-c.example/mouse/data_sources/gtex/base_sums/"
            "02/SRP002/SR/gtex.base_sums.SRP002_SRR2.ALL.bw",
        ]
        # No count resource for SRP003; no study across several projects;
        # blank, empty or missing external IDs; SRP004 has no organism.
        assert all(pd.isna(url) for url in urls[2:])
        assert list(out.columns) == [
            "external_id",
            "study",
            "score",
            "BigWigURL",
        ]
        assert list(out.index) == [f"row{i}" for i in range(8)]

    def test_matches_the_reference(self) -> None:
        bundle = self._multi_project_bundle()
        _assert_identical_frames(
            bundle._add_bigwig_urls(self._metadata()),
            _reference_add_bigwig_urls(bundle, self._metadata()),
        )

    def test_a_single_project_supplies_samples_without_a_matching_study(
        self,
    ) -> None:
        bundle = R3ResourceBundle(
            resources=[
                _count_resource("SRP001"),
                _count_resource("SRP001", genomic_unit="exon"),
            ]
        )
        frame = pd.DataFrame(
            {
                "external_id": ["SRR1", "SRR2", "SRR3"],
                "study": pd.array([pd.NA, "SRP999", "SRP001"], dtype="string"),
            }
        )
        out = bundle._add_bigwig_urls(frame)
        assert all(
            url.endswith(f"SRP001_SRR{i}.ALL.bw")
            for i, url in enumerate(out["BigWigURL"], start=1)
        )
        _assert_identical_frames(out, _reference_add_bigwig_urls(bundle, frame))

    def test_without_a_study_column_only_a_single_project_resolves(
        self,
    ) -> None:
        frame = pd.DataFrame({"external_id": ["SRR1"]})
        single = R3ResourceBundle(resources=[_count_resource("SRP001")])
        multi = R3ResourceBundle(
            resources=[_count_resource("SRP001"), _count_resource("SRP002")]
        )
        assert single._add_bigwig_urls(frame)["BigWigURL"].iloc[0]
        assert pd.isna(multi._add_bigwig_urls(frame)["BigWigURL"].iloc[0])

    @pytest.mark.parametrize(
        ("source", "expected"),
        [
            ("s3://bucket/gtex/", "gtex"),
            ("gs://bucket/data/tcga", "tcga"),
            ("tcga", "tcga"),
            ("/", None),
            ("", None),
            (42, "sra"),
            (None, "sra"),
        ],
    )
    def test_file_source_overrides_the_data_source(
        self, source: Any, expected: str | None
    ) -> None:
        bundle = R3ResourceBundle(resources=[_count_resource("SRP001")])
        frame = pd.DataFrame(
            {
                "external_id": ["SRR1"],
                "RECOUNT_PROJECT__FILE_SOURCE": pd.Series(
                    [source], dtype=object
                ),
                "other__file_source": ["unused"],
            }
        )
        url = bundle._add_bigwig_urls(frame)["BigWigURL"].iloc[0]
        if expected is None:
            assert pd.isna(url)
        else:
            assert f"/data_sources/{expected}/" in url
            assert f"/{expected}.base_sums." in url
        _assert_identical_frames(
            bundle._add_bigwig_urls(frame),
            _reference_add_bigwig_urls(bundle, frame),
        )

    def test_identifiers_take_the_frames_common_row_dtype(self) -> None:
        """Rows are read at the dtype all columns share, as iterrows did.

        In an all-numeric frame that dtype is float, so an integer sample
        ID renders as ``12.0``; with any text column it stays ``12``.
        """
        bundle = R3ResourceBundle(resources=[_count_resource("SRP001")])
        numeric = pd.DataFrame({"external_id": [11, 12], "study": [1.5, 2.5]})
        out = bundle._add_bigwig_urls(numeric)
        assert (
            out["BigWigURL"]
            .iloc[1]
            .endswith("/.0/sra.base_sums.SRP001_12.0.ALL.bw")
        )
        _assert_identical_frames(
            out, _reference_add_bigwig_urls(bundle, numeric)
        )
        mixed = numeric.assign(note=["a", "b"])
        out = bundle._add_bigwig_urls(mixed)
        assert (
            out["BigWigURL"]
            .iloc[1]
            .endswith("/12/sra.base_sums.SRP001_12.ALL.bw")
        )
        _assert_identical_frames(out, _reference_add_bigwig_urls(bundle, mixed))

    def test_an_empty_frame_gains_an_empty_column(self) -> None:
        bundle = R3ResourceBundle(resources=[_count_resource("SRP001")])
        frame = pd.DataFrame(
            {"external_id": pd.array([], dtype="string"), "score": []}
        )
        out = bundle._add_bigwig_urls(frame)
        assert list(out.columns) == ["external_id", "score", "BigWigURL"]
        _assert_identical_frames(out, _reference_add_bigwig_urls(bundle, frame))

    def test_the_input_is_copied_not_modified(self) -> None:
        bundle = self._multi_project_bundle()
        frame = self._metadata()
        before = frame.copy()
        out = bundle._add_bigwig_urls(frame)
        _assert_identical_frames(frame, before)
        out.loc["row0", "score"] = -1.0
        assert frame.loc["row0", "score"] == 0.0

    def test_a_non_text_column_label_fails_only_when_consulted(self) -> None:
        bundle = R3ResourceBundle(resources=[_count_resource("SRP001")])
        unreached = pd.DataFrame({"external_id": [pd.NA], 5: ["x"]})
        assert pd.isna(bundle._add_bigwig_urls(unreached)["BigWigURL"].iloc[0])
        reached = pd.DataFrame({"external_id": ["SRR1"], 5: ["x"]})
        with pytest.raises(AttributeError):
            bundle._add_bigwig_urls(reached)

    def test_each_resource_uses_its_own_configuration(self) -> None:
        bundle = R3ResourceBundle(
            resources=[
                _count_resource("SRP001"),
                _count_resource(
                    "SRP002", base_url="https://mirror-z.example/x/"
                ),
            ]
        )
        frame = pd.DataFrame(
            {"external_id": ["SRR1", "SRR2"], "study": ["SRP001", "SRP002"]}
        )
        urls = bundle._add_bigwig_urls(frame)["BigWigURL"].tolist()
        assert urls[0].startswith("https://mirror-a.example/recount3/")
        assert urls[1].startswith("https://mirror-z.example/x/")

    @pytest.mark.parametrize("seed", range(60))
    def test_matches_the_reference_on_random_metadata(self, seed: int) -> None:
        rng = np.random.default_rng(seed)
        projects = ["SRP001", "SRP002", "SRP003"]
        resources = [
            _count_resource(projects[int(i)])
            for i in rng.integers(0, 3, size=int(rng.integers(0, 4)))
        ]
        bundle = R3ResourceBundle(resources=resources)
        size = int(rng.integers(0, 8))
        choices_ext = ["SRR1", "SRR2", "", "  ", None, "ERR9"]
        choices_study = [*projects, "SRP999", None]
        frame = pd.DataFrame(
            {
                "external_id": pd.array(
                    [choices_ext[int(i)] for i in rng.integers(0, 6, size)],
                    dtype="string",
                ),
                "study": pd.array(
                    [choices_study[int(i)] for i in rng.integers(0, 5, size)],
                    dtype="string",
                ),
                "x__file_source": pd.Series(
                    [
                        ["a/sra", "gtex", "/", None][int(i)]
                        for i in rng.integers(0, 4, size)
                    ],
                    dtype=object,
                ),
            }
        )
        if rng.random() < 0.3:
            frame = frame.drop(columns="x__file_source")
        _assert_identical_frames(
            bundle._add_bigwig_urls(frame),
            _reference_add_bigwig_urls(bundle, frame),
        )


class TestAddBigwigUrlsCandidates:
    def test_junction_count_files_identify_their_project(self) -> None:
        junctions = R3Resource(
            R3ResourceDescription(
                resource_type="count_files_junctions",
                organism="human",
                data_source="sra",
                project="SRP001",
                junction_type="ALL",
                junction_extension="MM",
            ),
            config=_mirror_config("https://mirror-a.example/recount3/"),
        )
        frame = pd.DataFrame({"external_id": ["SRR1"], "study": ["SRP001"]})
        url = R3ResourceBundle(resources=[junctions])._add_bigwig_urls(frame)
        assert url["BigWigURL"].iloc[0].endswith("SRP001_SRR1.ALL.bw")

    def test_metadata_files_never_supply_a_project(self) -> None:
        metadata = R3Resource(
            R3ResourceDescription(
                resource_type="metadata_files",
                organism="human",
                data_source="sra",
                table_name="recount_project",
                project="SRP003",
            )
        )
        bundle = R3ResourceBundle(
            resources=[
                metadata,
                _count_resource("SRP001"),
                _count_resource("SRP002"),
            ]
        )
        frame = pd.DataFrame(
            {"external_id": ["SRR1", "SRR2"], "study": ["SRP003", "SRP001"]}
        )
        urls = bundle._add_bigwig_urls(frame)["BigWigURL"]
        assert pd.isna(urls.iloc[0])
        assert urls.iloc[1].endswith("SRP001_SRR2.ALL.bw")

    def test_a_missing_study_matches_no_project_name(self) -> None:
        bundle = R3ResourceBundle(
            resources=[_count_resource("nan"), _count_resource("SRP002")]
        )
        frame = pd.DataFrame(
            {"external_id": ["SRR1"], "study": [np.nan], "note": ["x"]}
        )
        assert pd.isna(bundle._add_bigwig_urls(frame)["BigWigURL"].iloc[0])

    def test_numeric_studies_match_projects_by_their_text(self) -> None:
        bundle = R3ResourceBundle(
            resources=[_count_resource("7"), _count_resource("SRP002")]
        )
        frame = pd.DataFrame(
            {"external_id": ["SRR1"], "study": [7], "note": ["x"]}
        )
        url = bundle._add_bigwig_urls(frame)["BigWigURL"].iloc[0]
        assert url.endswith("sra.base_sums.7_SRR1.ALL.bw")

    def test_only_slashes_are_trimmed_from_the_source(self) -> None:
        """A source is the last path component, other characters kept.

        An unknown source is an error, so trimming more than slashes from
        ``tcga_/`` would wrongly turn it into the valid ``tcga``.
        """
        bundle = R3ResourceBundle(resources=[_count_resource("SRP001")])
        frame = pd.DataFrame(
            {"external_id": ["SRR1"], "file_source": ["s3://b/tcga_/"]}
        )
        with pytest.raises(ValueError, match="Invalid data_source: 'tcga_'"):
            bundle._add_bigwig_urls(frame)


class TestNormalizeSampleMetadataLookup:
    """Sample identifiers resolve to metadata rows exactly as before."""

    @staticmethod
    def _bundle(external: list[Any], rail: list[Any]) -> R3ResourceBundle:
        frame = pd.DataFrame(
            {
                "external_id": pd.array(external, dtype="string"),
                "rail_id": pd.array(rail, dtype="string"),
                "study": pd.array(["SRP001"] * len(external), dtype="string"),
                "value": np.arange(len(external), dtype=float),
            }
        )
        return R3ResourceBundle(resources=[_metadata_resource(frame)])

    def test_the_first_conflict_in_row_order_is_reported(self) -> None:
        bundle = self._bundle(["X", "Y"], ["Y", "X"])
        with pytest.raises(
            ValueError, match=r"^Ambiguous sample identifier 'Y' in metadata\.$"
        ):
            bundle._normalize_sample_metadata(sample_ids=["X", "Y"])

    def test_a_rail_id_conflict_is_reported(self) -> None:
        bundle = self._bundle(["A", "C"], ["B", "A"])
        with pytest.raises(ValueError, match=r"identifier 'A'"):
            bundle._normalize_sample_metadata(sample_ids=["A", "C"])

    def test_one_row_may_use_the_same_identifier_twice(self) -> None:
        bundle = self._bundle(["Z", "W"], ["Z", "7"])
        out = bundle._normalize_sample_metadata(sample_ids=["7", "Z"])
        assert list(out.index) == ["7", "Z"]
        assert list(out["external_id"]) == ["W", "Z"]
        assert list(out["recount_qc__value"]) == [1.0, 0.0]

    def test_samples_resolve_through_either_identifier(self) -> None:
        bundle = self._bundle(["SRR1", "SRR2", "SRR3"], ["11", "12", "13"])
        out = bundle._normalize_sample_metadata(sample_ids=["13", "SRR1", "12"])
        assert list(out.index) == ["13", "SRR1", "12"]
        assert list(out["external_id"]) == ["SRR3", "SRR1", "SRR2"]
        assert list(out["rail_id"]) == ["13", "11", "12"]
        assert list(out["recount_qc__value"]) == [2.0, 0.0, 1.0]

    def test_missing_identifiers_are_never_matched(self) -> None:
        bundle = self._bundle(["SRR1", pd.NA], [pd.NA, "12"])
        out = bundle._normalize_sample_metadata(
            sample_ids=["SRR1", "12", "<NA>"], metadata_join="outer"
        )
        assert list(out.index) == ["SRR1", "12", "<NA>"]
        assert list(out["recount_qc__value"].iloc[:2]) == [0.0, 1.0]
        assert pd.isna(out["recount_qc__value"].iloc[2])
        assert list(out["external_id"]) == ["SRR1", "12", "<NA>"]

    def test_identifiers_shared_across_projects_are_ambiguous(self) -> None:
        frames = [
            pd.DataFrame(
                {
                    "external_id": pd.array(["SRR1"], dtype="string"),
                    "rail_id": pd.array([rail], dtype="string"),
                    "study": pd.array([project], dtype="string"),
                }
            )
            for project, rail in (("SRP001", "1"), ("SRP002", "SRR1"))
        ]
        bundle = R3ResourceBundle(
            resources=[
                _metadata_resource(frames[0], project="SRP001"),
                _metadata_resource(frames[1], project="SRP002"),
            ]
        )
        with pytest.raises(ValueError, match=r"identifier 'SRR1'"):
            bundle._normalize_sample_metadata(sample_ids=["SRR1"])


class TestNormalizeSampleMetadataMessages:
    @staticmethod
    def _bundle(*tables: pd.DataFrame) -> R3ResourceBundle:
        return R3ResourceBundle(
            resources=[
                _metadata_resource(table, table_name=f"t{number}")
                for number, table in enumerate(tables)
            ]
        )

    @staticmethod
    def _table(external: list[Any], rail: list[str]) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "external_id": pd.array(external, dtype="string"),
                "rail_id": pd.array(rail, dtype="string"),
                "study": pd.array(["SRP001"] * len(rail), dtype="string"),
            }
        )

    def test_metadata_rows_without_counts(self) -> None:
        bundle = self._bundle(self._table(["S1", "S2"], ["1", "2"]))
        with pytest.raises(
            ValueError,
            match=(
                r"^Metadata contains samples missing from the counts "
                r"matrix\.$"
            ),
        ):
            bundle._normalize_sample_metadata(sample_ids=["S1"])

    def test_matched_rows_without_external_ids(self) -> None:
        bundle = self._bundle(self._table([pd.NA], ["1"]))
        with pytest.raises(
            ValueError,
            match=r"^Matched metadata rows have missing external_id values\.$",
        ):
            bundle._normalize_sample_metadata(sample_ids=["1"])

    def test_an_inner_join_matching_nothing(self) -> None:
        bundle = self._bundle(
            self._table(["S1"], ["1"]), self._table(["S9"], ["9"])
        )
        with pytest.raises(ValueError) as raised:
            bundle._normalize_sample_metadata(sample_ids=["S1"])
        assert str(raised.value) == (
            "No count sample matched the merged sample metadata. The "
            "metadata tables for this project share no "
            "rail_id/external_id/study rows, so the inner join produced no "
            "samples. Check that the tables describe the same samples, or "
            "pass metadata_join='outer' to keep every count sample and leave "
            "the unmatched metadata missing."
        )

    def test_no_requested_samples_and_no_shared_rows_is_empty(self) -> None:
        bundle = self._bundle(
            self._table(["S1"], ["1"]), self._table(["S9"], ["9"])
        )
        out = bundle._normalize_sample_metadata(sample_ids=[])
        assert out.shape[0] == 0
        assert "external_id" in out.columns


class TestConstructionStaysLinear:
    """Guard the choices that keep SE/RSE construction linear in its data.

    Each test forbids the operation whose cost grows faster than the input:
    a sparse ``take`` expands a column to the full target length, pairwise
    unions rehash the growing union per project, per-project scans of the
    bundle are quadratic in projects, and ``iterrows`` builds every metadata
    column for every sample.
    """

    @staticmethod
    def _forbid(owner: Any, name: str) -> Any:
        def fail(*args: Any, **kwargs: Any) -> None:
            raise AssertionError(f"{name} must not be called here")

        return patch.object(owner, name, fail)

    @pytest.mark.parametrize("dtype", [np.int64, np.float64])
    def test_loader_sparse_frames_are_never_expanded(self, dtype: Any) -> None:
        frames = TestMergeCountFramesCharacterization._three_loader_frames(
            dtype=dtype
        )
        if frames[0].dtypes.iloc[0].fill_value != 0:
            pytest.skip("This pandas version uses NaN fill for this dtype")
        expected = {
            policy: _reference_merge_count_frames(frames, policy)
            for policy in ("inner", "outer")
        }
        with self._forbid(pd.arrays.SparseArray, "take"):
            merged = {
                policy: bmod._merge_count_frames(frames, policy)
                for policy in expected
            }
        for policy, want in expected.items():
            _assert_identical_frames(merged[policy], want)

    @pytest.mark.parametrize("shared", [[], ["shared"]], ids=["empty", "one"])
    def test_small_inner_join_skips_full_sparse_conversion(
        self, shared: list[str]
    ) -> None:
        frames = []
        for project in range(2):
            labels = [f"P{project}F{i}" for i in range(128)] + shared
            frames.append(
                _loader_frame(
                    np.arange(len(labels)).reshape(-1, 1),
                    labels,
                    [f"S{project}"],
                )
            )
        expected = _reference_merge_count_frames(frames, "inner")
        with self._forbid(pd.DataFrame.sparse, "to_coo"):
            merged = bmod._merge_count_frames(frames, "inner")
        _assert_identical_frames(merged, expected)

    def test_the_outer_join_is_not_built_pairwise(self) -> None:
        frames = TestMergeCountFramesCharacterization._three_loader_frames()
        with self._forbid(pd.Index, "union"):
            out = bmod._merge_count_frames(frames, "outer")
        assert list(out.index) == ["a", "b", "c", "d", "e"]

    @pytest.mark.parametrize(
        "dtype", [object, pd.StringDtype(storage="python")]
    )
    @pytest.mark.parametrize("missing_first", [False, True])
    def test_overlapping_outer_join_does_not_append_numpy_labels(
        self, dtype: Any, missing_first: bool
    ) -> None:
        shared = pd.Index(
            [f"F{i}" for i in range(128)], dtype=dtype, name="feature"
        )
        subset = shared[1:].rename(None)
        indexes = [shared] * 40
        indexes.insert(0 if missing_first else len(indexes), subset)
        expected = (
            subset.union(shared, sort=False)
            if missing_first
            else shared.union(subset, sort=False)
        )
        with (
            self._forbid(pd.Index, "append"),
            patch.object(
                bmod,
                "_has_duplicate_features",
                wraps=bmod._has_duplicate_features,
            ) as duplicate_check,
        ):
            joined = bmod._join_feature_index(indexes, "outer")
        assert duplicate_check.call_count == 2
        pd.testing.assert_index_equal(joined, expected, exact=True)

    def test_projects_are_grouped_in_one_pass(self) -> None:
        resources = []
        for number in range(40):
            resources.extend(
                _project_resources(f"SRP{number:03d}", [f"S{number}"])
            )
        bundle = R3ResourceBundle(resources=resources)
        calls = 0
        original = bmod._project_key

        def counting(res: Any) -> Any:
            nonlocal calls
            calls += 1
            return original(res)

        with patch.object(bmod, "_project_key", counting):
            counts, *_ = bundle._prepare_experiment(
                genomic_unit="gene",
                annotation_extension=None,
                join_policy="inner",
                metadata_join="inner",
                autoload=False,
            )
        assert counts.shape[1] == 40
        # A scan of the bundle per project would make 40 * 80 calls.
        assert calls <= 10 * len(resources)

    def test_sample_loops_read_only_the_columns_they_need(self) -> None:
        bundle = TestAddBigwigUrlsCharacterization._multi_project_bundle()
        metadata = TestAddBigwigUrlsCharacterization._metadata()
        expected = _reference_add_bigwig_urls(bundle, metadata)
        lookup = TestNormalizeSampleMetadataLookup._bundle(
            ["SRR1", "SRR2"], ["11", "12"]
        )
        with self._forbid(pd.DataFrame, "iterrows"):
            urls = bundle._add_bigwig_urls(metadata)
            out = lookup._normalize_sample_metadata(sample_ids=["12", "SRR1"])
        _assert_identical_frames(urls, expected)
        assert list(out["external_id"]) == ["SRR2", "SRR1"]

    def test_dense_merges_skip_sparse_reindexing(self) -> None:
        first = _gene_df(features=["G1", "G2"], samples=["S1"])
        second = _gene_df(features=["G2", "G3"], samples=["S2"])
        with patch.object(
            bmod, "_reindex_sparse_frame", side_effect=AssertionError
        ):
            out = bmod._merge_count_frames([first, second], "outer")
        assert list(out.index) == ["G1", "G2", "G3"]


def test_docstring_examples_run() -> None:
    import doctest

    results = doctest.testmod(bmod, raise_on_error=False)
    assert results.attempted > 0
    assert results.failed == 0
