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

from collections.abc import Sequence
import logging
import pathlib
import tracemalloc
from typing import Any
from unittest import mock
import warnings

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import recount3._utils as _utils
from recount3.se import _assay_frame
from recount3.se import _expand_sra_attributes_df
from recount3.se import _resolve_annotation_extension
from recount3.se import _resolve_sra_attributes_column
from recount3.se import build_ranged_summarized_experiment
from recount3.se import build_summarized_experiment
from recount3.se import compute_read_counts
from recount3.se import compute_scale_factors
from recount3.se import compute_tpm
from recount3.se import create_ranged_summarized_experiment
from recount3.se import create_rse
from recount3.se import expand_sra_attributes
from recount3.se import is_paired_end
from recount3.se import to_anndata
from recount3.se import transform_counts

_MIRROR = (
    pathlib.Path(__file__).parent / "data" / "recount3_mirror" / "recount3"
).resolve()
# as_uri() rather than "file://" + str(path), which misplaces a Windows
# drive letter into the URL authority.
_MIRROR_URL = _MIRROR.as_uri() + "/"


class _MockBiocFrame:
    def __init__(self, df: pd.DataFrame) -> None:
        self._df = df

    def to_pandas(self) -> pd.DataFrame:
        return self._df.copy()

    @classmethod
    def from_pandas(cls, df: pd.DataFrame) -> _MockBiocFrame:
        return cls(df)


class _NoPandasBiocFrame:
    """BiocFrame whose to_pandas() method is absent."""


class _MockSE:
    def __init__(
        self,
        assays: dict[str, np.ndarray],
        col_data_df: pd.DataFrame | None = None,
        row_names: list[str] | None = None,
        col_names: list[str] | None = None,
    ) -> None:
        self._assays = assays
        _df = col_data_df if col_data_df is not None else pd.DataFrame()
        self._bioc_col_data = _MockBiocFrame(_df)
        self.row_names = row_names
        self.col_names = col_names

    @property
    def assay_names(self) -> list[str]:
        return list(self._assays.keys())

    def assay(self, name: str) -> np.ndarray:
        return self._assays[name]

    def get_assay(self, name: str) -> np.ndarray:
        return self._assays[name]

    @property
    def col_data(self) -> _MockBiocFrame:
        return self._bioc_col_data

    @property
    def column_data(self) -> _MockBiocFrame:
        return self._bioc_col_data

    def get_column_data(self) -> _MockBiocFrame:
        return self._bioc_col_data

    def get_row_names(self) -> list[str] | None:
        return self.row_names

    def get_column_names(self) -> list[str] | None:
        return self.col_names

    def set_column_data(self, new_bioc_col_data: _MockBiocFrame) -> _MockSE:
        obj = _MockSE.__new__(_MockSE)
        obj._assays = self._assays
        obj._bioc_col_data = new_bioc_col_data
        obj.row_names = self.row_names
        obj.col_names = self.col_names
        return obj


class _MockRSE(_MockSE):
    def __init__(
        self,
        assays: dict[str, np.ndarray],
        col_data_df: pd.DataFrame | None = None,
        row_names: list[str] | None = None,
        col_names: list[str] | None = None,
        row_ranges_widths: list[int] | None = None,
    ) -> None:
        super().__init__(assays, col_data_df, row_names, col_names)
        if row_ranges_widths is not None:
            rr = mock.MagicMock()
            rr.width = row_ranges_widths
            self.row_ranges = rr

    @property
    def width(self) -> list[int]:
        rr = getattr(self, "row_ranges", None)
        if rr is None:
            raise AttributeError("width requires row_ranges to be set")
        return rr.width

    def set_column_data(self, new_bioc_col_data: _MockBiocFrame) -> _MockRSE:
        obj = _MockRSE.__new__(_MockRSE)
        obj._assays = self._assays
        obj._bioc_col_data = new_bioc_col_data
        obj.row_names = self.row_names
        obj.col_names = self.col_names
        if hasattr(self, "row_ranges"):
            obj.row_ranges = self.row_ranges
        return obj


class _MockSEWithBadColData(_MockSE):
    @property
    def column_data(self) -> _NoPandasBiocFrame:
        return _NoPandasBiocFrame()


def _patch_biocpy() -> mock._patch:
    return mock.patch.multiple(
        _utils,
        get_biocframe_class=mock.DEFAULT,
        get_summarizedexperiment_class=mock.DEFAULT,
        get_ranged_summarizedexperiment_class=mock.DEFAULT,
    )


_META = pd.DataFrame(
    {
        "external_id": ["S1", "S2", "S3"],
        # avg_mapped_length: S1=100 (SE), S2=200 (PE), S3=100 (SE)
        "recount_qc.star.average_mapped_length": [100.0, 200.0, 100.0],
        "recount_seq_qc.avg_len": [100.0, 100.0, 100.0],
        # auc
        "recount_qc.bc_auc.all_reads_all_bases": [4e9, 8e9, 4e9],
        # mapped reads
        "recount_qc.star.all_mapped_reads": [1e6, 2e6, 1e6],
    }
)

_RAW = np.array([[100.0, 200.0, 300.0], [400.0, 500.0, 600.0]], dtype=float)


def _make_rse(
    counts: np.ndarray | None = None,
    col_data_df: pd.DataFrame | None = None,
    row_names: list[str] | None = None,
    col_names: list[str] | None = None,
    widths: list[int] | None = None,
    assay_name: str = "raw_counts",
) -> _MockRSE:
    """Convenience factory for _MockRSE objects."""
    arr = counts if counts is not None else _RAW.copy()
    df = col_data_df if col_data_df is not None else _META.copy()
    return _MockRSE(
        assays={assay_name: arr},
        col_data_df=df,
        row_names=row_names,
        col_names=col_names,
        row_ranges_widths=widths,
    )


# The dense matrix [[20, 0], [0, 7], [30, 0]] with feature 0's coverage in
# sample 0 split across two duplicate entries.
_DUPLICATE_DENSE = np.array([[20.0, 0.0], [0.0, 7.0], [30.0, 0.0]])


def _duplicate_entry_matrices() -> list[Any]:
    data = np.array([10.0, 10.0, 30.0, 7.0])
    sorted_csc = sparse.csc_matrix(
        (data, np.array([0, 0, 2, 1]), np.array([0, 3, 4])), shape=(3, 2)
    )
    unsorted_csc = sparse.csc_matrix(
        (
            np.array([30.0, 10.0, 10.0, 7.0]),
            np.array([2, 0, 0, 1]),
            np.array([0, 3, 4]),
        ),
        shape=(3, 2),
    )
    csr = sparse.csr_matrix(
        (
            np.array([10.0, 10.0, 7.0, 30.0]),
            np.array([0, 0, 1, 0]),
            np.array([0, 2, 3, 4]),
        ),
        shape=(3, 2),
    )
    coo = sparse.coo_matrix(
        ([10.0, 10.0, 7.0, 30.0], ([0, 0, 1, 2], [0, 0, 1, 0])), shape=(3, 2)
    )
    return [
        sorted_csc,
        unsorted_csc,
        csr,
        coo,
        sparse.csc_array(unsorted_csc),
        sparse.csr_array(csr),
    ]


_DUPLICATE_IDS = [
    "sorted_csc",
    "unsorted_csc",
    "csr",
    "coo",
    "csc_array",
    "csr_array",
]


class TestResolveSraAttributesColumn:
    def test_exact_match(self) -> None:
        cols = ["sra.sample_attributes", "other"]
        assert (
            _resolve_sra_attributes_column(cols, "sra.sample_attributes")
            == "sra.sample_attributes"
        )

    def test_dot_to_underscore_variant(self) -> None:
        cols = ["sra__sample_attributes"]
        assert (
            _resolve_sra_attributes_column(cols, "sra.sample_attributes")
            == "sra__sample_attributes"
        )

    def test_underscore_to_dot_variant(self) -> None:
        cols = ["sra.sample_attributes"]
        assert (
            _resolve_sra_attributes_column(cols, "sra__sample_attributes")
            == "sra.sample_attributes"
        )

    def test_no_match_returns_none(self) -> None:
        assert (
            _resolve_sra_attributes_column(
                ["unrelated"], "sra__sample_attributes"
            )
            is None
        )


class TestExpandSraAttributesDf:
    def test_missing_column_returns_copy(self) -> None:
        df = pd.DataFrame({"a": [1, 2]})
        result = _expand_sra_attributes_df(df)
        assert list(result.columns) == ["a"]
        assert result is not df  # copy

    def test_empty_dataframe_returns_copy(self) -> None:
        """Empty DataFrame -> `rows` list stays empty -> early return."""
        df = pd.DataFrame({"sra.sample_attributes": pd.Series([], dtype=str)})
        result = _expand_sra_attributes_df(df)
        assert list(result.columns) == ["sra.sample_attributes"]
        assert len(result) == 0

    def test_non_string_value_produces_empty_row(self) -> None:
        """Non-string cell (e.g., NaN) produces an empty attribute dict."""
        df = pd.DataFrame({"sra.sample_attributes": [float("nan")]})
        result = _expand_sra_attributes_df(df)
        assert set(result.columns) == {"sra.sample_attributes"}

    def test_empty_string_value_produces_empty_row(self) -> None:
        df = pd.DataFrame({"sra.sample_attributes": [""]})
        result = _expand_sra_attributes_df(df)
        assert set(result.columns) == {"sra.sample_attributes"}

    def test_basic_parsing(self) -> None:
        attrs = "age;;30|disease;;Control|tissue;;Brain"
        df = pd.DataFrame({"sra.sample_attributes": [attrs]})
        result = _expand_sra_attributes_df(df)
        assert "sra_attribute.age" in result.columns
        assert result["sra_attribute.age"].iloc[0] == "30"
        assert result["sra_attribute.disease"].iloc[0] == "Control"
        assert result["sra_attribute.tissue"].iloc[0] == "Brain"

    def test_entry_without_double_semicolon_skipped(self) -> None:
        """Entries lacking ';;' separator are silently dropped."""
        attrs = "age;;30|no_separator_here|tissue;;Brain"
        df = pd.DataFrame({"sra.sample_attributes": [attrs]})
        result = _expand_sra_attributes_df(df)
        assert "sra_attribute.no_separator_here" not in result.columns

    def test_empty_entry_in_pipe_split_skipped(self) -> None:
        """Empty string produced by trailing/double pipes is dropped."""
        attrs = "age;;30||tissue;;Brain"
        df = pd.DataFrame({"sra.sample_attributes": [attrs]})
        result = _expand_sra_attributes_df(df)
        assert "sra_attribute.age" in result.columns
        assert "sra_attribute.tissue" in result.columns

    def test_empty_key_after_strip_skipped(self) -> None:
        """Key that is blank after stripping ';;' is discarded."""
        attrs = " ;;value1|age;;30"
        df = pd.DataFrame({"sra.sample_attributes": [attrs]})
        result = _expand_sra_attributes_df(df)
        assert "sra_attribute." not in [c for c in result.columns]
        assert "sra_attribute.age" in result.columns

    def test_spaces_in_key_replaced_by_underscore(self) -> None:
        attrs = "cell type;;neuron"
        df = pd.DataFrame({"sra.sample_attributes": [attrs]})
        result = _expand_sra_attributes_df(df)
        assert "sra_attribute.cell_type" in result.columns
        assert result["sra_attribute.cell_type"].iloc[0] == "neuron"

    def test_value_with_double_semicolon_uses_first_split_only(self) -> None:
        """Value may itself contain ';;'; only the first split matters."""
        attrs = "key;;val;;extra"
        df = pd.DataFrame({"sra.sample_attributes": [attrs]})
        result = _expand_sra_attributes_df(df)
        assert result["sra_attribute.key"].iloc[0] == "val;;extra"

    def test_original_column_preserved(self) -> None:
        attrs = "age;;25"
        df = pd.DataFrame({"sra.sample_attributes": [attrs]})
        result = _expand_sra_attributes_df(df)
        assert "sra.sample_attributes" in result.columns
        assert result["sra.sample_attributes"].iloc[0] == attrs

    def test_custom_column_and_prefix(self) -> None:
        df = pd.DataFrame({"my_col": ["x;;1|y;;2"]})
        result = _expand_sra_attributes_df(
            df,
            sra_attributes_column="my_col",
            attribute_column_prefix="attr.",
        )
        assert "attr.x" in result.columns
        assert result["attr.x"].iloc[0] == "1"

    def test_multiple_rows(self) -> None:
        df = pd.DataFrame(
            {
                "sra.sample_attributes": [
                    "age;;30|disease;;Cancer",
                    "age;;45|disease;;Control",
                ]
            }
        )
        result = _expand_sra_attributes_df(df)
        assert result["sra_attribute.age"].tolist() == ["30", "45"]
        assert result["sra_attribute.disease"].tolist() == ["Cancer", "Control"]

    def test_index_preserved(self) -> None:
        df = pd.DataFrame(
            {"sra.sample_attributes": ["age;;30"]}, index=["custom_idx"]
        )
        result = _expand_sra_attributes_df(df)
        assert list(result.index) == ["custom_idx"]

    def test_namespaced_column_auto_detected(self) -> None:
        """Expand a namespaced 'sra__sample_attributes' column.

        Requesting the default dotted name also expands the namespaced column.
        """
        df = pd.DataFrame(
            {"sra__sample_attributes": ["age;;30|disease;;Cancer"]}
        )
        result = _expand_sra_attributes_df(df)
        assert "sra_attribute.age" in result.columns
        assert "sra_attribute.disease" in result.columns
        assert "sra__sample_attributes" in result.columns

    def test_explicit_namespaced_request_resolves_to_dotted(self) -> None:
        """Find the dotted column when the namespaced name is requested."""
        df = pd.DataFrame({"sra.sample_attributes": ["age;;42"]})
        result = _expand_sra_attributes_df(
            df, sra_attributes_column="sra__sample_attributes"
        )
        assert result["sra_attribute.age"].iloc[0] == "42"


class TestResolveAnnotationExtension:
    """Tests for the private _resolve_annotation_extension helper."""

    def test_junction_unit_returns_none(self) -> None:
        result = _resolve_annotation_extension(
            organism="human",
            genomic_unit="junction",
            annotation_label=None,
            annotation_extension=None,
        )
        assert result is None

    def test_explicit_extension_wins(self) -> None:
        result = _resolve_annotation_extension(
            organism="human",
            genomic_unit="gene",
            annotation_label="gencode_v26",
            annotation_extension="  G029  ",
        )
        assert result == "G029"

    def test_explicit_extension_strips_whitespace(self) -> None:
        result = _resolve_annotation_extension(
            organism="human",
            genomic_unit="exon",
            annotation_label=None,
            annotation_extension="  M023  ",
        )
        assert result == "M023"

    def test_annotation_label_calls_annotation_ext(self) -> None:
        with mock.patch(
            "recount3.search.annotation_ext", return_value="G029"
        ) as m:
            result = _resolve_annotation_extension(
                organism="human",
                genomic_unit="gene",
                annotation_label="gencode_v29",
                annotation_extension=None,
            )
        m.assert_called_once_with(organism="human", annotation="gencode_v29")
        assert result == "G029"

    def test_human_default_annotation(self) -> None:
        result = _resolve_annotation_extension(
            organism="human",
            genomic_unit="gene",
            annotation_label=None,
            annotation_extension=None,
        )
        assert result == "G026"

    def test_human_with_whitespace(self) -> None:
        result = _resolve_annotation_extension(
            organism="  Human  ",
            genomic_unit="gene",
            annotation_label=None,
            annotation_extension=None,
        )
        assert result == "G026"

    def test_mouse_default_annotation(self) -> None:
        result = _resolve_annotation_extension(
            organism="mouse",
            genomic_unit="exon",
            annotation_label=None,
            annotation_extension=None,
        )
        assert result == "M023"

    def test_unsupported_organism_raises(self) -> None:
        with pytest.raises(ValueError, match="Unsupported organism"):
            _resolve_annotation_extension(
                organism="zebrafish",
                genomic_unit="gene",
                annotation_label=None,
                annotation_extension=None,
            )

    def test_non_gene_exon_unit_returns_none(self) -> None:
        """Any unit besides 'gene'/'exon' (e.g. junction) returns None."""
        result = _resolve_annotation_extension(
            organism="human",
            genomic_unit="junction",
            annotation_label="ignored",
            annotation_extension="ignored",
        )
        assert result is None


class TestBuildSummarizedExperiment:
    def test_delegates_to_bundle(self) -> None:
        bundle = mock.MagicMock()
        bundle.to_summarized_experiment.return_value = "SE_sentinel"

        result = build_summarized_experiment(
            bundle,
            genomic_unit="gene",
            annotation_extension="G026",
            assay_name="raw_counts",
            join_policy="inner",
            metadata_join="inner",
            autoload=True,
        )

        bundle.to_summarized_experiment.assert_called_once_with(
            genomic_unit="gene",
            annotation_extension="G026",
            assay_name="raw_counts",
            join_policy="inner",
            metadata_join="inner",
            autoload=True,
        )
        assert result == "SE_sentinel"

    def test_normalizes_genomic_unit(self) -> None:
        bundle = mock.MagicMock()
        build_summarized_experiment(bundle, genomic_unit="Gene")
        call_kwargs = bundle.to_summarized_experiment.call_args.kwargs
        assert call_kwargs["genomic_unit"] == "gene"

    def test_invalid_genomic_unit_raises(self) -> None:
        bundle = mock.MagicMock()
        with pytest.raises(ValueError, match="Invalid genomic_unit"):
            build_summarized_experiment(bundle, genomic_unit="transcript")


class TestBuildRangedSummarizedExperiment:
    def test_delegates_to_bundle(self) -> None:
        bundle = mock.MagicMock()
        bundle.to_ranged_summarized_experiment.return_value = "RSE_sentinel"

        result = build_ranged_summarized_experiment(
            bundle,
            genomic_unit="exon",
            annotation_extension="G026",
            prefer_rr_junction_coordinates=False,
            assay_name="raw_counts",
            join_policy="outer",
            metadata_join="inner",
            autoload=False,
            allow_fallback_to_se=True,
        )

        bundle.to_ranged_summarized_experiment.assert_called_once_with(
            genomic_unit="exon",
            annotation_extension="G026",
            prefer_rr_junction_coordinates=False,
            assay_name="raw_counts",
            join_policy="outer",
            metadata_join="inner",
            autoload=False,
            allow_fallback_to_se=True,
        )
        assert result == "RSE_sentinel"

    def test_normalizes_genomic_unit(self) -> None:
        bundle = mock.MagicMock()
        build_ranged_summarized_experiment(bundle, genomic_unit="Junction")
        call_kwargs = bundle.to_ranged_summarized_experiment.call_args.kwargs
        assert call_kwargs["genomic_unit"] == "junction"


class TestCreateRangedSummarizedExperiment:
    def _make_bundle_mock(self) -> mock.MagicMock:
        bundle = mock.MagicMock()
        bundle.to_ranged_summarized_experiment.return_value = "RSE"
        return bundle

    def test_default_junction_extensions_for_gene(self) -> None:
        """Gene/exon assemblies do not need RR; default stays ('MM',)."""
        bundle = self._make_bundle_mock()
        with mock.patch("recount3.bundle.R3ResourceBundle") as mock_bundle:
            mock_bundle.discover.return_value = bundle
            create_ranged_summarized_experiment(
                project="SRP009615",
                genomic_unit="gene",
                organism="human",
            )
        call_kwargs = mock_bundle.discover.call_args.kwargs
        assert call_kwargs["junction_exts"] == ("MM",)

    def test_default_junction_extensions_for_junction(self) -> None:
        """Include the RR sidecar for junction-level RSE coordinates.

        The default includes 'RR' so create_rse works out of the box.
        """
        bundle = self._make_bundle_mock()
        with mock.patch("recount3.bundle.R3ResourceBundle") as mock_bundle:
            mock_bundle.discover.return_value = bundle
            create_ranged_summarized_experiment(
                project="SRP009615",
                genomic_unit="junction",
                organism="human",
            )
        call_kwargs = mock_bundle.discover.call_args.kwargs
        assert call_kwargs["junction_exts"] == ("MM", "RR")

    def test_explicit_junction_extensions(self) -> None:
        """Provided junction_extensions are tupled with strip()."""
        bundle = self._make_bundle_mock()
        with mock.patch("recount3.bundle.R3ResourceBundle") as mock_bundle:
            mock_bundle.discover.return_value = bundle
            create_ranged_summarized_experiment(
                project="SRP009615",
                genomic_unit="junction",
                organism="human",
                junction_extensions=["MM", " RR "],
            )
        call_kwargs = mock_bundle.discover.call_args.kwargs
        assert call_kwargs["junction_exts"] == ("MM", "RR")

    def test_annotations_default_when_ann_ext_none(self) -> None:
        """Junction unit -> ann_ext=None -> annotations='default'."""
        bundle = self._make_bundle_mock()
        with mock.patch("recount3.bundle.R3ResourceBundle") as mock_bundle:
            mock_bundle.discover.return_value = bundle
            create_ranged_summarized_experiment(
                project="SRP009615",
                genomic_unit="junction",
                organism="human",
            )
        call_kwargs = mock_bundle.discover.call_args.kwargs
        assert call_kwargs["annotations"] == "default"

    def test_annotations_tuple_when_ann_ext_provided(self) -> None:
        """Gene unit with default human -> ann_ext='G026' -> annotations=('G026',)."""
        bundle = self._make_bundle_mock()
        with mock.patch("recount3.bundle.R3ResourceBundle") as mock_bundle:
            mock_bundle.discover.return_value = bundle
            create_ranged_summarized_experiment(
                project="SRP009615",
                genomic_unit="gene",
                organism="human",
            )
        call_kwargs = mock_bundle.discover.call_args.kwargs
        assert call_kwargs["annotations"] == ("G026",)

    def test_all_kwargs_forwarded_to_rse(self) -> None:
        bundle = self._make_bundle_mock()
        with mock.patch("recount3.bundle.R3ResourceBundle") as mock_bundle:
            mock_bundle.discover.return_value = bundle
            result = create_ranged_summarized_experiment(
                project="SRP009615",
                genomic_unit="gene",
                organism="human",
                data_source="sra",
                include_metadata=False,
                include_bigwig=True,
                assay_name="counts",
                join_policy="outer",
                autoload=False,
                allow_fallback_to_se=True,
                strict=False,
            )
        rse_kwargs = bundle.to_ranged_summarized_experiment.call_args.kwargs
        assert rse_kwargs["assay_name"] == "counts"
        assert rse_kwargs["join_policy"] == "outer"
        assert rse_kwargs["autoload"] is False
        assert rse_kwargs["allow_fallback_to_se"] is True
        assert result == "RSE"


class TestCreateRse:
    def test_delegates_to_create_rse(self) -> None:
        sentinel = object()
        with mock.patch(
            "recount3.se.create_ranged_summarized_experiment",
            return_value=sentinel,
        ) as m:
            result = create_rse(
                project="SRP009615",
                genomic_unit="gene",
                organism="human",
                data_source="sra",
                annotation_label="gencode_v26",
                annotation_extension="G026",
                junction_type="ALL",
                junction_extensions=["MM"],
                include_metadata=True,
                include_bigwig=False,
                assay_name="raw_counts",
                join_policy="inner",
                autoload=True,
                allow_fallback_to_se=False,
                strict=True,
            )
        m.assert_called_once()
        assert result is sentinel

    def test_passes_all_kwargs_through(self) -> None:
        with mock.patch(
            "recount3.se.create_ranged_summarized_experiment",
        ) as m:
            create_rse(project="P1", junction_extensions=["RR"])
        call_kwargs = m.call_args.kwargs
        assert call_kwargs["project"] == "P1"
        assert call_kwargs["junction_extensions"] == ["RR"]


class TestCreateRseAgainstLocalMirror:
    """End-to-end runs against the mirrored recount3 tree in tests/data."""

    @staticmethod
    def _use_mirror(
        monkeypatch: pytest.MonkeyPatch, cache_dir: pathlib.Path
    ) -> None:
        monkeypatch.setenv("RECOUNT3_URL", _MIRROR_URL)
        monkeypatch.setenv("RECOUNT3_CACHE_DIR", str(cache_dir))

    @staticmethod
    def _warnings(caplog: pytest.LogCaptureFixture) -> list[str]:
        return [
            rec.getMessage()
            for rec in caplog.records
            if rec.levelno >= logging.WARNING
        ]

    @pytest.mark.requires_biocpy
    def test_cold_cache_reports_only_the_real_reason(
        self,
        tmp_path: pathlib.Path,
        monkeypatch: pytest.MonkeyPatch,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """The regression: a cold cache used to warn about FileNotFoundError."""
        self._use_mirror(monkeypatch, tmp_path / "cache")

        with caplog.at_level(logging.DEBUG):
            experiment = create_rse(
                project="SRP014565",
                organism="human",
                annotation_extension="G026",
                allow_fallback_to_se=True,
            )

        assert experiment.shape[1] > 0
        messages = self._warnings(caplog)
        assert len(messages) == 1
        assert "Falling back to a plain SummarizedExperiment" in messages[0]
        assert "does not cover every counted feature" in messages[0]
        assert "FileNotFoundError" not in messages[0]

    @pytest.mark.requires_biocpy
    def test_warm_cache_reports_exactly_the_same_thing(
        self,
        tmp_path: pathlib.Path,
        monkeypatch: pytest.MonkeyPatch,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Behavior must not depend on whether the cache is already filled."""
        self._use_mirror(monkeypatch, tmp_path / "cache")

        kwargs = {
            "project": "SRP014565",
            "organism": "human",
            "annotation_extension": "G026",
            "allow_fallback_to_se": True,
        }
        with caplog.at_level(logging.DEBUG):
            create_rse(**kwargs)
        first = self._warnings(caplog)

        caplog.clear()
        with caplog.at_level(logging.DEBUG):
            create_rse(**kwargs)
        second = self._warnings(caplog)

        assert first == second

    @pytest.mark.requires_biocpy
    def test_without_fallback_the_error_explains_the_option(
        self, tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        self._use_mirror(monkeypatch, tmp_path / "cache")

        with pytest.raises(ValueError) as excinfo:
            create_rse(
                project="SRP014565",
                organism="human",
                annotation_extension="G026",
                allow_fallback_to_se=False,
            )

        message = str(excinfo.value)
        assert "Could not derive genomic ranges" in message
        assert "does not cover every counted feature" in message
        # The option changes the returned object; it is not a retry or a fix.
        assert "no genomic ranges" in message
        assert "neither retries" in message


class TestExpandSraAttributes:
    def test_dataframe_missing_column_returns_copy(self) -> None:
        df = pd.DataFrame({"other": [1, 2]})
        result = expand_sra_attributes(df)
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ["other"]

    def test_dataframe_with_column_expands(self) -> None:
        df = pd.DataFrame(
            {"sra.sample_attributes": ["age;;30|disease;;Cancer"]}
        )
        result = expand_sra_attributes(df)
        assert isinstance(result, pd.DataFrame)
        assert "sra_attribute.age" in result.columns

    def test_dataframe_custom_params(self) -> None:
        df = pd.DataFrame({"my_col": ["k;;v"]})
        result = expand_sra_attributes(
            df,
            sra_attributes_column="my_col",
            attribute_column_prefix="p.",
        )
        assert "p.k" in result.columns

    def _biocpy_patches(self) -> mock._patch:
        return mock.patch.multiple(
            _utils,
            get_biocframe_class=mock.MagicMock(return_value=_MockBiocFrame),
            get_summarizedexperiment_class=mock.MagicMock(return_value=_MockSE),
            get_ranged_summarizedexperiment_class=mock.MagicMock(
                return_value=_MockRSE
            ),
        )

    def test_se_with_column_returns_expanded_se(self) -> None:
        col_df = pd.DataFrame(
            {"sra.sample_attributes": ["age;;25|disease;;Control"]},
            index=["S1"],
        )
        se = _MockSE(assays={}, col_data_df=col_df)

        with self._biocpy_patches():
            result = expand_sra_attributes(se)

        assert isinstance(result, _MockSE)
        expanded = result.column_data.to_pandas()
        assert "sra_attribute.age" in expanded.columns
        assert expanded["sra_attribute.age"].iloc[0] == "25"

    def test_rse_with_column_returns_expanded_rse(self) -> None:
        col_df = pd.DataFrame({"sra.sample_attributes": ["tissue;;brain"]})
        rse = _MockRSE(assays={}, col_data_df=col_df)

        with self._biocpy_patches():
            result = expand_sra_attributes(rse)

        assert isinstance(result, _MockRSE)
        expanded = result.column_data.to_pandas()
        assert "sra_attribute.tissue" in expanded.columns

    def test_se_missing_column_returns_unchanged_with_warning(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        col_df = pd.DataFrame({"other_col": [1]})
        se = _MockSE(assays={}, col_data_df=col_df)

        with self._biocpy_patches():
            with caplog.at_level(logging.WARNING):
                result = expand_sra_attributes(se)

        assert result is se
        assert (
            "not present" in caplog.text.lower() or "not present" in caplog.text
        )

    def test_se_with_bad_col_data_raises_attribute_error(self) -> None:
        se = _MockSEWithBadColData(assays={})

        with self._biocpy_patches():
            with pytest.raises(AttributeError, match="to_pandas"):
                expand_sra_attributes(se)

    def test_invalid_type_raises_type_error(self) -> None:
        with self._biocpy_patches():
            with pytest.raises(TypeError, match="pandas DataFrame or a BiocPy"):
                expand_sra_attributes(42)  # type: ignore[arg-type]

    def test_custom_sra_column_and_prefix_forwarded(self) -> None:
        col_df = pd.DataFrame({"my_attrs": ["x;;1"]})
        se = _MockSE(assays={}, col_data_df=col_df)

        with self._biocpy_patches():
            result = expand_sra_attributes(
                se,
                sra_attributes_column="my_attrs",
                attribute_column_prefix="pfx.",
            )

        expanded = result.column_data.to_pandas()
        assert "pfx.x" in expanded.columns


class TestComputeReadCounts:
    def _patch_rse_class(self) -> mock._patch:
        return mock.patch.object(
            _utils,
            "get_ranged_summarizedexperiment_class",
            return_value=_MockRSE,
        )

    def test_non_rse_raises_type_error(self) -> None:
        with self._patch_rse_class():
            with pytest.raises(
                TypeError,
                match=(
                    r"^rse must be a RangedSummarizedExperiment to compute "
                    r"read counts\.$"
                ),
            ):
                compute_read_counts("not_an_rse")

    def test_round_to_integers_not_bool_raises_type_error(self) -> None:
        rse = _make_rse()
        with self._patch_rse_class():
            with pytest.raises(
                TypeError, match=r"^round_to_integers must be a bool\.$"
            ):
                compute_read_counts(rse, round_to_integers=1)  # type: ignore[arg-type]

    def test_missing_avg_length_column_raises_value_error(self) -> None:
        meta = pd.DataFrame({"external_id": ["S1"]})
        rse = _make_rse(counts=np.array([[10.0]]), col_data_df=meta)
        with self._patch_rse_class():
            with pytest.raises(
                ValueError,
                match=(
                    r"^Required metadata column not found in col_data: "
                    r"'recount_qc\.star\.average_mapped_length'\.$"
                ),
            ):
                compute_read_counts(rse)

    def test_rounds_by_default(self) -> None:
        rse = _make_rse()
        with self._patch_rse_class():
            result = compute_read_counts(rse)
        # S2: 500 / 200 = 2.5 rounds half to even.
        assert result.iloc[1, 1] == 2.0

    def test_string_in_avg_length_raises_value_error(self) -> None:
        meta = _META.copy()
        meta["recount_qc.star.average_mapped_length"] = ["abc", 100.0, 100.0]
        rse = _make_rse(col_data_df=meta)
        with self._patch_rse_class():
            with pytest.raises(ValueError, match="must be numeric"):
                compute_read_counts(rse)

    def test_bytes_in_avg_length_raises_value_error(self) -> None:
        meta = _META.copy()
        meta["recount_qc.star.average_mapped_length"] = [b"bad", 100.0, 100.0]
        rse = _make_rse(col_data_df=meta)
        with self._patch_rse_class():
            with pytest.raises(ValueError, match="must be numeric"):
                compute_read_counts(rse)

    def test_non_numeric_non_string_in_avg_length_raises(self) -> None:
        """Object that is not str/bytes/Real and not NaN triggers the last branch."""
        meta = _META.copy()
        meta["recount_qc.star.average_mapped_length"] = meta[
            "recount_qc.star.average_mapped_length"
        ].astype(object)
        meta.at[0, "recount_qc.star.average_mapped_length"] = complex(1, 2)
        rse = _make_rse(col_data_df=meta)
        with self._patch_rse_class():
            with pytest.raises(ValueError, match="must be numeric"):
                compute_read_counts(rse)

    def test_nan_in_avg_length_is_skipped_without_error(self) -> None:
        meta = _META.copy()
        meta = meta.copy()
        meta["recount_qc.star.average_mapped_length"] = meta[
            "recount_qc.star.average_mapped_length"
        ].astype(float)
        meta.at[0, "recount_qc.star.average_mapped_length"] = float("nan")
        rse = _make_rse(col_data_df=meta)
        with self._patch_rse_class():
            result = compute_read_counts(rse, round_to_integers=False)
        assert isinstance(result, pd.DataFrame)
        # NaN avg_length -> NaN read counts
        assert result.iloc[:, 0].isna().all()

    def test_1d_assay_raises_value_error(self) -> None:
        rse = _make_rse(counts=np.array([1.0, 2.0, 3.0]))  # 1D
        with self._patch_rse_class():
            with pytest.raises(ValueError, match="2D matrix"):
                compute_read_counts(rse)

    def test_shape_mismatch_raises_value_error(self) -> None:
        rse = _make_rse(counts=np.ones((2, 4)))
        with self._patch_rse_class():
            with pytest.raises(ValueError, match="Mismatch"):
                compute_read_counts(rse)

    def test_basic_computation_rounded(self) -> None:
        """read_counts = raw / avg_len, then rint."""
        rse = _make_rse()
        with self._patch_rse_class():
            result = compute_read_counts(rse, round_to_integers=True)

        assert isinstance(result, pd.DataFrame)
        assert result.shape == (2, 3)
        # S1: avg=100 -> [100/100, 400/100] = [1, 4]
        np.testing.assert_array_equal(result.iloc[:, 0].values, [1.0, 4.0])
        # S2: avg=200 -> [200/200, 500/200] = [1.0, 2.5] -> rint -> [1.0, 2.0]
        np.testing.assert_array_equal(result.iloc[:, 1].values, [1.0, 2.0])
        # S3: avg=100 -> [300/100, 600/100] = [3.0, 6.0]
        np.testing.assert_array_equal(result.iloc[:, 2].values, [3.0, 6.0])

    def test_basic_computation_unrounded(self) -> None:
        rse = _make_rse()
        with self._patch_rse_class():
            result = compute_read_counts(rse, round_to_integers=False)

        # S2: 500/200 = 2.5 (not rounded)
        assert result.iloc[1, 1] == pytest.approx(2.5)

    def test_with_row_and_col_names(self) -> None:
        rse = _make_rse(
            row_names=["GENE1", "GENE2"],
            col_names=["S1", "S2", "S3"],
        )
        with self._patch_rse_class():
            result = compute_read_counts(rse)

        assert list(result.index) == ["GENE1", "GENE2"]
        assert list(result.columns) == ["S1", "S2", "S3"]

    def test_without_row_and_col_names_uses_none_index(self) -> None:
        rse = _make_rse(row_names=None, col_names=None)
        with self._patch_rse_class():
            result = compute_read_counts(rse)

        assert isinstance(result.index, pd.RangeIndex)
        assert isinstance(result.columns, pd.RangeIndex)

    def test_fallback_assay_name(self) -> None:
        """'counts' assay accepted as fallback."""
        rse = _MockRSE(
            assays={"counts": _RAW.copy()},
            col_data_df=_META.copy(),
        )
        with self._patch_rse_class():
            result = compute_read_counts(rse)
        assert result.shape == (2, 3)


# ===========================================================================
# compute_tpm
# ===========================================================================

_READ_LENGTH_COLUMN = "recount_qc.star.average_mapped_length"


class _MockRowData:
    """Row metadata exposing the BiocFrame members compute_tpm relies on."""

    def __init__(self, columns: dict[str, Any]) -> None:
        self._columns = columns

    @property
    def column_names(self) -> list[str]:
        return list(self._columns)

    def __getitem__(self, name: str) -> Any:
        return self._columns[name]


class _MockRSEWithRowData(_MockRSE):
    """RSE mock whose ``get_row_data`` returns a configurable object."""

    def __init__(self, *args: Any, row_data: Any = None, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self._row_data = row_data

    def get_row_data(self) -> Any:
        return self._row_data


class _MockRSEWithBrokenRowData(_MockRSE):
    """RSE mock whose ``get_row_data`` raises like a half-built object."""

    def get_row_data(self) -> Any:
        raise AttributeError("row data is unavailable")


def _make_tpm_rse(
    counts: Any,
    read_lengths: Sequence[Any],
    *,
    widths: Sequence[Any] | None = None,
    bp_length: Sequence[Any] | None = None,
    row_data: Any = "auto",
    row_names: list[str] | None = None,
    col_names: list[str] | None = None,
    assay_name: str = "raw_counts",
    read_length_column: str = _READ_LENGTH_COLUMN,
) -> _MockRSE:
    """Build an RSE mock for compute_tpm.

    ``row_data="auto"`` attaches row metadata carrying ``bp_length`` when it
    is given and omits ``get_row_data`` entirely otherwise; pass an explicit
    object (or ``None``) to control what ``get_row_data`` returns.
    """
    col_data = pd.DataFrame({read_length_column: list(read_lengths)})
    kwargs: dict[str, Any] = dict(
        assays={assay_name: counts},
        col_data_df=col_data,
        row_names=row_names,
        col_names=col_names,
        row_ranges_widths=None if widths is None else list(widths),
    )
    if row_data == "auto":
        if bp_length is None:
            return _MockRSE(**kwargs)
        row_data = _MockRowData({"bp_length": list(bp_length)})
    return _MockRSEWithRowData(**kwargs, row_data=row_data)


def _reference_tpm(
    counts: Any,
    widths: Sequence[Any],
    read_lengths: Sequence[Any],
) -> np.ndarray:
    """Evaluate the documented TPM recipe literally, step by step.

    Column sums skip missing values, as pandas does by default.
    """
    dense = counts.toarray() if sparse.issparse(counts) else counts
    dense = np.asarray(dense, dtype=float)
    widths_arr = np.asarray(widths, dtype=float)
    lengths_arr = np.asarray(read_lengths, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        reads = dense / lengths_arr[np.newaxis, :]
        rpk = reads / (widths_arr[:, np.newaxis] / 1000.0)
        scale = np.nansum(rpk, axis=0) / 1_000_000.0
        return rpk / scale[np.newaxis, :]


def _assert_tpm_equal(result: pd.DataFrame, expected: np.ndarray) -> None:
    """Compare TPM values, treating NaN as equal to NaN."""
    np.testing.assert_allclose(
        result.to_numpy(dtype=float), expected, rtol=1e-12, atol=0.0
    )


def _random_tpm_case(
    seed: int,
    n_features: int = 40,
    n_samples: int = 7,
    density: float = 0.3,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Draw integer coverage sums, feature lengths, and read lengths."""
    rng = np.random.default_rng(seed)
    counts = rng.integers(1, 10_000, size=(n_features, n_samples))
    counts = counts * (rng.random((n_features, n_samples)) < density)
    widths = rng.integers(50, 20_000, size=n_features).astype(float)
    read_lengths = rng.uniform(40.0, 300.0, size=n_samples)
    return counts.astype(np.int64), widths, read_lengths


def _patch_tpm_rse_class() -> mock._patch:
    return mock.patch.object(
        _utils,
        "get_ranged_summarizedexperiment_class",
        return_value=_MockRSE,
    )


class TestComputeTpmValidation:
    """compute_tpm rejects inputs it cannot normalize."""

    def test_non_rse_raises_type_error(self) -> None:
        with _patch_tpm_rse_class():
            with pytest.raises(
                TypeError,
                match=(
                    r"^compute_tpm requires a RangedSummarizedExperiment with "
                    r"defined rowRanges to determine feature lengths\.$"
                ),
            ):
                compute_tpm("not_rse")

    def test_type_check_uses_the_ranged_experiment_class(self) -> None:
        """A plain SE mock is not an instance of the RSE class."""
        se_obj = _MockSE(
            assays={"raw_counts": _RAW.copy()}, col_data_df=_META.copy()
        )
        with _patch_tpm_rse_class():
            with pytest.raises(TypeError):
                compute_tpm(se_obj)

    def test_missing_row_ranges_raises_value_error(self) -> None:
        rse = _make_tpm_rse(_RAW.copy(), [100.0, 200.0, 100.0])
        with _patch_tpm_rse_class():
            with pytest.raises(
                ValueError,
                match=r"^Could not determine feature widths from rowRanges\.$",
            ) as excinfo:
                compute_tpm(rse)
        assert isinstance(excinfo.value.__cause__, AttributeError)

    def test_row_data_attribute_error_raises_value_error(self) -> None:
        rse = _MockRSEWithBrokenRowData(
            assays={"raw_counts": _RAW.copy()},
            col_data_df=_META.copy(),
            row_ranges_widths=[1000, 2000],
        )
        with _patch_tpm_rse_class():
            with pytest.raises(ValueError, match="feature widths") as excinfo:
                compute_tpm(rse)
        assert isinstance(excinfo.value.__cause__, AttributeError)

    @pytest.mark.parametrize("n_widths", [1, 3])
    def test_width_count_mismatch_raises_value_error(
        self, n_widths: int
    ) -> None:
        rse = _make_tpm_rse(
            _RAW.copy(), [100.0, 200.0, 100.0], widths=[1000] * n_widths
        )
        with _patch_tpm_rse_class():
            with pytest.raises(
                ValueError,
                match=(
                    r"^Expected one feature length per row of the "
                    rf"'raw_counts' assay \(2\); got {n_widths}\.$"
                ),
            ):
                compute_tpm(rse)

    def test_bp_length_count_mismatch_raises_value_error(self) -> None:
        rse = _make_tpm_rse(
            _RAW.copy(), [100.0, 200.0, 100.0], bp_length=[1000.0]
        )
        with _patch_tpm_rse_class():
            with pytest.raises(ValueError, match="got 1"):
                compute_tpm(rse)

    def test_two_dimensional_lengths_raise_value_error(self) -> None:
        rse = _make_tpm_rse(
            _RAW.copy(), [100.0, 200.0, 100.0], bp_length=[[1000.0, 2000.0]]
        )
        with _patch_tpm_rse_class():
            with pytest.raises(ValueError, match="got 2"):
                compute_tpm(rse)

    def test_missing_counts_assay_raises_value_error(self) -> None:
        rse = _make_tpm_rse(
            _RAW.copy(),
            [100.0, 200.0, 100.0],
            widths=[1000, 2000],
            assay_name="tpm",
        )
        with _patch_tpm_rse_class():
            with pytest.raises(ValueError, match="'raw_counts' assay"):
                compute_tpm(rse)

    def test_missing_read_length_column_raises_value_error(self) -> None:
        rse = _make_tpm_rse(
            _RAW.copy(),
            [100.0, 200.0, 100.0],
            widths=[1000, 2000],
            read_length_column="unrelated",
        )
        with _patch_tpm_rse_class():
            with pytest.raises(
                ValueError, match="Required metadata column not found"
            ):
                compute_tpm(rse)

    @pytest.mark.parametrize(
        ("bad", "message"),
        [
            ("100", r"must be numeric; found str value '100'\.$"),
            (b"100", r"must be numeric; found bytes value b'100'\.$"),
            # pandas stores a Python complex as NumPy's complex128.
            (complex(1, 2), r"must be numeric; found complex128\.$"),
        ],
    )
    def test_non_numeric_read_length_raises_value_error(
        self, bad: Any, message: str
    ) -> None:
        rse = _make_tpm_rse(
            _RAW.copy(), [bad, 200.0, 100.0], widths=[1000, 2000]
        )
        with _patch_tpm_rse_class():
            with pytest.raises(
                ValueError,
                match=(
                    r"^Metadata column "
                    r"'recount_qc\.star\.average_mapped_length' " + message
                ),
            ):
                compute_tpm(rse)

    def test_non_numeric_read_length_after_a_missing_one_is_rejected(
        self,
    ) -> None:
        rse = _make_tpm_rse(
            _RAW.copy(), [np.nan, "abc", 100.0], widths=[1000, 2000]
        )
        with _patch_tpm_rse_class():
            with pytest.raises(ValueError, match="found str value 'abc'"):
                compute_tpm(rse)

    def test_sample_count_mismatch_raises_value_error(self) -> None:
        rse = _make_tpm_rse(_RAW.copy(), [100.0, 200.0], widths=[1000, 2000])
        with _patch_tpm_rse_class():
            with pytest.raises(ValueError, match="Mismatch"):
                compute_tpm(rse)

    def test_one_dimensional_assay_raises_value_error(self) -> None:
        rse = _make_tpm_rse(
            np.array([1.0, 2.0, 3.0]), [100.0, 200.0, 100.0], widths=[1000]
        )
        with _patch_tpm_rse_class():
            with pytest.raises(ValueError, match="2D matrix"):
                compute_tpm(rse)


class TestComputeTpmValues:
    """compute_tpm follows the documented four-step recipe."""

    def test_matches_hand_computed_values(self) -> None:
        # Reads = coverage / read length, RPK = reads / kb, TPM = RPK share.
        counts = np.array([[100.0, 0.0, 300.0], [400.0, 500.0, 50.0]])
        rse = _make_tpm_rse(counts, [100.0, 200.0, 50.0], widths=[1000, 2000])
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)

        # S1: RPK = [1.0, 2.0] -> [1/3, 2/3] of a million.
        # S2: RPK = [0.0, 1.25] -> [0, 1e6].
        # S3: RPK = [6.0, 0.5] -> [12/13, 1/13] of a million.
        expected = np.array(
            [
                [1e6 / 3, 0.0, 12e6 / 13],
                [2e6 / 3, 1e6, 1e6 / 13],
            ]
        )
        _assert_tpm_equal(result, expected)

    def test_returns_float64_dataframe(self) -> None:
        rse = _make_tpm_rse(
            _RAW.astype(np.int64), [100, 200, 100], widths=[1000, 2000]
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        assert isinstance(result, pd.DataFrame)
        assert result.shape == (2, 3)
        assert all(dtype == np.float64 for dtype in result.dtypes)

    @pytest.mark.parametrize("seed", range(5))
    def test_matches_reference_on_random_dense_input(self, seed: int) -> None:
        counts, widths, read_lengths = _random_tpm_case(seed)
        rse = _make_tpm_rse(counts, read_lengths, widths=widths)
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(result, _reference_tpm(counts, widths, read_lengths))

    @pytest.mark.parametrize("seed", range(3))
    def test_every_sample_sums_to_one_million(self, seed: int) -> None:
        counts, widths, read_lengths = _random_tpm_case(seed, density=0.8)
        rse = _make_tpm_rse(counts, read_lengths, widths=widths)
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        np.testing.assert_allclose(
            result.sum(axis=0).to_numpy(), 1_000_000.0, rtol=1e-12
        )

    def test_values_are_relative_to_feature_length(self) -> None:
        """Equal coverage on a feature twice as long gets half the TPM."""
        counts = np.array([[500.0], [500.0], [0.0]])
        rse = _make_tpm_rse(counts, [100.0], widths=[1000, 2000, 10])
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(result, np.array([[2e6 / 3], [1e6 / 3], [0.0]]))

    def test_read_length_does_not_change_finite_results(self) -> None:
        """Per-sample read length cancels in the normalization."""
        counts, widths, read_lengths = _random_tpm_case(11)
        base = _make_tpm_rse(counts, read_lengths, widths=widths)
        rescaled = _make_tpm_rse(
            counts, read_lengths * np.linspace(0.01, 7.0, 7), widths=widths
        )
        with _patch_tpm_rse_class():
            np.testing.assert_allclose(
                compute_tpm(base).to_numpy(),
                compute_tpm(rescaled).to_numpy(),
                rtol=1e-12,
            )

    def test_negative_read_length_cancels_like_a_positive_one(self) -> None:
        counts = _RAW.copy()
        rse = _make_tpm_rse(counts, [-100.0, 200.0, 100.0], widths=[1000, 2000])
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(
            result, _reference_tpm(counts, [1000, 2000], [100.0, 200.0, 100.0])
        )

    def test_integer_read_lengths_are_accepted(self) -> None:
        rse = _make_tpm_rse(_RAW.copy(), [100, 200, 100], widths=[1000, 2000])
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(
            result, _reference_tpm(_RAW, [1000, 2000], [100, 200, 100])
        )

    def test_names_are_preserved(self) -> None:
        rse = _make_tpm_rse(
            _RAW.copy(),
            [100.0, 200.0, 100.0],
            widths=[1000, 2000],
            row_names=["GENE1", "GENE2"],
            col_names=["S1", "S2", "S3"],
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        assert list(result.index) == ["GENE1", "GENE2"]
        assert list(result.columns) == ["S1", "S2", "S3"]

    def test_duplicate_names_keep_positional_values(self) -> None:
        counts = np.array([[1.0, 2.0, 3.0], [4.0, 6.0, 9.0]])
        rse = _make_tpm_rse(
            counts,
            [100.0, 100.0, 100.0],
            widths=[1000, 3000],
            row_names=["G", "G"],
            col_names=["S", "S", "T"],
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        assert list(result.columns) == ["S", "S", "T"]
        _assert_tpm_equal(
            result, _reference_tpm(counts, [1000, 3000], [100.0] * 3)
        )

    def test_missing_names_use_range_indexes(self) -> None:
        rse = _make_tpm_rse(
            _RAW.copy(), [100.0, 200.0, 100.0], widths=[1000, 2000]
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        assert isinstance(result.index, pd.RangeIndex)
        assert isinstance(result.columns, pd.RangeIndex)

    def test_legacy_counts_assay_is_used_with_a_warning(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        rse = _make_tpm_rse(
            _RAW.copy(),
            [100.0, 200.0, 100.0],
            widths=[1000, 2000],
            assay_name="counts",
        )
        with _patch_tpm_rse_class(), caplog.at_level(logging.WARNING):
            result = compute_tpm(rse)
        assert "falling back to legacy assay 'counts'" in caplog.text
        _assert_tpm_equal(
            result, _reference_tpm(_RAW, [1000, 2000], [100.0, 200.0, 100.0])
        )

    @pytest.mark.parametrize(
        "column",
        [
            "recount_qc__star.average_mapped_length",
            "RECOUNT_QC.STAR.AVERAGE_MAPPED_LENGTH",
        ],
    )
    def test_read_length_column_spellings_resolve(self, column: str) -> None:
        rse = _make_tpm_rse(
            _RAW.copy(),
            [100.0, 200.0, 100.0],
            widths=[1000, 2000],
            read_length_column=column,
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(
            result, _reference_tpm(_RAW, [1000, 2000], [100.0, 200.0, 100.0])
        )

    def test_list_assay_is_accepted(self) -> None:
        rse = _make_tpm_rse(
            _RAW.tolist(), [100.0, 200.0, 100.0], widths=[1000, 2000]
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(
            result, _reference_tpm(_RAW, [1000, 2000], [100.0, 200.0, 100.0])
        )

    @pytest.mark.parametrize("dtype", [np.float64, np.float32, np.int64])
    def test_input_assay_is_not_modified(self, dtype: Any) -> None:
        counts = _RAW.astype(dtype)
        original = counts.copy()
        rse = _make_tpm_rse(counts, [100.0, 200.0, 100.0], widths=[1000, 2000])
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        np.testing.assert_array_equal(counts, original)
        result.iloc[0, 0] = -1.0
        np.testing.assert_array_equal(counts, original)

    def test_read_only_assay_is_accepted(self) -> None:
        counts = _RAW.copy()
        counts.setflags(write=False)
        rse = _make_tpm_rse(counts, [100.0, 200.0, 100.0], widths=[1000, 2000])
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(
            result, _reference_tpm(_RAW, [1000, 2000], [100.0, 200.0, 100.0])
        )

    def test_no_features(self) -> None:
        rse = _make_tpm_rse(
            np.zeros((0, 3)),
            [100.0, 200.0, 100.0],
            widths=[],
            col_names=["S1", "S2", "S3"],
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        assert result.shape == (0, 3)
        assert list(result.columns) == ["S1", "S2", "S3"]

    def test_no_samples(self) -> None:
        rse = _make_tpm_rse(
            np.zeros((2, 0)), [], widths=[1000, 2000], row_names=["G1", "G2"]
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        assert result.shape == (2, 0)
        assert list(result.index) == ["G1", "G2"]


class TestComputeTpmFeatureLengths:
    """compute_tpm prefers annotated lengths over genomic spans."""

    def test_bp_length_takes_precedence_over_widths(self) -> None:
        rse = _make_tpm_rse(
            _RAW.copy(),
            [100.0, 200.0, 100.0],
            widths=[1000, 1000],
            bp_length=[250.0, 4000.0],
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(
            result, _reference_tpm(_RAW, [250.0, 4000.0], [100.0, 200.0, 100.0])
        )

    def test_bp_length_alone_is_enough(self) -> None:
        """Annotated lengths need no rowRanges widths."""
        rse = _make_tpm_rse(
            _RAW.copy(), [100.0, 200.0, 100.0], bp_length=["250", 4000]
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(
            result, _reference_tpm(_RAW, [250.0, 4000.0], [100.0, 200.0, 100.0])
        )

    def test_row_data_without_bp_length_falls_back_to_widths(self) -> None:
        rse = _make_tpm_rse(
            _RAW.copy(),
            [100.0, 200.0, 100.0],
            widths=[1000, 3000],
            row_data=_MockRowData({"gene_name": ["A", "B"]}),
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(
            result, _reference_tpm(_RAW, [1000, 3000], [100.0, 200.0, 100.0])
        )

    def test_absent_row_data_falls_back_to_widths(self) -> None:
        rse = _make_tpm_rse(
            _RAW.copy(),
            [100.0, 200.0, 100.0],
            widths=[1000, 3000],
            row_data=None,
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(
            result, _reference_tpm(_RAW, [1000, 3000], [100.0, 200.0, 100.0])
        )

    def test_widths_from_a_numpy_array(self) -> None:
        rse = _make_tpm_rse(
            _RAW.copy(),
            [100.0, 200.0, 100.0],
            widths=np.array([1000, 3000], dtype=np.int32),
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(
            result, _reference_tpm(_RAW, [1000, 3000], [100.0, 200.0, 100.0])
        )


class TestComputeTpmMissingValues:
    """Missing or degenerate inputs propagate as NaN, never as warnings."""

    def test_missing_range_width_blanks_the_feature_only(self) -> None:
        counts = np.array([[5.0, 1.0], [10.0, 3.0]])
        rse = _make_tpm_rse(counts, [100.0, 100.0], widths=[None, 1000])
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(result, np.array([[np.nan, np.nan], [1e6, 1e6]]))

    @pytest.mark.parametrize("bad_length", [np.nan, 0.0, np.inf, -np.inf])
    def test_unusable_read_length_blanks_only_that_sample(
        self, bad_length: float
    ) -> None:
        counts = _RAW.copy()
        rse = _make_tpm_rse(
            counts, [100.0, bad_length, 100.0], widths=[1000, 2000]
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        values = result.to_numpy()
        assert np.isnan(values[:, 1]).all()
        expected = _reference_tpm(counts, [1000, 2000], [100.0, 100.0, 100.0])
        np.testing.assert_allclose(values[:, [0, 2]], expected[:, [0, 2]])

    def test_missing_read_length_as_none_blanks_the_sample(self) -> None:
        col_data = pd.DataFrame({_READ_LENGTH_COLUMN: [100.0, None, 100.0]})
        rse = _MockRSE(
            assays={"raw_counts": _RAW.copy()},
            col_data_df=col_data.astype(object),
            row_ranges_widths=[1000, 2000],
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        assert np.isnan(result.to_numpy()[:, 1]).all()
        assert not np.isnan(result.to_numpy()[:, [0, 2]]).any()

    def test_all_zero_sample_is_nan(self) -> None:
        counts = np.array([[0.0, 10.0], [0.0, 30.0]])
        rse = _make_tpm_rse(counts, [100.0, 100.0], widths=[1000, 1000])
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(
            result, np.array([[np.nan, 250_000.0], [np.nan, 750_000.0]])
        )

    @pytest.mark.parametrize("missing", [np.nan, None])
    def test_missing_feature_length_blanks_the_feature_only(
        self, missing: Any
    ) -> None:
        counts = np.array([[5.0, 0.0], [10.0, 20.0], [30.0, 0.0]])
        rse = _make_tpm_rse(
            counts, [100.0, 100.0], bp_length=[missing, 1000.0, 3000.0]
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        # The remaining features still share a full million per sample.
        _assert_tpm_equal(
            result,
            np.array([[np.nan, np.nan], [5e5, 1e6], [5e5, 0.0]]),
        )

    def test_zero_feature_length(self) -> None:
        """Division by a zero length follows IEEE arithmetic."""
        counts = np.array([[100.0, 0.0], [400.0, 500.0], [0.0, 0.0]])
        rse = _make_tpm_rse(
            counts, [100.0, 100.0], bp_length=[0.0, 1000.0, 500.0]
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        # Sample 1: an infinite RPK makes the scale infinite, so every finite
        # RPK becomes 0 and the infinite one becomes NaN.  Sample 2: 0 / 0 is
        # NaN and is skipped when summing.
        _assert_tpm_equal(
            result,
            np.array([[np.nan, np.nan], [0.0, 1e6], [0.0, 0.0]]),
        )

    def test_infinite_feature_length_contributes_zero(self) -> None:
        counts = np.array([[100.0], [400.0]])
        rse = _make_tpm_rse(counts, [100.0], bp_length=[np.inf, 1000.0])
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(result, np.array([[0.0], [1e6]]))

    def test_missing_count_is_skipped_by_the_sample_sum(self) -> None:
        counts = np.array([[np.nan, 1.0], [10.0, 1.0], [30.0, 2.0]])
        rse = _make_tpm_rse(counts, [100.0, 100.0], widths=[1000] * 3)
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(
            result,
            np.array(
                [
                    [np.nan, 250_000.0],
                    [250_000.0, 250_000.0],
                    [750_000.0, 500_000.0],
                ]
            ),
        )

    def test_infinite_count_is_nan_and_zeroes_its_sample(self) -> None:
        counts = np.array([[np.inf, 1.0], [10.0, 3.0]])
        rse = _make_tpm_rse(counts, [100.0, 100.0], widths=[1000, 1000])
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(
            result, np.array([[np.nan, 250_000.0], [0.0, 750_000.0]])
        )

    def test_cancelling_negative_counts_follow_ieee_arithmetic(self) -> None:
        """Signed values summing to zero divide by a zero scale factor."""
        counts = np.array([[5.0, 1.0], [-5.0, 1.0], [0.0, 2.0]])
        rse = _make_tpm_rse(counts, [100.0, 100.0], widths=[1000] * 3)
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(
            result,
            np.array(
                [[np.inf, 250_000.0], [-np.inf, 250_000.0], [np.nan, 500_000.0]]
            ),
        )

    def test_degenerate_inputs_emit_no_warnings(self) -> None:
        counts = np.array([[0.0, 1.0, 0.0], [0.0, 3.0, np.nan]])
        rse = _make_tpm_rse(
            counts, [100.0, 0.0, np.nan], bp_length=[0.0, np.nan]
        )
        with _patch_tpm_rse_class(), warnings.catch_warnings():
            warnings.simplefilter("error")
            result = compute_tpm(rse)
        assert np.isnan(result.to_numpy()).all()


_SPARSE_FORMATS = [
    sparse.csc_matrix,
    sparse.csr_matrix,
    sparse.coo_matrix,
    sparse.csc_array,
    sparse.csr_array,
]


class TestComputeTpmSparse:
    """Sparse assays give the dense result."""

    @pytest.mark.parametrize("to_sparse", _SPARSE_FORMATS)
    @pytest.mark.parametrize("seed", range(3))
    def test_matches_dense_result(self, to_sparse: Any, seed: int) -> None:
        counts, widths, read_lengths = _random_tpm_case(seed, density=0.1)
        sparse_rse = _make_tpm_rse(
            to_sparse(counts), read_lengths, widths=widths
        )
        dense_rse = _make_tpm_rse(counts, read_lengths, widths=widths)
        with _patch_tpm_rse_class():
            sparse_result = compute_tpm(sparse_rse)
            dense_result = compute_tpm(dense_rse)
        assert all(isinstance(d, pd.SparseDtype) for d in sparse_result.dtypes)
        _assert_tpm_equal(sparse_result, dense_result.to_numpy())
        _assert_tpm_equal(
            sparse_result, _reference_tpm(counts, widths, read_lengths)
        )

    def test_degenerate_values_match_dense_result(self) -> None:
        """NaN/zero lengths, zero samples and stored zeros, all at once."""
        counts = np.array(
            [
                [0.0, 0.0, 7.0, 0.0, 0.0, 2.0],
                [3.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 5.0, 0.0, 0.0, 0.0],
                [9.0, 0.0, 0.0, 4.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 6.0, 0.0],
            ]
        )
        widths = [1000.0, np.nan, 0.0, 500.0, 2000.0]
        read_lengths = [100.0, 100.0, 100.0, np.nan, 100.0, 0.0]
        # Store explicit zeros, which must behave like implicit ones: one in
        # an all-zero sample and one in a valid feature and sample.
        rows, cols = np.nonzero(counts)
        matrix = sparse.coo_matrix(
            (
                np.append(counts[rows, cols], [0.0, 0.0]),
                (np.append(rows, [0, 3]), np.append(cols, [1, 4])),
            ),
            shape=counts.shape,
        ).tocsc()
        assert matrix.nnz == np.count_nonzero(counts) + 2
        sparse_rse = _make_tpm_rse(matrix, read_lengths, bp_length=widths)
        with _patch_tpm_rse_class():
            result = compute_tpm(sparse_rse)
        _assert_tpm_equal(result, _reference_tpm(counts, widths, read_lengths))

    def test_signed_values_match_dense_result(self) -> None:
        counts = np.array([[5.0, 0.0, 1.0], [-5.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
        widths = [1000.0, 1000.0, 1000.0]
        rse = _make_tpm_rse(
            sparse.csc_matrix(counts), [100.0] * 3, widths=widths
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(result, _reference_tpm(counts, widths, [100.0] * 3))

    def test_integer_sparse_input(self) -> None:
        counts, widths, read_lengths = _random_tpm_case(3)
        matrix = sparse.csc_matrix(counts.astype(np.int32))
        rse = _make_tpm_rse(matrix, read_lengths, widths=widths)
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(result, _reference_tpm(counts, widths, read_lengths))

    @pytest.mark.parametrize("to_sparse", _SPARSE_FORMATS)
    def test_input_matrix_is_not_modified(self, to_sparse: Any) -> None:
        counts, widths, read_lengths = _random_tpm_case(4)
        matrix = to_sparse(counts.astype(float))
        original = matrix.copy()
        rse = _make_tpm_rse(matrix, read_lengths, widths=widths)
        with _patch_tpm_rse_class():
            compute_tpm(rse)
        assert (matrix != original).nnz == 0
        np.testing.assert_array_equal(matrix.data, original.data)

    def test_names_are_preserved(self) -> None:
        rse = _make_tpm_rse(
            sparse.csc_matrix(_RAW),
            [100.0, 200.0, 100.0],
            widths=[1000, 2000],
            row_names=["GENE1", "GENE2"],
            col_names=["S1", "S2", "S3"],
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        assert list(result.index) == ["GENE1", "GENE2"]
        assert list(result.columns) == ["S1", "S2", "S3"]

    def test_empty_sparse_matrix(self) -> None:
        rse = _make_tpm_rse(
            sparse.csc_matrix((0, 2)), [100.0, 100.0], widths=[]
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        assert result.shape == (0, 2)


class _MockColumnFrame:
    """Column data exposing BiocFrame's per-column accessors."""

    def __init__(
        self, columns: dict[str, Any], flattened: pd.DataFrame
    ) -> None:
        self._columns = columns
        self._flattened = flattened
        self.to_pandas_calls = 0

    @property
    def column_names(self) -> list[str]:
        return list(self._columns)

    def get_column(self, name: str) -> Any:
        return self._columns[name]

    def to_pandas(self) -> pd.DataFrame:
        self.to_pandas_calls += 1
        return self._flattened.copy()


class _FrameWithoutColumnNames(_MockBiocFrame):
    """Column data with ``get_column`` but no ``column_names``."""

    def get_column(self, name: str) -> Any:
        raise AssertionError("get_column must not be used")


class _FrameWithoutGetColumn(_MockBiocFrame):
    """Column data with ``column_names`` but no ``get_column``."""

    @property
    def column_names(self) -> list[str]:
        return list(self._df.columns)


class _MockRSEWithColumnFrame(_MockRSE):
    def __init__(self, column_frame: _MockColumnFrame, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self._bioc_col_data = column_frame  # type: ignore[assignment]


class TestComputeTpmSampleMetadata:
    """compute_tpm reads only the read-length column when it can."""

    def _rse(self, column_frame: _MockColumnFrame) -> _MockRSE:
        return _MockRSEWithColumnFrame(
            column_frame,
            assays={"raw_counts": _RAW.copy()},
            row_ranges_widths=[1000, 2000],
        )

    @pytest.mark.parametrize(
        "column",
        [
            _READ_LENGTH_COLUMN,
            "recount_qc__star.average_mapped_length",
            "Recount_QC.Star.Average_Mapped_Length",
        ],
    )
    def test_reads_the_column_without_converting_the_frame(
        self, column: str
    ) -> None:
        frame = _MockColumnFrame(
            {"external_id": ["S1", "S2", "S3"], column: [100.0, None, 50.0]},
            flattened=pd.DataFrame(),
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(self._rse(frame))
        assert frame.to_pandas_calls == 0
        values = result.to_numpy()
        assert np.isnan(values[:, 1]).all()
        np.testing.assert_allclose(
            values[:, [0, 2]],
            _reference_tpm(_RAW, [1000, 2000], [100.0, 1.0, 50.0])[:, [0, 2]],
            rtol=1e-12,
        )

    def test_values_are_validated_like_the_pandas_path(self) -> None:
        frame = _MockColumnFrame(
            {_READ_LENGTH_COLUMN: [100.0, "abc", 50.0]},
            flattened=pd.DataFrame(),
        )
        with _patch_tpm_rse_class():
            with pytest.raises(
                ValueError, match=r"must be numeric; found str value 'abc'"
            ):
                compute_tpm(self._rse(frame))

    def test_flattened_name_falls_back_to_pandas(self) -> None:
        """A name produced by flattening a nested frame still resolves."""
        frame = _MockColumnFrame(
            {"recount_qc": object()},
            flattened=pd.DataFrame({_READ_LENGTH_COLUMN: [100.0, 200.0, 50.0]}),
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(self._rse(frame))
        assert frame.to_pandas_calls == 1
        _assert_tpm_equal(
            result, _reference_tpm(_RAW, [1000, 2000], [100.0, 200.0, 50.0])
        )

    def test_nested_frame_column_falls_back_to_pandas(self) -> None:
        nested = _MockColumnFrame({"x": [1, 2, 3]}, flattened=pd.DataFrame())
        frame = _MockColumnFrame(
            {_READ_LENGTH_COLUMN: nested},
            flattened=pd.DataFrame(
                {"recount_qc__star.average_mapped_length": [10.0, 20.0, 5.0]}
            ),
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(self._rse(frame))
        assert frame.to_pandas_calls == 1
        _assert_tpm_equal(
            result, _reference_tpm(_RAW, [1000, 2000], [10.0, 20.0, 5.0])
        )

    @pytest.mark.parametrize(
        "frame_cls", [_FrameWithoutColumnNames, _FrameWithoutGetColumn]
    )
    def test_partial_frame_interface_uses_pandas(
        self, frame_cls: type[_MockBiocFrame]
    ) -> None:
        frame = frame_cls(
            pd.DataFrame({_READ_LENGTH_COLUMN: [100.0, 200.0, 50.0]})
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(self._rse(frame))  # type: ignore[arg-type]
        _assert_tpm_equal(
            result, _reference_tpm(_RAW, [1000, 2000], [100.0, 200.0, 50.0])
        )

    def test_missing_column_raises_value_error(self) -> None:
        frame = _MockColumnFrame(
            {"external_id": ["S1", "S2", "S3"]},
            flattened=pd.DataFrame({"external_id": ["S1", "S2", "S3"]}),
        )
        with _patch_tpm_rse_class():
            with pytest.raises(
                ValueError,
                match=(
                    r"^Required metadata column not found in col_data: "
                    r"'recount_qc\.star\.average_mapped_length'\.$"
                ),
            ):
                compute_tpm(self._rse(frame))


def _degenerate_tpm_case(
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Random sparse-ish counts with missing and zero lengths mixed in."""
    rng = np.random.default_rng(seed)
    counts, widths, read_lengths = _random_tpm_case(
        seed, n_features=30, n_samples=9, density=0.2
    )
    counts = counts.astype(float)
    widths[rng.choice(30, size=3, replace=False)] = rng.choice(
        [np.nan, 0.0, np.inf], size=3
    )
    read_lengths[rng.choice(9, size=2, replace=False)] = rng.choice(
        [np.nan, 0.0, np.inf], size=2
    )
    counts[:, rng.integers(9)] = 0.0
    return counts, widths, read_lengths


class TestComputeTpmDenseLayout:
    """Dense results use pandas' own column layout."""

    @pytest.mark.parametrize("order", ["C", "F"])
    def test_each_sample_is_contiguous(self, order: str) -> None:
        counts, widths, read_lengths = _random_tpm_case(12, 50, 6)
        rse = _make_tpm_rse(
            np.asarray(counts, order=order), read_lengths, widths=widths
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        for column in range(result.shape[1]):
            assert result.iloc[:, column].to_numpy().flags.c_contiguous
        _assert_tpm_equal(result, _reference_tpm(counts, widths, read_lengths))

    @pytest.mark.parametrize(
        "dtype", [np.bool_, np.uint16, np.int32, np.float16, np.float32]
    )
    def test_other_real_dtypes_are_computed_in_float64(
        self, dtype: Any
    ) -> None:
        counts = np.array([[1, 0, 3], [2, 1, 1]], dtype=dtype)
        rse = _make_tpm_rse(counts, [100.0] * 3, widths=[1000, 3000])
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        assert all(d == np.float64 for d in result.dtypes)
        _assert_tpm_equal(
            result,
            _reference_tpm(
                counts.astype(np.float64), [1000, 3000], [100.0] * 3
            ),
        )

    def test_object_assay_is_converted_to_float(self) -> None:
        counts = np.array([[1, None], ["2", 3.0]], dtype=object)
        rse = _make_tpm_rse(counts, [100.0] * 2, widths=[1000, 1000])
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(result, np.array([[1e6 / 3, np.nan], [2e6 / 3, 1e6]]))


class TestComputeTpmSparseStorage:
    """Sparse results stay sparse and match the dense computation."""

    @pytest.mark.parametrize("to_sparse", _SPARSE_FORMATS)
    @pytest.mark.parametrize("seed", range(6))
    def test_degenerate_random_input_matches_reference(
        self, to_sparse: Any, seed: int
    ) -> None:
        counts, widths, read_lengths = _degenerate_tpm_case(seed)
        rse = _make_tpm_rse(to_sparse(counts), read_lengths, bp_length=widths)
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(result, _reference_tpm(counts, widths, read_lengths))

    def test_only_nonzero_counts_are_stored(self) -> None:
        counts, widths, read_lengths = _random_tpm_case(
            21, n_features=2000, n_samples=20, density=0.01
        )
        matrix = sparse.csc_matrix(counts)
        rse = _make_tpm_rse(matrix, read_lengths, widths=widths)
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        assert all(
            dtype == pd.SparseDtype(np.float64, 0.0) for dtype in result.dtypes
        )
        stored = sum(result[c].array.npoints for c in result.columns)
        assert stored == matrix.nnz
        _assert_tpm_equal(result, _reference_tpm(counts, widths, read_lengths))

    def test_nan_lines_are_stored_once_each(self) -> None:
        """Flagged features and samples store each entry exactly once."""
        counts = np.array(
            [
                [1.0, 0.0, 0.0, 2.0],
                [0.0, 3.0, 0.0, 0.0],
                [4.0, 0.0, 0.0, 5.0],
                [0.0, 2.0, 0.0, 6.0],
            ]
        )
        # Feature 1 has no length; sample 2 has no coverage; sample 3 has no
        # read length.
        rse = _make_tpm_rse(
            sparse.csc_matrix(counts),
            [100.0, 100.0, 100.0, np.nan],
            bp_length=[1000.0, np.nan, 1000.0, 1000.0],
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        # Samples 2 and 3 store every feature; samples 0 and 1 store their
        # counts plus feature 1, whose count in sample 1 is already stored.
        stored = [result[c].array.npoints for c in result.columns]
        assert stored == [3, 2, 4, 4]
        _assert_tpm_equal(
            result,
            np.array(
                [
                    [2e5, 0.0, np.nan, np.nan],
                    [np.nan, np.nan, np.nan, np.nan],
                    [8e5, 0.0, np.nan, np.nan],
                    [0.0, 1e6, np.nan, np.nan],
                ]
            ),
        )

    def test_several_flagged_features_in_normal_samples(self) -> None:
        counts = np.array(
            [
                [1.0, 0.0, 2.0],
                [0.0, 7.0, 0.0],
                [3.0, 0.0, 0.0],
                [0.0, 0.0, 6.0],
                [5.0, 4.0, 2.0],
            ]
        )
        widths = [np.nan, 1000.0, 0.0, np.nan, 1000.0]
        rse = _make_tpm_rse(
            sparse.csc_matrix(counts), [100.0] * 3, bp_length=widths
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        # Feature 2 has no length but coverage in sample 0, so that sample's
        # scale is 0 and its finite values become zeros that are not stored.
        assert [result[c].array.npoints for c in result.columns] == [3, 5, 4]
        _assert_tpm_equal(result, _reference_tpm(counts, widths, [100.0] * 3))

    def test_duplicate_entries_are_summed(self) -> None:
        matrix = sparse.coo_matrix(
            ([100.0, 200.0, 300.0, 50.0], ([0, 0, 1, 1], [0, 0, 0, 1])),
            shape=(2, 2),
        )
        rse = _make_tpm_rse(matrix, [100.0, 100.0], widths=[1000, 1000])
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(result, np.array([[5e5, 0.0], [5e5, 1e6]]))

    def test_unsorted_and_duplicate_csc_input(self) -> None:
        matrix = sparse.csc_matrix(
            (
                np.array([30.0, 10.0, 10.0, 7.0]),
                np.array([2, 0, 0, 1]),
                np.array([0, 3, 4]),
            ),
            shape=(3, 2),
        )
        assert not matrix.has_canonical_format
        original = matrix.copy()
        rse = _make_tpm_rse(matrix, [100.0, 100.0], widths=[1000] * 3)
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        _assert_tpm_equal(
            result, np.array([[4e5, 0.0], [0.0, 1e6], [6e5, 0.0]])
        )
        np.testing.assert_array_equal(matrix.indices, original.indices)

    def test_all_samples_unusable(self) -> None:
        rse = _make_tpm_rse(
            sparse.csc_matrix(_RAW), [np.nan] * 3, widths=[1000, 2000]
        )
        with _patch_tpm_rse_class():
            result = compute_tpm(rse)
        assert np.isnan(result.to_numpy(dtype=float)).all()
        assert [result[c].array.npoints for c in result.columns] == [2, 2, 2]


class TestComputeTpmMemory:
    """compute_tpm allocates no more than one output-sized array."""

    _SHAPE = (4000, 50)
    # NumPy's ufunc buffers and other fixed-size scratch space; a second
    # matrix-sized array would be 1.6 MB.
    _OVERHEAD = 512 * 1024

    @staticmethod
    def _peak_bytes(rse: _MockRSE) -> int:
        tracemalloc.start()
        try:
            with _patch_tpm_rse_class():
                result = compute_tpm(rse)
            peak = tracemalloc.get_traced_memory()[1]
        finally:
            tracemalloc.stop()
        assert result.shape == TestComputeTpmMemory._SHAPE
        return peak

    @pytest.mark.parametrize("dtype", [np.int64, np.int32, np.float64])
    def test_dense_peak_is_one_output_array(self, dtype: Any) -> None:
        counts, widths, read_lengths = _random_tpm_case(
            1, *self._SHAPE, density=0.9
        )
        rse = _make_tpm_rse(counts.astype(dtype), read_lengths, widths=widths)
        output_bytes = counts.size * 8
        assert self._peak_bytes(rse) < output_bytes + self._OVERHEAD

    def test_dense_missing_lengths_add_at_most_a_mask(self) -> None:
        counts, widths, read_lengths = _random_tpm_case(
            2, *self._SHAPE, density=0.9
        )
        widths[0] = np.nan
        rse = _make_tpm_rse(counts, read_lengths, widths=widths)
        output_bytes = counts.size * 8
        # One boolean mask, an eighth of the float output, is allowed.
        mask_bytes = counts.size
        assert (
            self._peak_bytes(rse) < output_bytes + mask_bytes + self._OVERHEAD
        )

    def test_sparse_peak_scales_with_stored_values(self) -> None:
        counts, widths, read_lengths = _random_tpm_case(
            3, *self._SHAPE, density=0.01
        )
        matrix = sparse.csc_matrix(counts)
        rse = _make_tpm_rse(matrix, read_lengths, widths=widths)
        dense_bytes = counts.size * 8
        # Dense float output would need 1.6 MB; stored values need ~16 kB.
        assert self._peak_bytes(rse) < dense_bytes / 10


@pytest.mark.requires_biocpy
class TestComputeTpmWithBiocPy:
    """compute_tpm on real BiocPy experiments."""

    @staticmethod
    def _experiment(
        counts: Any,
        read_lengths: Sequence[Any],
        *,
        widths: Sequence[int],
        bp_length: Sequence[Any] | None = None,
        extra_column_data: dict[str, Any] | None = None,
    ) -> Any:
        biocframe = pytest.importorskip("biocframe")
        genomicranges = pytest.importorskip("genomicranges")
        iranges = pytest.importorskip("iranges")
        summarizedexperiment = pytest.importorskip("summarizedexperiment")

        n_features, n_samples = counts.shape
        features = [f"ENSG{i:011d}" for i in range(n_features)]
        samples = [f"SRR{j:06d}" for j in range(n_samples)]
        row_data: dict[str, Any] = {"gene_id": features}
        if bp_length is not None:
            row_data["bp_length"] = list(bp_length)
        column_data: dict[str, Any] = {
            "external_id": samples,
            _READ_LENGTH_COLUMN: list(read_lengths),
        }
        column_data.update(extra_column_data or {})
        return summarizedexperiment.RangedSummarizedExperiment(
            assays={"raw_counts": counts},
            row_ranges=genomicranges.GenomicRanges(
                seqnames=["chr1"] * n_features,
                ranges=iranges.IRanges(
                    start=[1] * n_features, width=list(widths)
                ),
                names=features,
            ),
            row_data=biocframe.BiocFrame(row_data, row_names=features),
            column_data=biocframe.BiocFrame(column_data, row_names=samples),
        )

    def test_uses_bp_length_when_annotated(self) -> None:
        counts, widths, read_lengths = _random_tpm_case(5)
        bp_length = widths[::-1].copy()
        rse = self._experiment(
            counts,
            read_lengths,
            widths=widths.astype(int),
            bp_length=bp_length,
        )
        result = compute_tpm(rse)
        assert list(result.index) == list(rse.get_row_names())
        assert list(result.columns) == list(rse.get_column_names())
        _assert_tpm_equal(
            result, _reference_tpm(counts, bp_length, read_lengths)
        )

    def test_uses_range_widths_without_bp_length(self) -> None:
        counts, widths, read_lengths = _random_tpm_case(6)
        rse = self._experiment(counts, read_lengths, widths=widths.astype(int))
        _assert_tpm_equal(
            compute_tpm(rse), _reference_tpm(counts, widths, read_lengths)
        )

    def test_sparse_assay(self) -> None:
        counts, widths, read_lengths = _random_tpm_case(7, density=0.05)
        rse = self._experiment(
            sparse.csc_matrix(counts),
            read_lengths,
            widths=widths.astype(int),
        )
        result = compute_tpm(rse)
        assert all(isinstance(d, pd.SparseDtype) for d in result.dtypes)
        _assert_tpm_equal(result, _reference_tpm(counts, widths, read_lengths))

    def test_missing_read_length_in_column_data(self) -> None:
        counts, widths, read_lengths = _random_tpm_case(8)
        lengths: list[Any] = list(read_lengths)
        lengths[2] = None
        rse = self._experiment(counts, lengths, widths=widths.astype(int))
        result = compute_tpm(rse).to_numpy()
        assert np.isnan(result[:, 2]).all()
        assert not np.isnan(np.delete(result, 2, axis=1)).any()

    def test_namespaced_read_length_column(self) -> None:
        counts, widths, read_lengths = _random_tpm_case(9)
        rse = self._experiment(counts, read_lengths, widths=widths.astype(int))
        renamed = rse.set_column_data(
            rse.get_column_data().set_column_names(
                ["external_id", "recount_qc__star.average_mapped_length"]
            )
        )
        _assert_tpm_equal(
            compute_tpm(renamed), _reference_tpm(counts, widths, read_lengths)
        )

    def test_plain_summarized_experiment_is_rejected(self) -> None:
        biocframe = pytest.importorskip("biocframe")
        summarizedexperiment = pytest.importorskip("summarizedexperiment")
        experiment = summarizedexperiment.SummarizedExperiment(
            assays={"raw_counts": _RAW.copy()},
            column_data=biocframe.BiocFrame(
                {_READ_LENGTH_COLUMN: [100.0, 200.0, 100.0]}
            ),
        )
        with pytest.raises(TypeError, match="compute_tpm requires"):
            compute_tpm(experiment)


class TestIsPairedEnd:
    """Tests for is_paired_end."""

    def _meta(
        self,
        mapped: list,
        avg_len: list,
        external_ids: list | None = None,
    ) -> pd.DataFrame:
        ids = external_ids or [f"S{i+1}" for i in range(len(mapped))]
        return pd.DataFrame(
            {
                "external_id": ids,
                "recount_qc.star.average_mapped_length": mapped,
                "recount_seq_qc.avg_len": avg_len,
            }
        )

    def test_single_end_ratio_one(self) -> None:
        meta = self._meta([100.0], [100.0])
        result = is_paired_end(meta)
        assert result.dtype == "boolean"
        assert result.iloc[0] == False

    def test_paired_end_ratio_two(self) -> None:
        meta = self._meta([200.0], [100.0])
        result = is_paired_end(meta)
        assert result.iloc[0] == True

    def test_mixed_se_and_pe(self) -> None:
        meta = self._meta([100.0, 200.0], [100.0, 100.0])
        result = is_paired_end(meta)
        assert result.iloc[0] == False
        assert result.iloc[1] == True

    def test_invalid_ratio_becomes_na_with_warning(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        meta = self._meta([300.0], [100.0])
        with caplog.at_level(logging.WARNING):
            result = is_paired_end(meta)
        assert pd.isna(result.iloc[0])
        assert "paired" in caplog.text.lower()

    def test_mixed_with_invalid_ratio(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        meta = self._meta([100.0, 300.0, 200.0], [100.0, 100.0, 100.0])
        with caplog.at_level(logging.WARNING):
            result = is_paired_end(meta)
        assert result.iloc[0] == False
        assert pd.isna(result.iloc[1])
        assert result.iloc[2] == True

    def test_indexed_by_external_id(self) -> None:
        meta = self._meta(
            [100.0, 200.0], [100.0, 100.0], external_ids=["ERR001", "ERR002"]
        )
        result = is_paired_end(meta)
        assert list(result.index) == ["ERR001", "ERR002"]

    def test_nan_in_mapped_length_propagates(self) -> None:
        meta = self._meta([float("nan"), 200.0], [100.0, 100.0])
        result = is_paired_end(meta)
        assert pd.isna(result.iloc[0])
        assert result.iloc[1] == True

    def test_se_like_object_accepted(self) -> None:
        """Objects with col_data.to_pandas() are coerced automatically."""
        meta = self._meta([100.0], [100.0])
        se_like = mock.MagicMock()
        se_like.col_data.to_pandas.return_value = meta
        result = is_paired_end(se_like)
        assert result.iloc[0] == False

    def test_missing_avg_mapped_length_column_raises(self) -> None:
        meta = pd.DataFrame(
            {"external_id": ["S1"], "recount_seq_qc.avg_len": [100.0]}
        )
        with pytest.raises(ValueError, match="not found"):
            is_paired_end(meta)

    def test_missing_avg_read_length_column_raises(self) -> None:
        meta = pd.DataFrame(
            {
                "external_id": ["S1"],
                "recount_qc.star.average_mapped_length": [100.0],
            }
        )
        with pytest.raises(ValueError, match="not found"):
            is_paired_end(meta)

    def test_custom_column_names(self) -> None:
        meta = pd.DataFrame(
            {
                "external_id": ["S1"],
                "my.mapped_len": [100.0],
                "my.read_len": [100.0],
            }
        )
        result = is_paired_end(
            meta,
            avg_mapped_read_length_column="my.mapped_len",
            avg_read_length_column="my.read_len",
        )
        assert result.iloc[0] == False


class TestAucScalingWithoutSequencingQc:
    """AUC scaling must not depend on the ``recount_seq_qc`` table.

    Paired-end status cancels out of the AUC formula. R reaches that
    conclusion through lazy evaluation -- ``compute_scale_factors()``
    defaults ``paired_end`` to an ``is_paired_end()`` call that the AUC
    branch never forces -- so AUC scaling keeps working on projects whose
    ``recount_seq_qc`` table is empty or was never published. recount3
    documents exactly such a project in its own suite (SRP103067).

    Executed against the R package's own ``R/transform_counts.R``, R gives
    ``compute_scale_factors(meta, by="auc")`` -> ``0.04, 0.02`` on the
    metadata below, and errors only for ``by="mapped_reads"``.
    """

    @staticmethod
    def _meta_without_avg_len() -> pd.DataFrame:
        return pd.DataFrame(
            {
                "external_id": ["s1", "s2"],
                "recount_qc.bc_auc.all_reads_all_bases": [1e9, 2e9],
                "recount_qc.star.average_mapped_length": [100.0, 100.0],
                "recount_qc.star.all_mapped_reads": [1e7, 2e7],
            }
        )

    def test_auc_does_not_need_the_average_read_length(self) -> None:
        result = compute_scale_factors(
            self._meta_without_avg_len(), by="auc", target_read_count=4e7
        )
        assert result["s1"] == pytest.approx(0.04)
        assert result["s2"] == pytest.approx(0.02)

    def test_mapped_reads_still_needs_the_average_read_length(self) -> None:
        with pytest.raises(ValueError, match="recount_seq_qc.avg_len"):
            compute_scale_factors(
                self._meta_without_avg_len(), by="mapped_reads"
            )

    def test_auc_ignores_a_supplied_paired_end_status(self) -> None:
        """It cancels out, so even an unusable value cannot matter."""
        result = compute_scale_factors(
            self._meta_without_avg_len(),
            by="auc",
            paired_end_status=[True, False, True],
        )
        assert list(result.index) == ["s1", "s2"]

    def test_auc_matches_the_value_computed_with_the_table_present(
        self,
    ) -> None:
        without = compute_scale_factors(self._meta_without_avg_len(), by="auc")
        with_table = compute_scale_factors(
            self._meta_without_avg_len().assign(
                **{"recount_seq_qc.avg_len": [100.0, 100.0]}
            ),
            by="auc",
        )
        pd.testing.assert_series_equal(without, with_table)

    def test_paired_end_is_not_resolved_for_auc(self) -> None:
        with mock.patch(
            "recount3.se.is_paired_end",
            side_effect=AssertionError("must not be called for by='auc'"),
        ) as never:
            compute_scale_factors(self._meta_without_avg_len(), by="auc")
        never.assert_not_called()


class TestComputeScaleFactors:
    """Tests for compute_scale_factors."""

    def test_invalid_by_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="'auc' or 'mapped_reads'"):
            compute_scale_factors(_META.copy(), by="rpkm")

    def test_target_read_count_is_bool_raises_type_error(self) -> None:
        with pytest.raises(TypeError, match="target_read_count"):
            compute_scale_factors(_META.copy(), target_read_count=True)

    def test_target_read_length_is_bool_raises_type_error(self) -> None:
        with pytest.raises(TypeError, match="target_read_length_bp"):
            compute_scale_factors(_META.copy(), target_read_length_bp=True)

    def test_target_read_count_is_string_raises_type_error(self) -> None:
        with pytest.raises(TypeError, match="target_read_count"):
            compute_scale_factors(_META.copy(), target_read_count="40000000")

    def test_target_read_length_is_string_raises_type_error(self) -> None:
        with pytest.raises(TypeError, match="target_read_length_bp"):
            compute_scale_factors(_META.copy(), target_read_length_bp="100")

    def test_missing_auc_column_raises_value_error(self) -> None:
        meta = pd.DataFrame(
            {
                "external_id": ["S1"],
                "recount_qc.star.average_mapped_length": [100.0],
                "recount_seq_qc.avg_len": [100.0],
                "recount_qc.star.all_mapped_reads": [1e6],
                # auc column missing
            }
        )
        with pytest.raises(ValueError):
            compute_scale_factors(meta, by="auc")

    def test_auc_scale_factors_correct(self) -> None:
        result = compute_scale_factors(
            _META.copy(), by="auc", target_read_count=4e7
        )
        assert result.name == "scale_factor"
        assert result["S1"] == pytest.approx(4e7 / 4e9)
        assert result["S2"] == pytest.approx(4e7 / 8e9)
        assert result["S3"] == pytest.approx(4e7 / 4e9)

    def test_auc_scale_indexed_by_external_id(self) -> None:
        result = compute_scale_factors(_META.copy(), by="auc")
        assert list(result.index) == ["S1", "S2", "S3"]

    def test_mapped_reads_with_provided_paired_end_status(self) -> None:
        pe = [False, False, False]
        result = compute_scale_factors(
            _META.copy(),
            by="mapped_reads",
            target_read_count=4e7,
            target_read_length_bp=100.0,
            paired_end_status=pe,
        )
        assert result.name == "scale_factor"
        expected_s1 = 4e7 * 100.0 * 1.0 / (1e6 * (100.0**2))
        assert result.iloc[0] == pytest.approx(expected_s1)

    def test_mapped_reads_inferred_paired_end(self) -> None:
        """S2 has ratio=2 (PE), S1/S3 have ratio=1 (SE)."""
        result = compute_scale_factors(
            _META.copy(),
            by="mapped_reads",
            target_read_count=4e7,
            target_read_length_bp=100.0,
        )
        s1 = 4e7 * 100.0 * 1.0 / (1e6 * (100.0**2))
        s2 = 4e7 * 100.0 * 2.0 / (2e6 * (200.0**2))
        assert result.iloc[0] == pytest.approx(s1)
        assert result.iloc[1] == pytest.approx(s2)

    def test_mapped_reads_unknown_pe_status_gives_nan(self) -> None:
        """Samples with unknown PE status -> NaN scale factor."""
        pe: list = [pd.NA, False, False]
        result = compute_scale_factors(
            _META.copy(),
            by="mapped_reads",
            paired_end_status=pe,
        )
        assert pd.isna(result.iloc[0])
        assert not pd.isna(result.iloc[1])

    def test_result_is_series_named_scale_factor(self) -> None:
        result = compute_scale_factors(_META.copy(), by="auc")
        assert isinstance(result, pd.Series)
        assert result.name == "scale_factor"

    def test_se_like_input_accepted(self) -> None:
        se_like = mock.MagicMock()
        se_like.col_data.to_pandas.return_value = _META.copy()
        result = compute_scale_factors(se_like, by="auc")
        assert len(result) == 3


class TestTransformCounts:
    """Tests for transform_counts."""

    def _patch_rse_class(self) -> mock._patch:
        return mock.patch.object(
            _utils,
            "get_ranged_summarizedexperiment_class",
            return_value=_MockRSE,
        )

    def test_non_rse_raises_type_error(self) -> None:
        with self._patch_rse_class():
            with pytest.raises(TypeError, match="RangedSummarizedExperiment"):
                transform_counts("not_rse")

    def test_round_not_bool_raises_type_error(self) -> None:
        rse = _make_rse()
        with self._patch_rse_class():
            with pytest.raises(TypeError, match="bool"):
                transform_counts(rse, round_to_integers=1)  # type: ignore[arg-type]

    def test_1d_assay_raises_value_error(self) -> None:
        rse = _make_rse(counts=np.array([1.0, 2.0, 3.0]))
        with self._patch_rse_class():
            with pytest.raises(ValueError, match="2D"):
                transform_counts(rse, by="auc")

    def test_shape_mismatch_raises_value_error(self) -> None:
        rse = _make_rse(counts=np.ones((2, 4)))
        with self._patch_rse_class():
            with pytest.raises(ValueError, match="Mismatch"):
                transform_counts(rse, by="auc")

    @pytest.mark.parametrize(
        "matrix", _duplicate_entry_matrices(), ids=_DUPLICATE_IDS
    )
    def test_sparse_duplicate_entries_are_summed(self, matrix: Any) -> None:
        meta = _META.iloc[:2].copy()
        sparse_rse = _MockRSE(assays={"raw_counts": matrix}, col_data_df=meta)
        dense_rse = _MockRSE(
            assays={"raw_counts": _DUPLICATE_DENSE.copy()}, col_data_df=meta
        )
        with self._patch_rse_class():
            result = transform_counts(sparse_rse, round_to_integers=False)
            expected = transform_counts(dense_rse, round_to_integers=False)
        assert all(isinstance(d, pd.SparseDtype) for d in result.dtypes)
        np.testing.assert_allclose(
            result.to_numpy(dtype=float), expected.to_numpy(dtype=float)
        )

    def test_basic_auc_scaling_rounded(self) -> None:
        rse = _make_rse()
        with self._patch_rse_class():
            result = transform_counts(
                rse, by="auc", target_read_count=4e7, round_to_integers=True
            )

        assert isinstance(result, pd.DataFrame)
        assert result.shape == (2, 3)
        assert result.iloc[0, 0] == pytest.approx(1.0)

    def test_basic_auc_scaling_unrounded(self) -> None:
        rse = _make_rse()
        with self._patch_rse_class():
            result = transform_counts(
                rse, by="auc", target_read_count=4e7, round_to_integers=False
            )

        assert isinstance(result, pd.DataFrame)
        assert result.iloc[0, 0] == pytest.approx(1.0)

    def test_with_row_and_col_names(self) -> None:
        rse = _make_rse(
            row_names=["G1", "G2"],
            col_names=["S1", "S2", "S3"],
        )
        with self._patch_rse_class():
            result = transform_counts(rse, by="auc")

        assert list(result.index) == ["G1", "G2"]
        assert list(result.columns) == ["S1", "S2", "S3"]

    def test_without_row_and_col_names(self) -> None:
        rse = _make_rse(row_names=None, col_names=None)
        with self._patch_rse_class():
            result = transform_counts(rse, by="auc")

        assert isinstance(result.index, pd.RangeIndex)
        assert isinstance(result.columns, pd.RangeIndex)

    def test_by_mapped_reads(self) -> None:
        rse = _make_rse()
        with self._patch_rse_class():
            result = transform_counts(
                rse,
                by="mapped_reads",
                target_read_count=4e7,
                target_read_length_bp=100.0,
            )
        assert isinstance(result, pd.DataFrame)
        assert result.shape == (2, 3)

    def test_kwargs_forwarded_to_compute_scale_factors(self) -> None:
        rse = _make_rse()
        pe = [False, False, False]
        with self._patch_rse_class():
            result = transform_counts(
                rse,
                by="mapped_reads",
                paired_end_status=pe,
            )
        assert isinstance(result, pd.DataFrame)

    def test_invalid_by_propagates_from_scale_factors(self) -> None:
        rse = _make_rse()
        with self._patch_rse_class():
            with pytest.raises(ValueError, match="'auc' or 'mapped_reads'"):
                transform_counts(rse, by="bad_method")


class TestAssayFrameSparse:
    """_assay_frame hands pandas a canonical matrix."""

    @pytest.mark.parametrize(
        "matrix", _duplicate_entry_matrices(), ids=_DUPLICATE_IDS
    )
    def test_duplicate_entries_are_summed(self, matrix: Any) -> None:
        frame = _assay_frame(matrix, None, None)
        assert all(
            dtype == pd.SparseDtype(np.float64, 0.0) for dtype in frame.dtypes
        )
        np.testing.assert_array_equal(
            frame.to_numpy(dtype=float), _DUPLICATE_DENSE
        )
        assert [frame[c].array.npoints for c in frame.columns] == [2, 1]

    @pytest.mark.parametrize(
        "matrix", _duplicate_entry_matrices(), ids=_DUPLICATE_IDS
    )
    def test_input_matrix_is_not_modified(self, matrix: Any) -> None:
        original = matrix.copy()
        _assay_frame(matrix, None, None)
        for name in ("data", "indices", "indptr", "row", "col"):
            if hasattr(matrix, name):
                np.testing.assert_array_equal(
                    getattr(matrix, name), getattr(original, name)
                )

    def test_canonical_matrix_is_not_copied(self) -> None:
        matrix = sparse.csc_matrix(_DUPLICATE_DENSE)
        frame = _assay_frame(matrix, ["a", "b", "c"], ["s1", "s2"])
        assert np.shares_memory(frame["s1"].array.sp_values, matrix.data)
        assert list(frame.index) == ["a", "b", "c"]


class TestComputeReadCountsSparseAssay:
    @pytest.mark.parametrize(
        "matrix", _duplicate_entry_matrices(), ids=_DUPLICATE_IDS
    )
    def test_duplicate_entries_are_summed(self, matrix: Any) -> None:
        meta = _META.iloc[:2].copy()  # read lengths 100 and 200
        rse = _MockRSE(assays={"raw_counts": matrix}, col_data_df=meta)
        with mock.patch.object(
            _utils,
            "get_ranged_summarizedexperiment_class",
            return_value=_MockRSE,
        ):
            result = compute_read_counts(rse, round_to_integers=False)
        np.testing.assert_allclose(
            result.to_numpy(dtype=float), _DUPLICATE_DENSE / [100.0, 200.0]
        )

    def test_an_unsorted_sparse_assay_is_not_modified(self) -> None:
        """Converting to a sparse frame sorts indices on a private copy."""
        matrix = sparse.csc_matrix(
            (
                np.array([400.0, 100.0, 500.0]),
                np.array([1, 0, 1]),
                np.array([0, 2, 3, 3]),
            ),
            shape=(2, 3),
        )
        assert not matrix.has_sorted_indices
        original_indices = matrix.indices.copy()
        rse = _MockRSE(assays={"raw_counts": matrix}, col_data_df=_META.copy())
        with mock.patch.object(
            _utils,
            "get_ranged_summarizedexperiment_class",
            return_value=_MockRSE,
        ):
            result = compute_read_counts(rse, round_to_integers=False)
        np.testing.assert_array_equal(matrix.indices, original_indices)
        np.testing.assert_allclose(
            result.to_numpy(dtype=float), [[1.0, 0.0, 0.0], [4.0, 2.5, 0.0]]
        )

    def test_a_sparse_assay_is_kept_sparse_and_matches_the_dense_result(
        self,
    ) -> None:
        """Sparse assays must not be densified on the way to read counts."""
        sparse_rse = _MockRSE(
            assays={"raw_counts": sparse.csc_matrix(_RAW.copy())},
            col_data_df=_META.copy(),
        )
        patch_rse_class = mock.patch.object(
            _utils,
            "get_ranged_summarizedexperiment_class",
            return_value=_MockRSE,
        )

        with patch_rse_class:
            result = compute_read_counts(sparse_rse, round_to_integers=False)
            expected = compute_read_counts(_make_rse(), round_to_integers=False)

        assert all(isinstance(dtype, pd.SparseDtype) for dtype in result.dtypes)
        np.testing.assert_allclose(
            result.to_numpy(dtype=float), expected.to_numpy(dtype=float)
        )


# ===========================================================================
# to_anndata
# ===========================================================================


class TestToAnndata:
    """The package's own converter, standing in for BiocPy's broken one."""

    @pytest.mark.requires_anndata
    def test_round_trips_orientation_assays_and_provenance(self) -> None:
        """BiocPy's to_anndata() raises on this metadata; ours must not."""
        anndata = pytest.importorskip("anndata")
        biocframe = pytest.importorskip("biocframe")
        summarizedexperiment = pytest.importorskip("summarizedexperiment")

        feature_ids = ["ENSG00000000003.14", "ENSG00000000005.5"]
        samples = ["SRR387777", "SRR387778", "SRR387779"]
        counts = np.array([[10, 0, 7], [0, 3, 0]], dtype=np.int64)
        experiment = summarizedexperiment.SummarizedExperiment(
            assays={"raw_counts": counts},
            row_data=biocframe.BiocFrame(
                {"feature_id": feature_ids}, row_names=feature_ids
            ),
            column_data=biocframe.BiocFrame(
                {"sample": samples}, row_names=samples
            ),
            metadata={
                "project": "SRP009615",
                "resource_urls": ("http://example.org/a.gz",),
                "metadata_columns": {"qc__auc": ("recount_qc", "auc")},
            },
        )

        with pytest.raises(ValueError, match="mutable mapping"):
            experiment.to_anndata()

        adata = to_anndata(experiment)

        assert isinstance(adata, anndata.AnnData)
        # AnnData is observation-major: samples become rows.
        assert adata.shape == (len(samples), len(feature_ids))
        assert list(adata.obs_names) == samples
        assert list(adata.var_names) == feature_ids
        np.testing.assert_array_equal(
            np.asarray(adata.layers["raw_counts"]), counts.T
        )
        assert adata.uns["project"] == "SRP009615"
        assert adata.uns["resource_urls"] == ["http://example.org/a.gz"]
        assert adata.uns["metadata_columns"] == {
            "qc__auc": ["recount_qc", "auc"]
        }

    @pytest.mark.requires_anndata
    def test_sanitize_for_hdf5_renames_columns_and_uns_keys(
        self, tmp_path: pathlib.Path
    ) -> None:
        """A '/' in a QC field name is a path separator to HDF5."""
        anndata = pytest.importorskip("anndata")
        biocframe = pytest.importorskip("biocframe")
        summarizedexperiment = pytest.importorskip("summarizedexperiment")

        slashed = "recount_qc.star.number_of_splices:_gt/ag"
        experiment = summarizedexperiment.SummarizedExperiment(
            assays={"raw_counts": np.array([[1, 2]], dtype=np.int64)},
            row_data=biocframe.BiocFrame(
                {"feature_id": ["g1"]}, row_names=["g1"]
            ),
            column_data=biocframe.BiocFrame(
                {slashed: [1, 2]}, row_names=["s1", "s2"]
            ),
            metadata={"metadata_columns": {slashed: ("recount_qc", slashed)}},
        )

        adata = to_anndata(experiment, sanitize_for_hdf5=True)

        safe = slashed.replace("/", "_")
        assert safe in adata.obs.columns
        assert slashed not in adata.obs.columns
        assert list(adata.uns["metadata_columns"]) == [safe]
        # The unsanitized name survives in the provenance value.
        assert adata.uns["metadata_columns"][safe][1] == slashed

        out = tmp_path / "sanitized.h5ad"
        adata.write_h5ad(out)
        assert anndata.read_h5ad(out).uns["metadata_columns"][safe][1] == (
            slashed
        )

    @pytest.mark.requires_anndata
    def test_without_sanitizing_the_names_are_left_verbatim(self) -> None:
        biocframe = pytest.importorskip("biocframe")
        summarizedexperiment = pytest.importorskip("summarizedexperiment")

        slashed = "recount_qc.star.number_of_splices:_gt/ag"
        experiment = summarizedexperiment.SummarizedExperiment(
            assays={"raw_counts": np.array([[1, 2]], dtype=np.int64)},
            row_data=biocframe.BiocFrame(
                {"feature_id": ["g1"]}, row_names=["g1"]
            ),
            column_data=biocframe.BiocFrame(
                {slashed: [1, 2]}, row_names=["s1", "s2"]
            ),
            metadata={"metadata_columns": {slashed: ("recount_qc", slashed)}},
        )

        adata = to_anndata(experiment)

        assert slashed in adata.obs.columns
        assert list(adata.uns["metadata_columns"]) == [slashed]

    def test_requires_the_anndata_extra(self) -> None:
        with mock.patch.object(
            _utils,
            "ensure_anndata_support",
            side_effect=ImportError("install recount3[anndata]"),
        ):
            with pytest.raises(ImportError, match="recount3\\[anndata\\]"):
                to_anndata(mock.MagicMock())
