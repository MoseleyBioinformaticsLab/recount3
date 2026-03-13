from __future__ import annotations

import logging
from unittest import mock

import numpy as np
import pandas as pd
import pytest

import recount3._utils as _utils
from recount3.se import (
    _expand_sra_attributes_df,
    _resolve_annotation_extension,
    build_summarized_experiment,
    build_ranged_summarized_experiment,
    create_ranged_summarized_experiment,
    create_rse,
    expand_sra_attributes,
    compute_read_counts,
    compute_tpm,
    is_paired_end,
    compute_scale_factors,
    transform_counts,
)


class _MockBiocFrame:
    def __init__(self, df: pd.DataFrame) -> None:
        self._df = df

    def to_pandas(self) -> pd.DataFrame:
        return self._df.copy()

    @classmethod
    def from_pandas(cls, df: pd.DataFrame) -> "_MockBiocFrame":
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

    @property
    def col_data(self) -> _MockBiocFrame:
        return self._bioc_col_data

    @property
    def column_data(self) -> _MockBiocFrame:
        return self._bioc_col_data

    def set_column_data(self, new_bioc_col_data: _MockBiocFrame) -> "_MockSE":
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

    def set_column_data(self, new_bioc_col_data: _MockBiocFrame) -> "_MockRSE":
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


def _patch_biocpy():
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
            {"sra.sample_attributes": ["age;;30|disease;;Cancer", "age;;45|disease;;Control"]}
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
        with mock.patch("recount3.se.annotation_ext", return_value="G029") as m:
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
            autoload=True,
        )

        bundle.to_summarized_experiment.assert_called_once_with(
            genomic_unit="gene",
            annotation_extension="G026",
            assay_name="raw_counts",
            join_policy="inner",
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
            autoload=False,
            allow_fallback_to_se=True,
        )

        bundle.to_ranged_summarized_experiment.assert_called_once_with(
            genomic_unit="exon",
            annotation_extension="G026",
            prefer_rr_junction_coordinates=False,
            assay_name="raw_counts",
            join_policy="outer",
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

    def test_default_junction_extensions(self) -> None:
        """When junction_extensions is None, ('MM',) is used."""
        bundle = self._make_bundle_mock()
        with mock.patch("recount3.se.R3ResourceBundle") as MockBundle:
            MockBundle.discover.return_value = bundle
            create_ranged_summarized_experiment(
                project="SRP009615",
                genomic_unit="gene",
                organism="human",
            )
        call_kwargs = MockBundle.discover.call_args.kwargs
        assert call_kwargs["junction_exts"] == ("MM",)

    def test_explicit_junction_extensions(self) -> None:
        """Provided junction_extensions are tupled with strip()."""
        bundle = self._make_bundle_mock()
        with mock.patch("recount3.se.R3ResourceBundle") as MockBundle:
            MockBundle.discover.return_value = bundle
            create_ranged_summarized_experiment(
                project="SRP009615",
                genomic_unit="junction",
                organism="human",
                junction_extensions=["MM", " RR "],
            )
        call_kwargs = MockBundle.discover.call_args.kwargs
        assert call_kwargs["junction_exts"] == ("MM", "RR")

    def test_annotations_default_when_ann_ext_none(self) -> None:
        """Junction unit -> ann_ext=None -> annotations='default'."""
        bundle = self._make_bundle_mock()
        with mock.patch("recount3.se.R3ResourceBundle") as MockBundle:
            MockBundle.discover.return_value = bundle
            create_ranged_summarized_experiment(
                project="SRP009615",
                genomic_unit="junction",
                organism="human",
            )
        call_kwargs = MockBundle.discover.call_args.kwargs
        assert call_kwargs["annotations"] == "default"

    def test_annotations_tuple_when_ann_ext_provided(self) -> None:
        """Gene unit with default human -> ann_ext='G026' -> annotations=('G026',)."""
        bundle = self._make_bundle_mock()
        with mock.patch("recount3.se.R3ResourceBundle") as MockBundle:
            MockBundle.discover.return_value = bundle
            create_ranged_summarized_experiment(
                project="SRP009615",
                genomic_unit="gene",
                organism="human",
            )
        call_kwargs = MockBundle.discover.call_args.kwargs
        assert call_kwargs["annotations"] == ("G026",)

    def test_all_kwargs_forwarded_to_rse(self) -> None:
        bundle = self._make_bundle_mock()
        with mock.patch("recount3.se.R3ResourceBundle") as MockBundle:
            MockBundle.discover.return_value = bundle
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


class TestExpandSraAttributes:
    def test_dataframe_missing_column_returns_copy(self) -> None:
        df = pd.DataFrame({"other": [1, 2]})
        result = expand_sra_attributes(df)
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ["other"]

    def test_dataframe_with_column_expands(self) -> None:
        df = pd.DataFrame({"sra.sample_attributes": ["age;;30|disease;;Cancer"]})
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

    def _biocpy_patches(self):
        return mock.patch.multiple(
            _utils,
            get_biocframe_class=mock.MagicMock(return_value=_MockBiocFrame),
            get_summarizedexperiment_class=mock.MagicMock(return_value=_MockSE),
            get_ranged_summarizedexperiment_class=mock.MagicMock(return_value=_MockRSE),
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
        assert "not present" in caplog.text.lower() or "not present" in caplog.text

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
    def _patch_rse_class(self):
        return mock.patch.object(
            _utils,
            "get_ranged_summarizedexperiment_class",
            return_value=_MockRSE,
        )

    def test_non_rse_raises_value_error(self) -> None:
        with self._patch_rse_class():
            with pytest.raises(ValueError, match="RangedSummarizedExperiment"):
                compute_read_counts("not_an_rse")

    def test_round_to_integers_not_bool_raises_type_error(self) -> None:
        rse = _make_rse()
        with self._patch_rse_class():
            with pytest.raises(TypeError, match="bool"):
                compute_read_counts(rse, round_to_integers=1)  # type: ignore[arg-type]

    def test_missing_avg_length_column_raises_value_error(self) -> None:
        meta = pd.DataFrame({"external_id": ["S1"]})
        rse = _make_rse(counts=np.array([[10.0]]), col_data_df=meta)
        with self._patch_rse_class():
            with pytest.raises(ValueError, match="Required metadata column"):
                compute_read_counts(rse)

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


class TestComputeTpm:
    """Tests for compute_tpm."""

    def _patch_rse_class(self):
        return mock.patch.object(
            _utils,
            "get_ranged_summarizedexperiment_class",
            return_value=_MockRSE,
        )

    def test_non_rse_raises_type_error(self) -> None:
        with self._patch_rse_class():
            with pytest.raises(TypeError, match="RangedSummarizedExperiment"):
                compute_tpm("not_rse")

    def test_missing_row_ranges_raises_value_error(self) -> None:
        """No row_ranges attribute -> AttributeError -> wrapped ValueError."""
        rse = _make_rse()
        with self._patch_rse_class():
            with pytest.raises(ValueError, match="feature widths"):
                compute_tpm(rse)

    def test_basic_tpm_computation(self) -> None:
        """Verify TPM = RPK / (sum(RPK)/1e6) per sample, feature-wise."""
        widths = [1000, 2000]
        rse = _make_rse(widths=widths)

        with self._patch_rse_class():
            result = compute_tpm(rse)

        assert isinstance(result, pd.DataFrame)
        assert result.shape == (2, 3)
        for col in result.columns:
            total = result[col].sum()
            assert total == pytest.approx(1_000_000.0, rel=1e-5)

    def test_tpm_values_non_negative(self) -> None:
        widths = [500, 1000, 1500]
        counts = np.array([[10.0, 20.0, 30.0], [40.0, 50.0, 60.0], [70.0, 80.0, 90.0]])
        meta = _META.copy()
        rse = _make_rse(counts=counts, col_data_df=meta, widths=widths)

        with self._patch_rse_class():
            result = compute_tpm(rse)

        assert (result.values >= 0).all()


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

    def test_mixed_with_invalid_ratio(self, caplog) -> None:
        meta = self._meta([100.0, 300.0, 200.0], [100.0, 100.0, 100.0])
        with caplog.at_level(logging.WARNING):
            result = is_paired_end(meta)
        assert result.iloc[0] == False
        assert pd.isna(result.iloc[1])
        assert result.iloc[2] == True

    def test_indexed_by_external_id(self) -> None:
        meta = self._meta([100.0, 200.0], [100.0, 100.0], external_ids=["ERR001", "ERR002"])
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
        result = compute_scale_factors(_META.copy(), by="auc", target_read_count=4e7)
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
        expected_s1 = 4e7 * 100.0 * 1.0 / (1e6 * (100.0 ** 2))
        assert result.iloc[0] == pytest.approx(expected_s1)

    def test_mapped_reads_inferred_paired_end(self) -> None:
        """S2 has ratio=2 (PE), S1/S3 have ratio=1 (SE)."""
        result = compute_scale_factors(
            _META.copy(),
            by="mapped_reads",
            target_read_count=4e7,
            target_read_length_bp=100.0,
        )
        s1 = 4e7 * 100.0 * 1.0 / (1e6 * (100.0 ** 2))
        s2 = 4e7 * 100.0 * 2.0 / (2e6 * (200.0 ** 2))
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

    def _patch_rse_class(self):
        return mock.patch.object(
            _utils,
            "get_ranged_summarizedexperiment_class",
            return_value=_MockRSE,
        )

    def test_non_rse_raises_value_error(self) -> None:
        with self._patch_rse_class():
            with pytest.raises(ValueError, match="RangedSummarizedExperiment"):
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

    def test_basic_auc_scaling_rounded(self) -> None:
        rse = _make_rse()
        with self._patch_rse_class():
            result = transform_counts(rse, by="auc", target_read_count=4e7, round_to_integers=True)

        assert isinstance(result, pd.DataFrame)
        assert result.shape == (2, 3)
        assert result.iloc[0, 0] == pytest.approx(1.0)

    def test_basic_auc_scaling_unrounded(self) -> None:
        rse = _make_rse()
        with self._patch_rse_class():
            result = transform_counts(rse, by="auc", target_read_count=4e7, round_to_integers=False)

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
