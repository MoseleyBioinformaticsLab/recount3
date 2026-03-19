# pylint: disable=redefined-outer-name

from __future__ import annotations

import shutil
from pathlib import Path
from unittest import mock

import pandas as pd
import pytest

from recount3._descriptions import R3ResourceDescription
from recount3.config import Config
from recount3.resource import R3Resource
from recount3.search import (
    _ANN_EXT_HUMAN,
    _ANN_EXT_MOUSE,
    _ANNOTATION_NAME_TO_EXT_HUMAN,
    _ANNOTATION_NAME_TO_EXT_MOUSE,
    _build_param_grid,
    _make_resources,
    _normalize_organism_name,
    _posix_basename,
    _posix_dirname,
    _resolve_annotation_exts,
    _strip_metadata_column_prefix,
    annotation_ext,
    annotation_options,
    _normalize_to_tuple,
    available_projects,
    available_samples,
    create_sample_project_lists,
    match_spec,
    project_homes,
    samples_for_project,
    search_annotations,
    search_bigwig_files,
    search_count_files_gene_or_exon,
    search_count_files_junctions,
    search_data_source_metadata,
    search_data_sources,
    search_metadata_files,
    search_project_all,
)

_DATA_DIR = Path(__file__).parent / "data"
_MIRROR = _DATA_DIR / "recount3_mirror" / "recount3"
_META_HUMAN_MD_GZ = (
    _MIRROR
    / "human"
    / "data_sources"
    / "sra"
    / "metadata"
    / "sra.recount_project.MD.gz"
)
_META_MOUSE_MD_GZ = (
    _MIRROR
    / "mouse"
    / "data_sources"
    / "sra"
    / "metadata"
    / "sra.recount_project.MD.gz"
)
_BASE_URL = "https://example.org/recount3/"


@pytest.fixture()
def cfg(tmp_path: Path) -> Config:
    """Minimal Config pointing at a temporary cache directory."""
    return Config(
        base_url=_BASE_URL,
        timeout=30,
        insecure_ssl=False,
        max_retries=1,
        user_agent="test/1.0",
        cache_dir=tmp_path / "cache",
        cache_disabled=False,
        chunk_size=1024,
    )


def _seed_desc(desc_params: dict, cfg: Config, src: Path) -> None:
    """Seed the cache for the resource described by desc_params."""
    desc = R3ResourceDescription(**desc_params)
    res = R3Resource(desc, config=cfg)
    dest = res._cached_path()
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)


def _seed_human_sra(cfg: Config) -> None:
    _seed_desc(
        {
            "resource_type": "data_source_metadata",
            "organism": "human",
            "data_source": "sra",
        },
        cfg,
        _META_HUMAN_MD_GZ,
    )


def _seed_mouse_sra(cfg: Config) -> None:
    _seed_desc(
        {
            "resource_type": "data_source_metadata",
            "organism": "mouse",
            "data_source": "sra",
        },
        cfg,
        _META_MOUSE_MD_GZ,
    )


class TestAsTuple:
    def test_string_becomes_single_element_tuple(self) -> None:
        assert _normalize_to_tuple("human") == ("human",)

    def test_list_becomes_tuple(self) -> None:
        assert _normalize_to_tuple(["a", "b", "c"]) == ("a", "b", "c")

    def test_generator_becomes_tuple(self) -> None:
        assert _normalize_to_tuple(x for x in ["x", "y"]) == ("x", "y")

    def test_empty_list_becomes_empty_tuple(self) -> None:
        assert _normalize_to_tuple([]) == ()


class TestMatchSpec:
    def test_none_spec_always_matches(self) -> None:
        assert match_spec("anything", None) is True

    def test_callable_called_with_value_true(self) -> None:
        assert match_spec("hello", lambda v: v == "hello") is True

    def test_callable_called_with_value_false(self) -> None:
        assert match_spec("world", lambda v: v == "hello") is False

    def test_callable_returns_truthy_coerced_to_bool(self) -> None:
        assert match_spec("x", lambda v: v) is True

    def test_callable_returns_falsy_coerced_to_bool(self) -> None:
        assert match_spec("", lambda v: v) is False

    def test_string_spec_value_in_tuple(self) -> None:
        assert match_spec("sra", "sra") is True

    def test_string_spec_value_not_in_tuple(self) -> None:
        assert match_spec("gtex", "sra") is False

    def test_iterable_spec_value_present(self) -> None:
        assert match_spec("gtex", ["sra", "gtex"]) is True

    def test_iterable_spec_value_absent(self) -> None:
        assert match_spec("tcga", ["sra", "gtex"]) is False

    def test_empty_iterable_spec_never_matches(self) -> None:
        assert match_spec("anything", []) is False

    def test_value_none_with_callable(self) -> None:
        assert match_spec(None, lambda v: v is None) is True


class TestBuildParamGrid:
    def test_single_value_each_produces_one_combo(self) -> None:
        grid = _build_param_grid(
            "annotations", organism="human", genomic_unit="gene"
        )
        assert grid == [
            {
                "resource_type": "annotations",
                "organism": "human",
                "genomic_unit": "gene",
            }
        ]

    def test_multiple_values_cartesian_product(self) -> None:
        grid = _build_param_grid(
            "annotations", organism=["human", "mouse"], genomic_unit="gene"
        )
        assert len(grid) == 2
        organisms = {g["organism"] for g in grid}
        assert organisms == {"human", "mouse"}

    def test_two_multi_value_fields(self) -> None:
        grid = _build_param_grid(
            "count_files_gene_or_exon",
            organism=["human", "mouse"],
            genomic_unit=["gene", "exon"],
        )
        assert len(grid) == 4
        assert all(
            g["resource_type"] == "count_files_gene_or_exon" for g in grid
        )

    def test_resource_type_always_in_dict(self) -> None:
        grid = _build_param_grid("data_sources", organism="human")
        assert grid[0]["resource_type"] == "data_sources"


class TestMakeResources:
    def test_valid_params_return_r3resource_list(self) -> None:
        params = [{"resource_type": "data_sources", "organism": "human"}]
        result = _make_resources(params)
        assert len(result) == 1
        assert isinstance(result[0], R3Resource)

    def test_invalid_params_strict_raises(self) -> None:
        params = [
            {"resource_type": "data_sources", "organism": "invalid_organism"}
        ]
        with pytest.raises(ValueError):
            _make_resources(params, strict=True)

    def test_invalid_params_not_strict_skips(self) -> None:
        params = [
            {"resource_type": "data_sources", "organism": "invalid_organism"}
        ]
        result = _make_resources(params, strict=False)
        assert result == []

    def test_deduplication_same_url_skips_second(self) -> None:
        params = [
            {"resource_type": "data_sources", "organism": "human"},
            {"resource_type": "data_sources", "organism": "human"},
        ]
        result = _make_resources(params, deduplicate=True)
        assert len(result) == 1

    def test_no_deduplication_allows_duplicates(self) -> None:
        params = [
            {"resource_type": "data_sources", "organism": "human"},
            {"resource_type": "data_sources", "organism": "human"},
        ]
        result = _make_resources(params, deduplicate=False)
        assert len(result) == 2

    def test_url_none_resources_pass_through_dedup(self) -> None:
        """url=None bypasses dedup; two such resources pass through."""
        mock_res1 = mock.MagicMock(spec=R3Resource)
        mock_res1.url = None
        mock_res2 = mock.MagicMock(spec=R3Resource)
        mock_res2.url = None

        params = [
            {"resource_type": "data_sources", "organism": "human"},
            {"resource_type": "data_sources", "organism": "human"},
        ]
        with (
            mock.patch("recount3.search.R3ResourceDescription"),
            mock.patch(
                "recount3.search.R3Resource",
                side_effect=[mock_res1, mock_res2],
            ),
        ):
            result = _make_resources(params, deduplicate=True)
        assert len(result) == 2

    def test_empty_input_returns_empty_list(self) -> None:
        assert _make_resources([]) == []


class TestNormalizeOrganismName:
    def test_lowercase_human_accepted(self) -> None:
        assert _normalize_organism_name("human") == "human"

    def test_uppercase_normalized(self) -> None:
        assert _normalize_organism_name("HUMAN") == "human"

    def test_mixed_case_with_whitespace(self) -> None:
        assert _normalize_organism_name("  Mouse  ") == "mouse"

    def test_invalid_organism_raises(self) -> None:
        with pytest.raises(ValueError, match="Unsupported organism"):
            _normalize_organism_name("rat")


class TestStripMdPrefix:
    def test_strips_dot_prefix(self) -> None:
        df = pd.DataFrame({"sra.study": ["a"], "sra.project": ["b"]})
        result = _strip_metadata_column_prefix(df)
        assert list(result.columns) == ["study", "project"]

    def test_no_prefix_unchanged(self) -> None:
        df = pd.DataFrame({"study": ["a"], "project": ["b"]})
        result = _strip_metadata_column_prefix(df)
        assert list(result.columns) == ["study", "project"]

    def test_original_frame_not_mutated(self) -> None:
        df = pd.DataFrame({"sra.col": [1]})
        _strip_metadata_column_prefix(df)
        assert list(df.columns) == ["sra.col"]

    def test_multiple_dots_keeps_last_segment(self) -> None:
        df = pd.DataFrame({"a.b.c": [1]})
        result = _strip_metadata_column_prefix(df)
        assert list(result.columns) == ["c"]


class TestPathHelpers:
    def test_basename_empty_string(self) -> None:
        assert _posix_basename("") == ""

    def test_basename_normal_path(self) -> None:
        assert _posix_basename("data_sources/sra") == "sra"

    def test_basename_trailing_slash_stripped(self) -> None:
        assert _posix_basename("data_sources/sra/") == "sra"

    def test_basename_single_component(self) -> None:
        assert _posix_basename("human") == "human"

    def test_dirname_empty_string(self) -> None:
        assert _posix_dirname("") == ""

    def test_dirname_normal_path(self) -> None:
        assert _posix_dirname("data_sources/sra") == "data_sources"

    def test_dirname_trailing_slash_stripped(self) -> None:
        assert _posix_dirname("data_sources/sra/") == "data_sources"

    def test_dirname_single_component_is_empty(self) -> None:
        assert _posix_dirname("human") == ""


class TestSearchAnnotations:
    def test_returns_r3resource_list(self) -> None:
        result = search_annotations(
            organism="human", genomic_unit="gene", annotation_extension="G026"
        )
        assert len(result) == 1
        assert isinstance(result[0], R3Resource)

    def test_url_contains_correct_path(self) -> None:
        result = search_annotations(
            organism="human", genomic_unit="exon", annotation_extension="G026"
        )
        assert "exon_sums" in (result[0].url or "")

    def test_multiple_annotations_cartesian_product(self) -> None:
        result = search_annotations(
            organism="human",
            genomic_unit="gene",
            annotation_extension=["G026", "G029"],
        )
        assert len(result) == 2

    def test_strict_false_skips_invalid(self) -> None:
        result = search_annotations(
            organism="invalid",
            genomic_unit="gene",
            annotation_extension="G026",
            strict=False,
        )
        assert result == []

    def test_deduplicate_false_allows_dupes(self) -> None:
        result = search_annotations(
            organism="human",
            genomic_unit="gene",
            annotation_extension=["G026", "G026"],
            deduplicate=False,
        )
        assert len(result) == 2


class TestSearchCountFilesGeneOrExon:
    def test_basic_gene_resource(self) -> None:
        result = search_count_files_gene_or_exon(
            organism="human",
            data_source="sra",
            genomic_unit="gene",
            project="SRP009615",
        )
        assert len(result) == 1
        assert "gene_sums" in (result[0].url or "")

    def test_default_annotation_extension(self) -> None:
        result = search_count_files_gene_or_exon(
            organism="human",
            data_source="sra",
            genomic_unit="gene",
            project="SRP009615",
        )
        assert "G026" in (result[0].url or "")

    def test_multiple_projects_multiple_resources(self) -> None:
        result = search_count_files_gene_or_exon(
            organism="human",
            data_source="sra",
            genomic_unit="gene",
            project=["SRP009615", "SRP014565"],
        )
        assert len(result) == 2

    def test_strict_false_skips_invalid(self) -> None:
        result = search_count_files_gene_or_exon(
            organism="bad",
            data_source="sra",
            genomic_unit="gene",
            project="SRP009615",
            strict=False,
        )
        assert result == []


class TestSearchCountFilesJunctions:
    def test_basic_junction_resource(self) -> None:
        result = search_count_files_junctions(
            organism="human", data_source="sra", project="SRP009615"
        )
        assert len(result) == 1
        assert "junctions" in (result[0].url or "")

    def test_multiple_extensions(self) -> None:
        result = search_count_files_junctions(
            organism="human",
            data_source="sra",
            project="SRP009615",
            junction_extension=["MM", "RR"],
        )
        assert len(result) == 2

    def test_strict_false_skips_invalid(self) -> None:
        result = search_count_files_junctions(
            organism="bad", data_source="sra", project="P", strict=False
        )
        assert result == []


class TestSearchMetadataFiles:
    def test_returns_r3resource(self) -> None:
        result = search_metadata_files(
            organism="human",
            data_source="sra",
            table_name="recount_project",
            project="SRP009615",
        )
        assert len(result) == 1
        assert "metadata" in (result[0].url or "")

    def test_multiple_table_names(self) -> None:
        result = search_metadata_files(
            organism="human",
            data_source="sra",
            table_name=["recount_project", "recount_qc"],
            project="SRP009615",
        )
        assert len(result) == 2


class TestSearchBigwigFiles:
    def test_returns_r3resource(self) -> None:
        result = search_bigwig_files(
            organism="human",
            data_source="sra",
            project="SRP009615",
            sample="SRR387777",
        )
        assert len(result) == 1
        assert "base_sums" in (result[0].url or "")

    def test_multiple_samples(self) -> None:
        result = search_bigwig_files(
            organism="human",
            data_source="sra",
            project="SRP009615",
            sample=["SRR387777", "SRR387778"],
        )
        assert len(result) == 2


class TestSearchDataSources:
    def test_returns_homes_index_resource(self) -> None:
        result = search_data_sources(organism="human")
        assert len(result) == 1
        assert "homes_index" in (result[0].url or "")

    def test_multiple_organisms(self) -> None:
        result = search_data_sources(organism=["human", "mouse"])
        assert len(result) == 2

    def test_strict_false_skips_invalid(self) -> None:
        result = search_data_sources(organism="invalid", strict=False)
        assert result == []

    def test_deduplicate_param_exposed(self) -> None:
        result = search_data_sources(
            organism=["human", "human"], deduplicate=False
        )
        assert len(result) == 2


class TestSearchDataSourceMetadata:
    def test_returns_r3resource(self) -> None:
        result = search_data_source_metadata(
            organism="human", data_source="sra"
        )
        assert len(result) == 1
        assert "recount_project.MD.gz" in (result[0].url or "")

    def test_multiple_sources_multiple_resources(self) -> None:
        result = search_data_source_metadata(
            organism="human", data_source=["sra", "gtex"]
        )
        assert len(result) == 2

    def test_strict_false_skips_invalid(self) -> None:
        result = search_data_source_metadata(
            organism="bad", data_source="sra", strict=False
        )
        assert result == []


class TestAnnotationOptions:
    def test_human_returns_all_annotations(self) -> None:
        opts = annotation_options("human")
        assert opts == dict(_ANNOTATION_NAME_TO_EXT_HUMAN)

    def test_mouse_returns_mouse_annotations(self) -> None:
        opts = annotation_options("mouse")
        assert opts == dict(_ANNOTATION_NAME_TO_EXT_MOUSE)

    def test_case_insensitive(self) -> None:
        assert annotation_options("Human") == annotation_options("human")

    def test_invalid_organism_raises(self) -> None:
        with pytest.raises(ValueError, match="Unsupported organism"):
            annotation_options("rat")

    def test_returns_new_dict_each_call(self) -> None:
        a = annotation_options("human")
        b = annotation_options("human")
        assert a is not b


class TestAnnotationExt:
    def test_canonical_name_human(self) -> None:
        assert annotation_ext("human", "gencode_v26") == "G026"

    def test_all_human_canonical_names(self) -> None:
        for name, ext in _ANNOTATION_NAME_TO_EXT_HUMAN.items():
            assert annotation_ext("human", name) == ext

    def test_extension_code_passthrough_human(self) -> None:
        assert annotation_ext("human", "G026") == "G026"

    def test_all_human_extension_codes(self) -> None:
        for ext in _ANN_EXT_HUMAN:
            assert annotation_ext("human", ext) == ext

    def test_canonical_name_mouse(self) -> None:
        assert annotation_ext("mouse", "gencode_v23") == "M023"

    def test_extension_code_mouse(self) -> None:
        assert annotation_ext("mouse", "M023") == "M023"

    def test_case_insensitive_name(self) -> None:
        assert annotation_ext("human", "GENCODE_V26") == "G026"

    def test_case_insensitive_organism(self) -> None:
        assert annotation_ext("HUMAN", "G026") == "G026"

    def test_empty_annotation_raises(self) -> None:
        with pytest.raises(ValueError, match="non-empty"):
            annotation_ext("human", "")

    def test_whitespace_only_annotation_raises(self) -> None:
        with pytest.raises(ValueError, match="non-empty"):
            annotation_ext("human", "   ")

    def test_unknown_annotation_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown annotation"):
            annotation_ext("human", "ZZZ999")

    def test_invalid_organism_raises(self) -> None:
        with pytest.raises(ValueError, match="Unsupported organism"):
            annotation_ext("rat", "G026")

    def test_error_message_lists_valid_options(self) -> None:
        with pytest.raises(ValueError, match="gencode_v26"):
            annotation_ext("human", "bad_name")


class TestResolveAnnotationExts:
    def test_explicit_exts_win(self) -> None:
        result = _resolve_annotation_exts("human", "all", ["G029", "F006"])
        assert result == ("G029", "F006")

    def test_explicit_exts_strips_whitespace(self) -> None:
        result = _resolve_annotation_exts("human", None, [" G026 "])
        assert result == ("G026",)

    def test_annotations_none_human_default(self) -> None:
        assert _resolve_annotation_exts("human", None, None) == ("G026",)

    def test_annotations_default_string_human(self) -> None:
        assert _resolve_annotation_exts("human", "default", None) == ("G026",)

    def test_annotations_none_mouse_default(self) -> None:
        assert _resolve_annotation_exts("mouse", None, None) == ("M023",)

    def test_annotations_all_human(self) -> None:
        assert _resolve_annotation_exts("human", "all", None) == _ANN_EXT_HUMAN

    def test_annotations_all_mouse(self) -> None:
        assert _resolve_annotation_exts("mouse", "all", None) == _ANN_EXT_MOUSE

    def test_annotations_comma_separated_string(self) -> None:
        result = _resolve_annotation_exts("human", "G026, G029", None)
        assert result == ("G026", "G029")

    def test_annotations_iterable(self) -> None:
        result = _resolve_annotation_exts("human", ["G026", "G029"], None)
        assert result == ("G026", "G029")

    def test_invalid_organism_raises(self) -> None:
        with pytest.raises(ValueError, match="Unsupported organism"):
            _resolve_annotation_exts("rat", None, None)

    def test_empty_explicit_exts_falls_through(self) -> None:
        """Empty explicit_exts is falsy -> falls through to annotations."""
        result = _resolve_annotation_exts("human", None, [])
        assert result == ("G026",)


class TestAvailableSamples:
    def test_invalid_organism_raises(self) -> None:
        with pytest.raises(ValueError, match="Unsupported organism"):
            available_samples(organism="dolphin")

    def test_invalid_data_source_raises(self) -> None:
        with pytest.raises(ValueError, match="Unsupported data_sources"):
            available_samples(organism="human", data_sources=["bogus"])

    def test_empty_data_sources_strict_raises(self) -> None:
        with pytest.raises(ValueError, match="No data_sources"):
            available_samples(organism="human", data_sources=[])

    def test_empty_data_sources_not_strict_returns_empty(self) -> None:
        df = available_samples(organism="human", data_sources=[], strict=False)
        assert df.empty

    def test_no_resources_not_strict_returns_empty(self) -> None:
        """search_data_source_metadata returns [] with strict=False -> empty."""
        with mock.patch(
            "recount3.search.search_data_source_metadata", return_value=[]
        ):
            df = available_samples(
                organism="human", data_sources=["sra"], strict=False
            )
        assert df.empty

    def test_no_resources_strict_raises(self) -> None:
        with mock.patch(
            "recount3.search.search_data_source_metadata", return_value=[]
        ):
            with pytest.raises(ValueError, match="No data-source metadata"):
                available_samples(organism="human", data_sources=["sra"])

    def test_load_error_raises_runtime(self, cfg: Config) -> None:
        """All resources fail to load -> RuntimeError."""
        mock_res = mock.MagicMock()
        mock_res.url = "https://example.org/fake"
        mock_res.load.side_effect = OSError("network error")
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            with pytest.raises(RuntimeError, match="Failed to load"):
                available_samples(organism="human", data_sources=["sra"])

    def test_load_returns_non_dataframe_skipped_strict_raises(self) -> None:
        """Non-DataFrame load() skipped; empty result with strict raises."""
        mock_res = mock.MagicMock()
        mock_res.url = "https://example.org/fake"
        mock_res.load.return_value = "not a dataframe"
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            with pytest.raises(ValueError, match="no usable"):
                available_samples(organism="human", data_sources=["sra"])

    def test_load_returns_non_dataframe_not_strict_returns_empty(self) -> None:
        mock_res = mock.MagicMock()
        mock_res.url = "https://example.org/fake"
        mock_res.load.return_value = "not a dataframe"
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            df = available_samples(
                organism="human", data_sources=["sra"], strict=False
            )
        assert df.empty

    def test_returns_dataframe_with_expected_columns(self, cfg: Config) -> None:
        _seed_human_sra(cfg)
        with mock.patch("recount3.resource.default_config", return_value=cfg):
            df = available_samples(organism="human", data_sources=["sra"])
        assert isinstance(df, pd.DataFrame)
        assert not df.empty
        assert "project" in df.columns
        assert "organism" in df.columns

    def test_organism_normalized_homo_sapiens_to_human(
        self, cfg: Config
    ) -> None:
        _seed_human_sra(cfg)
        with mock.patch("recount3.resource.default_config", return_value=cfg):
            df = available_samples(organism="human", data_sources=["sra"])
        assert (df["organism"] == "human").all()

    def test_organism_column_absent_filled_from_arg(self) -> None:
        """No 'organism' column in metadata -> fill from the organism arg."""
        sample_df = pd.DataFrame(
            {
                "external_id": ["SRR1"],
                "study": ["SRP1"],
                "project": ["SRP1"],
            }
        )
        mock_res = mock.MagicMock()
        mock_res.url = "https://example.org/fake"
        mock_res.load.return_value = sample_df
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            df = available_samples(organism="mouse", data_sources=["sra"])
        assert (df["organism"] == "mouse").all()

    def test_study_renamed_to_project_when_project_absent(self) -> None:
        sample_df = pd.DataFrame(
            {
                "external_id": ["SRR1"],
                "study": ["SRP1"],
                "organism": ["human"],
            }
        )
        mock_res = mock.MagicMock()
        mock_res.url = "https://example.org/fake"
        mock_res.load.return_value = sample_df
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            df = available_samples(organism="human", data_sources=["sra"])
        assert "project" in df.columns
        assert "study" not in df.columns

    def test_project_and_study_identical_study_dropped(self) -> None:
        sample_df = pd.DataFrame(
            {
                "external_id": ["SRR1"],
                "project": ["SRP1"],
                "study": ["SRP1"],
                "organism": ["human"],
            }
        )
        mock_res = mock.MagicMock()
        mock_res.url = "https://example.org/fake"
        mock_res.load.return_value = sample_df
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            df = available_samples(organism="human", data_sources=["sra"])
        assert "study" not in df.columns

    def test_project_study_mismatch_raises(self) -> None:
        sample_df = pd.DataFrame(
            {
                "external_id": ["SRR1"],
                "project": ["SRP1"],
                "study": ["SRP_DIFFERENT"],
                "organism": ["human"],
            }
        )
        mock_res = mock.MagicMock()
        mock_res.url = "https://example.org/fake"
        mock_res.load.return_value = sample_df
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            with pytest.raises(ValueError, match="not identical"):
                available_samples(organism="human", data_sources=["sra"])

    def test_file_source_normalized_to_basename(self) -> None:
        sample_df = pd.DataFrame(
            {
                "external_id": ["SRR1"],
                "project": ["SRP1"],
                "organism": ["human"],
                "file_source": ["data_sources/sra"],
            }
        )
        mock_res = mock.MagicMock()
        mock_res.url = "https://example.org/fake"
        mock_res.load.return_value = sample_df
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            df = available_samples(organism="human", data_sources=["sra"])
        assert df["file_source"].iloc[0] == "sra"

    def test_metadata_source_renamed_to_project_home(self) -> None:
        sample_df = pd.DataFrame(
            {
                "external_id": ["SRR1"],
                "project": ["SRP1"],
                "organism": ["human"],
                "metadata_source": ["data_sources/sra"],
            }
        )
        mock_res = mock.MagicMock()
        mock_res.url = "https://example.org/fake"
        mock_res.load.return_value = sample_df
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            df = available_samples(organism="human", data_sources=["sra"])
        assert "project_home" in df.columns
        assert "metadata_source" not in df.columns

    def test_project_home_derives_project_type(self) -> None:
        sample_df = pd.DataFrame(
            {
                "external_id": ["SRR1"],
                "project": ["SRP1"],
                "organism": ["human"],
                "project_home": ["data_sources/sra"],
            }
        )
        mock_res = mock.MagicMock()
        mock_res.url = "https://example.org/fake"
        mock_res.load.return_value = sample_df
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            df = available_samples(organism="human", data_sources=["sra"])
        assert "project_type" in df.columns
        assert df["project_type"].iloc[0] == "data_sources"

    def test_rail_id_dropped(self) -> None:
        sample_df = pd.DataFrame(
            {
                "rail_id": [123],
                "external_id": ["SRR1"],
                "project": ["SRP1"],
                "organism": ["human"],
            }
        )
        mock_res = mock.MagicMock()
        mock_res.url = "https://example.org/fake"
        mock_res.load.return_value = sample_df
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            df = available_samples(organism="human", data_sources=["sra"])
        assert "rail_id" not in df.columns

    def test_preferred_columns_ordered_first(self) -> None:
        sample_df = pd.DataFrame(
            {
                "extra_col": ["x"],
                "external_id": ["SRR1"],
                "project": ["SRP1"],
                "organism": ["human"],
            }
        )
        mock_res = mock.MagicMock()
        mock_res.url = "https://example.org/fake"
        mock_res.load.return_value = sample_df
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            df = available_samples(organism="human", data_sources=["sra"])
        assert df.columns[0] == "external_id"
        assert "extra_col" in df.columns

    def test_data_sources_none_uses_all(self, cfg: Config) -> None:
        _seed_human_sra(cfg)
        # gtex and tcga not seeded -> strict=False
        with mock.patch("recount3.resource.default_config", return_value=cfg):
            df = available_samples(
                organism="human", data_sources=None, strict=False
            )
        assert isinstance(df, pd.DataFrame)


class TestAvailableProjects:
    def test_empty_samples_returns_empty_dataframe(self) -> None:
        with mock.patch(
            "recount3.search.available_samples", return_value=pd.DataFrame()
        ):
            df = available_projects(organism="human", strict=False)
        assert list(df.columns) == [
            "project",
            "organism",
            "file_source",
            "project_home",
            "project_type",
            "n_samples",
        ]

    def test_one_row_per_project(self, cfg: Config) -> None:
        _seed_human_sra(cfg)
        with mock.patch("recount3.resource.default_config", return_value=cfg):
            df = available_projects(organism="human", data_sources=["sra"])
        assert isinstance(df, pd.DataFrame)
        assert df["project"].nunique() == len(df)

    def test_n_samples_correctly_counted(self) -> None:
        samples = pd.DataFrame(
            {
                "external_id": ["SRR1", "SRR2", "SRR3"],
                "project": ["SRP1", "SRP1", "SRP2"],
                "organism": ["human", "human", "human"],
                "project_home": [
                    "data_sources/sra",
                    "data_sources/sra",
                    "data_sources/sra",
                ],
            }
        )
        with mock.patch(
            "recount3.search.available_samples", return_value=samples
        ):
            df = available_projects(organism="human")
        row_srp1 = df[df["project"] == "SRP1"].iloc[0]
        assert row_srp1["n_samples"] == 2

    def test_project_type_derived_from_project_home(self) -> None:
        samples = pd.DataFrame(
            {
                "external_id": ["SRR1"],
                "project": ["SRP1"],
                "organism": ["human"],
                "project_home": ["data_sources/sra"],
            }
        )
        with mock.patch(
            "recount3.search.available_samples", return_value=samples
        ):
            df = available_projects(organism="human")
        assert df["project_type"].iloc[0] == "data_sources"

    def test_no_key_cols_n_samples_undefined(self) -> None:
        """No key columns present -> n_samples column is populated with NA."""
        samples = pd.DataFrame(
            {
                "extra": ["foo"],
            }
        )
        with mock.patch(
            "recount3.search.available_samples", return_value=samples
        ):
            df = available_projects(organism="human")
        assert "n_samples" in df.columns


class TestProjectHomes:
    def test_empty_projects_returns_correct_columns(self) -> None:
        with mock.patch(
            "recount3.search.available_projects", return_value=pd.DataFrame()
        ):
            df = project_homes(organism="human", strict=False)
        assert "project_home" in df.columns
        assert "n_projects" in df.columns

    def test_project_home_column_absent_returns_empty(self) -> None:
        projects_df = pd.DataFrame({"project": ["SRP1"], "organism": ["human"]})
        with mock.patch(
            "recount3.search.available_projects", return_value=projects_df
        ):
            df = project_homes(organism="human")
        assert "project_home" in df.columns
        assert df.empty

    def test_n_projects_counts_correctly(self) -> None:
        projects_df = pd.DataFrame(
            {
                "project": ["SRP1", "SRP2"],
                "organism": ["human", "human"],
                "project_home": ["data_sources/sra", "data_sources/sra"],
                "project_type": ["data_sources", "data_sources"],
                "file_source": ["sra", "sra"],
            }
        )
        with mock.patch(
            "recount3.search.available_projects", return_value=projects_df
        ):
            df = project_homes(organism="human")
        assert df["n_projects"].iloc[0] == 2

    def test_project_type_derived_when_absent(self) -> None:
        projects_df = pd.DataFrame(
            {
                "project": ["SRP1"],
                "organism": ["human"],
                "project_home": ["data_sources/sra"],
                "file_source": ["sra"],
            }
        )
        with mock.patch(
            "recount3.search.available_projects", return_value=projects_df
        ):
            df = project_homes(organism="human")
        assert "project_type" in df.columns

    def test_file_source_added_as_na_when_absent(self) -> None:
        projects_df = pd.DataFrame(
            {
                "project": ["SRP1"],
                "organism": ["human"],
                "project_home": ["data_sources/sra"],
                "project_type": ["data_sources"],
            }
        )
        with mock.patch(
            "recount3.search.available_projects", return_value=projects_df
        ):
            df = project_homes(organism="human")
        assert "file_source" in df.columns


class TestSamplesForProject:
    def test_returns_sorted_sample_ids(self, cfg: Config) -> None:
        _seed_human_sra(cfg)
        with mock.patch("recount3.resource.default_config", return_value=cfg):
            ids = samples_for_project(
                organism="human", data_source="sra", project="SRP009615"
            )
        assert isinstance(ids, list)
        assert ids == sorted(ids)

    def test_no_ds_meta_raises(self) -> None:
        with mock.patch(
            "recount3.search.search_data_source_metadata", return_value=[]
        ):
            with pytest.raises(ValueError, match="No data-source metadata"):
                samples_for_project(
                    organism="human", data_source="sra", project="SRP1"
                )

    def test_empty_dataframe_raises(self) -> None:
        mock_res = mock.MagicMock()
        mock_res.load.return_value = pd.DataFrame()
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            with pytest.raises(ValueError, match="Empty"):
                samples_for_project(
                    organism="human", data_source="sra", project="SRP1"
                )

    def test_non_dataframe_raises(self) -> None:
        mock_res = mock.MagicMock()
        mock_res.load.return_value = "not a df"
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            with pytest.raises(ValueError, match="Empty"):
                samples_for_project(
                    organism="human", data_source="sra", project="SRP1"
                )

    def test_project_column_not_found_raises(self) -> None:
        df = pd.DataFrame({"unknown_col": ["x"]})
        mock_res = mock.MagicMock()
        mock_res.load.return_value = df
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            with pytest.raises(ValueError, match="Project column not found"):
                samples_for_project(
                    organism="human", data_source="sra", project="SRP1"
                )

    def test_project_not_in_metadata_raises(self) -> None:
        df = pd.DataFrame({"study": ["SRP_OTHER"], "external_id": ["SRR1"]})
        mock_res = mock.MagicMock()
        mock_res.load.return_value = df
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            with pytest.raises(ValueError, match="not found in metadata"):
                samples_for_project(
                    organism="human", data_source="sra", project="SRP1"
                )

    def test_no_sample_col_no_candidates_returns_empty(self) -> None:
        df = pd.DataFrame({"study": ["SRP1"], "other": ["x"]})
        mock_res = mock.MagicMock()
        mock_res.load.return_value = df
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            result = samples_for_project(
                organism="human", data_source="sra", project="SRP1"
            )
        assert result == []

    def test_no_sample_col_falls_back_to_candidate_col(self) -> None:
        """Column matching 'run' or 'sample' used as fallback samp_col."""
        df = pd.DataFrame({"study": ["SRP1"], "run_count": ["SRR999"]})
        mock_res = mock.MagicMock()
        mock_res.load.return_value = df
        with mock.patch(
            "recount3.search.search_data_source_metadata",
            return_value=[mock_res],
        ):
            result = samples_for_project(
                organism="human", data_source="sra", project="SRP1"
            )
        assert result == ["SRR999"]

    def test_md_prefix_stripped_before_column_lookup(self, cfg: Config) -> None:
        """Columns like 'sra.study' must be stripped before proj_keys lookup."""
        _seed_human_sra(cfg)
        with mock.patch("recount3.resource.default_config", return_value=cfg):
            ids = samples_for_project(
                organism="human", data_source="sra", project="SRP014565"
            )
        assert "SRR527082" in ids


class TestCreateSampleProjectLists:
    def test_empty_organism_loops_all(self, cfg: Config) -> None:
        _seed_human_sra(cfg)
        _seed_mouse_sra(cfg)
        with mock.patch("recount3.resource.default_config", return_value=cfg):
            samples, projects = create_sample_project_lists(organism="")
        assert isinstance(samples, list)
        assert isinstance(projects, list)

    def test_specific_organism_filters(self, cfg: Config) -> None:
        _seed_human_sra(cfg)
        with mock.patch("recount3.resource.default_config", return_value=cfg):
            samples, projects = create_sample_project_lists(organism="human")
        assert len(projects) > 0
        assert isinstance(projects[0], str)

    def test_sample_id_column_priority(self) -> None:
        """'external_id' should be picked before fallback columns."""
        samples_df = pd.DataFrame(
            {
                "external_id": ["SRR1", "SRR2"],
                "run": ["RUN1", "RUN2"],
                "project": ["SRP1", "SRP1"],
            }
        )
        projects_df = pd.DataFrame({"project": ["SRP1"]})
        with (
            mock.patch(
                "recount3.search.available_samples", return_value=samples_df
            ),
            mock.patch(
                "recount3.search.available_projects", return_value=projects_df
            ),
        ):
            samples, _ = create_sample_project_lists(organism="human")
        assert "SRR1" in samples
        assert "SRR2" in samples

    def test_no_samples_column_falls_back(self) -> None:
        """Falls back to 'run' when 'external_id' is absent."""
        samples_df = pd.DataFrame(
            {
                "run": ["RUN1"],
                "project": ["SRP1"],
            }
        )
        projects_df = pd.DataFrame({"project": ["SRP1"]})
        with (
            mock.patch(
                "recount3.search.available_samples", return_value=samples_df
            ),
            mock.patch(
                "recount3.search.available_projects", return_value=projects_df
            ),
        ):
            samples, _ = create_sample_project_lists(organism="human")
        assert "RUN1" in samples

    def test_no_recognized_sample_column_yields_empty_samples(self) -> None:
        """For-loop exhausts all priority columns without hit -> empty list."""
        samples_df = pd.DataFrame(
            {"mystery_col": ["SRR1"], "project": ["SRP1"]}
        )
        projects_df = pd.DataFrame({"project": ["SRP1"]})
        with (
            mock.patch(
                "recount3.search.available_samples", return_value=samples_df
            ),
            mock.patch(
                "recount3.search.available_projects", return_value=projects_df
            ),
        ):
            samples, _ = create_sample_project_lists(organism="human")
        assert samples == []

    def test_no_project_column_yields_empty_projects(self) -> None:
        samples_df = pd.DataFrame({"external_id": ["SRR1"], "other": ["x"]})
        projects_df = pd.DataFrame({"other": ["x"]})
        with (
            mock.patch(
                "recount3.search.available_samples", return_value=samples_df
            ),
            mock.patch(
                "recount3.search.available_projects", return_value=projects_df
            ),
        ):
            _, projects = create_sample_project_lists(organism="human")
        assert projects == []

    def test_empty_samples_and_projects(self) -> None:
        with (
            mock.patch(
                "recount3.search.available_samples", return_value=pd.DataFrame()
            ),
            mock.patch(
                "recount3.search.available_projects",
                return_value=pd.DataFrame(),
            ),
        ):
            samples, projects = create_sample_project_lists(organism="human")
        assert samples == []
        assert projects == []


class TestSearchProjectAll:
    def _patch_samples(self, sample_ids: list[str]):
        return mock.patch(
            "recount3.search.samples_for_project", return_value=sample_ids
        )

    def test_returns_resources_sorted(self) -> None:
        with self._patch_samples(["SRR387777"]):
            result = search_project_all(
                organism="human",
                data_source="sra",
                project="SRP009615",
                include_bigwig=False,
            )
        assert isinstance(result, list)
        assert len(result) > 0

    def test_includes_gene_and_exon_resources(self) -> None:
        with self._patch_samples([]):
            result = search_project_all(
                organism="human",
                data_source="sra",
                project="SRP009615",
                genomic_units=["gene", "exon"],
                include_bigwig=False,
                include_metadata=False,
                junction_extension=(),
            )
        resource_types = {r.description.resource_type for r in result}
        assert "count_files_gene_or_exon" in resource_types
        assert "annotations" in resource_types

    def test_no_genomic_units_skips_counts(self) -> None:
        with self._patch_samples([]):
            result = search_project_all(
                organism="human",
                data_source="sra",
                project="SRP009615",
                genomic_units=[],
                include_metadata=False,
                junction_extension=(),
                include_bigwig=False,
            )
        resource_types = {r.description.resource_type for r in result}
        assert "count_files_gene_or_exon" not in resource_types
        assert "annotations" not in resource_types

    def test_include_metadata_adds_metadata_resources(self) -> None:
        with self._patch_samples([]):
            result = search_project_all(
                organism="human",
                data_source="sra",
                project="SRP009615",
                include_metadata=True,
                genomic_units=[],
                junction_extension=(),
                include_bigwig=False,
            )
        resource_types = {r.description.resource_type for r in result}
        assert "metadata_files" in resource_types

    def test_exclude_metadata(self) -> None:
        with self._patch_samples([]):
            result = search_project_all(
                organism="human",
                data_source="sra",
                project="SRP009615",
                include_metadata=False,
                genomic_units=[],
                junction_extension=(),
                include_bigwig=False,
            )
        resource_types = {r.description.resource_type for r in result}
        assert "metadata_files" not in resource_types

    def test_include_bigwig_adds_bigwig_resources(self) -> None:
        with self._patch_samples(["SRR387777", "SRR387778"]):
            result = search_project_all(
                organism="human",
                data_source="sra",
                project="SRP009615",
                include_bigwig=True,
                include_metadata=False,
                genomic_units=[],
                junction_extension=(),
            )
        resource_types = {r.description.resource_type for r in result}
        assert "bigwig_files" in resource_types

    def test_include_bigwig_empty_samples_no_bigwig(self) -> None:
        with self._patch_samples([]):
            result = search_project_all(
                organism="human",
                data_source="sra",
                project="SRP009615",
                include_bigwig=True,
                include_metadata=False,
                genomic_units=[],
                junction_extension=(),
            )
        resource_types = {r.description.resource_type for r in result}
        assert "bigwig_files" not in resource_types

    def test_empty_junction_extension_skips_junctions(self) -> None:
        with self._patch_samples([]):
            result = search_project_all(
                organism="human",
                data_source="sra",
                project="SRP009615",
                junction_extension=(),
                include_metadata=False,
                genomic_units=[],
                include_bigwig=False,
            )
        resource_types = {r.description.resource_type for r in result}
        assert "count_files_junctions" not in resource_types

    def test_annotations_all_expands(self) -> None:
        with self._patch_samples([]):
            result = search_project_all(
                organism="human",
                data_source="sra",
                project="SRP009615",
                annotations="all",
                genomic_units=["gene"],
                include_metadata=False,
                junction_extension=(),
                include_bigwig=False,
            )
        ann_resources = [
            r for r in result if r.description.resource_type == "annotations"
        ]
        assert len(ann_resources) == len(_ANN_EXT_HUMAN)

    def test_invalid_project_raises(self) -> None:
        with mock.patch(
            "recount3.search.samples_for_project",
            side_effect=ValueError("Project 'NOPE' not found"),
        ):
            with pytest.raises(ValueError, match="not found"):
                search_project_all(
                    organism="human", data_source="sra", project="NOPE"
                )

    def test_result_is_deterministically_sorted(self) -> None:
        with self._patch_samples([]):
            r1 = search_project_all(
                organism="human",
                data_source="sra",
                project="SRP009615",
                include_bigwig=False,
            )
            r2 = search_project_all(
                organism="human",
                data_source="sra",
                project="SRP009615",
                include_bigwig=False,
            )
        assert [r.url for r in r1] == [r.url for r in r2]

    def test_gene_only_genomic_unit(self) -> None:
        with self._patch_samples([]):
            result = search_project_all(
                organism="human",
                data_source="sra",
                project="SRP009615",
                genomic_units=["gene"],
                include_metadata=False,
                junction_extension=(),
                include_bigwig=False,
            )
        count_resources = [
            r
            for r in result
            if r.description.resource_type == "count_files_gene_or_exon"
        ]
        genomic_units_found = {
            getattr(r.description, "genomic_unit", None)
            for r in count_resources
        }
        assert "gene" in genomic_units_found
        assert "exon" not in genomic_units_found

    def test_exon_only_genomic_unit(self) -> None:
        with self._patch_samples([]):
            result = search_project_all(
                organism="human",
                data_source="sra",
                project="SRP009615",
                genomic_units=["exon"],
                include_metadata=False,
                junction_extension=(),
                include_bigwig=False,
            )
        count_resources = [
            r
            for r in result
            if r.description.resource_type == "count_files_gene_or_exon"
        ]
        genomic_units_found = {
            getattr(r.description, "genomic_unit", None)
            for r in count_resources
        }
        assert "exon" in genomic_units_found
        assert "gene" not in genomic_units_found


class TestAnnotationConstants:
    """Sanity-check the module-level annotation constants."""

    def test_human_extensions_non_empty(self) -> None:
        assert len(_ANN_EXT_HUMAN) > 0

    def test_mouse_extensions_non_empty(self) -> None:
        assert len(_ANN_EXT_MOUSE) > 0

    def test_human_name_to_ext_values_in_human_exts(self) -> None:
        for ext in _ANNOTATION_NAME_TO_EXT_HUMAN.values():
            assert ext in _ANN_EXT_HUMAN

    def test_mouse_name_to_ext_values_in_mouse_exts(self) -> None:
        for ext in _ANNOTATION_NAME_TO_EXT_MOUSE.values():
            assert ext in _ANN_EXT_MOUSE
