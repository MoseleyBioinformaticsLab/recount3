from __future__ import annotations

import recount3


def test_package_imports_without_error() -> None:
    """Importing recount3 must not raise."""
    assert recount3 is not None


def test_version_accessible_at_top_level() -> None:
    """__version__ must be directly accessible from the package."""
    assert isinstance(recount3.__version__, str)
    assert recount3.__version__


def test_all_declared_exports_are_accessible() -> None:
    """Every name listed in __all__ must resolve on the package."""
    missing = [name for name in recount3.__all__ if not hasattr(recount3, name)]
    assert missing == [], f"Names in __all__ but not accessible: {missing}"


def test_all_is_a_non_empty_list() -> None:
    """__all__ must be a non-empty list of strings."""
    assert isinstance(recount3.__all__, list)
    assert len(recount3.__all__) > 0
    assert all(isinstance(name, str) for name in recount3.__all__)


def test_core_classes_are_exported() -> None:
    """Spot-check that the core public classes are present and are types."""
    for name in (
        "R3Resource",
        "R3ResourceDescription",
        "R3ResourceBundle",
        "BigWigFile",
        "R3Annotations",
        "R3GeneOrExonCounts",
        "R3JunctionCounts",
        "R3ProjectMetadata",
        "R3BigWig",
        "R3DataSources",
        "R3DataSourceMetadata",
    ):
        obj = getattr(recount3, name)
        assert isinstance(obj, type), f"{name} should be a class"


def test_search_functions_are_callable() -> None:
    """Spot-check that search helpers are callable."""
    for name in (
        "search_project_all",
        "search_annotations",
        "search_count_files_gene_or_exon",
        "search_count_files_junctions",
        "search_metadata_files",
        "search_bigwig_files",
        "search_data_sources",
        "search_data_source_metadata",
        "samples_for_project",
        "create_sample_project_lists",
        "annotation_ext",
        "annotation_options",
    ):
        assert callable(getattr(recount3, name)), f"{name} should be callable"


def test_error_classes_are_exported() -> None:
    """Error classes must be present and be exception subclasses."""
    for name in (
        "Recount3Error",
        "ConfigurationError",
        "DownloadError",
        "LoadError",
        "CompatibilityError",
    ):
        obj = getattr(recount3, name)
        assert isinstance(obj, type)
        assert issubclass(obj, Exception), f"{name} should be an Exception subclass"


def test_config_helpers_are_exported() -> None:
    """Config class and cache helpers must be callable."""
    assert isinstance(recount3.Config, type)
    for name in ("default_config", "recount3_cache", "recount3_cache_files", "recount3_cache_rm"):
        assert callable(getattr(recount3, name)), f"{name} should be callable"


def test_se_builders_are_callable() -> None:
    """SummarizedExperiment builders must be callable."""
    assert callable(recount3.build_summarized_experiment)
    assert callable(recount3.build_ranged_summarized_experiment)
    assert callable(recount3.create_rse)
