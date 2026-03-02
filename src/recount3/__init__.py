"""Typed, minimal API for recount3 data.

This package provides a small, typed API for interacting with the recount3
repository. The public surface is intentionally flat: core classes and helper
functions are re-exported here for discovery and convenience.

Typical usage example:

  from recount3 import (
      R3Resource, R3ResourceDescription, R3Annotations,
      search_count_files_gene_or_exon, R3ResourceBundle,
  )

  desc = R3Annotations(
      organism="human", genomic_unit="gene", annotation_extension="G026")
  res = R3Resource(desc)
  res.download(path=None, cache_mode="enable")
  df = res.load()
"""

from __future__ import annotations

from recount3._bigwig import BigWigFile
from recount3._descriptions import (
    R3Annotations,
    R3BigWig,
    R3DataSourceMetadata,
    R3DataSources,
    R3GeneOrExonCounts,
    R3JunctionCounts,
    R3ProjectMetadata,
    R3ResourceDescription,
)
from recount3.bundle import R3ResourceBundle
from recount3.config import (
    Config,
    default_config,
    recount3_cache,
    recount3_cache_files,
    recount3_cache_rm,
)
from recount3.errors import (
    CompatibilityError,
    ConfigurationError,
    DownloadError,
    LoadError,
    Recount3Error,
)
from recount3.resource import R3Resource
from recount3.se import (  # noqa: F401
    build_ranged_summarized_experiment,
    build_summarized_experiment,
)
from recount3.search import (
    annotation_ext,
    annotation_options,
    create_sample_project_lists,
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
from recount3.types import (
    CacheMode,
    CompatibilityMode,
    FieldSpec,
    StringOrIterable,
)
from recount3.version import __version__

__all__ = [
    "__version__",
    # Config
    "Config",
    "default_config",
    "recount3_cache",
    "recount3_cache_files",
    "recount3_cache_rm",
    # Errors
    "Recount3Error",
    "ConfigurationError",
    "DownloadError",
    "LoadError",
    "CompatibilityError",
    # Types
    "CacheMode",
    "CompatibilityMode",
    "FieldSpec",
    "StringOrIterable",
    # Descriptions
    "R3ResourceDescription",
    "R3Annotations",
    "R3GeneOrExonCounts",
    "R3JunctionCounts",
    "R3ProjectMetadata",
    "R3BigWig",
    "R3DataSources",
    "R3DataSourceMetadata",
    # Core
    "R3Resource",
    "BigWigFile",
    "R3ResourceBundle",
    # Builders
    "build_ranged_summarized_experiment",
    "build_summarized_experiment",
    # Search
    "samples_for_project",
    "search_project_all",
    "search_annotations",
    "search_count_files_gene_or_exon",
    "search_count_files_junctions",
    "search_metadata_files",
    "search_bigwig_files",
    "search_data_sources",
    "search_data_source_metadata",
    "create_sample_project_lists",
    "annotation_ext",
    "annotation_options",
]
