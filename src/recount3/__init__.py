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
      organism="human", genomic_unit="gene", annotation_file_extension="G026")
  res = R3Resource(desc)
  res.download(path=None, cache_mode="enable")
  df = res.load()
"""

from __future__ import annotations

from .version import __version__

from .config import Config, default_config
from .errors import (
    CompatibilityError,
    ConfigurationError,
    DownloadError,
    LoadError,
    Recount3Error,
)
from .types import CacheMode, CompatibilityMode, FieldSpec, StringOrIterable

from ._descriptions import (
    R3Annotations,
    R3BigWig,
    R3DataSourceMetadata,
    R3DataSources,
    R3GeneOrExonCounts,
    R3JunctionCounts,
    R3ProjectMetadata,
    R3ResourceDescription,
)

from .resource import R3Resource
from ._bigwig import BigWigFile
from .bundle import R3ResourceBundle
from .se import (  # noqa: F401
    build_summarized_experiment,
    build_ranged_summarized_experiment,
)
from .search import (
    samples_for_project,
    search_project_all,
    create_sample_project_lists,
    search_annotations,
    search_bigwig_files,
    search_count_files_gene_or_exon,
    search_count_files_junctions,
    search_data_source_metadata,
    search_data_sources,
    search_metadata_files,
    annotation_ext,
    annotation_options,

)

__all__ = [
    "__version__",
    # Config
    "Config",
    "default_config",
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
