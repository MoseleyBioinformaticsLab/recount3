"""Typed Python API for the recount3 RNA-seq data repository.

`recount3 <https://rna.recount.bio/docs/>`_ is a large-scale uniformly
processed RNA-seq resource covering tens of thousands of human and mouse
samples across multiple data sources (SRA, GTEx, TCGA). This package
provides a small, typed Python API for discovering, downloading, and loading
recount3 data files.

The public surface is intentionally flat: core classes and helper functions
are re-exported here for discovery and convenience.

Typical usage example: high-level (recommended for most cases)::

  from recount3 import create_rse

  rse = create_rse(
      project="SRP009615",
      organism="human",
      annotation_label="gencode_v26",
  )

Typical usage example: lower-level (multi-project or custom workflows)::

  from recount3 import R3ResourceBundle

  bundle = R3ResourceBundle.discover(
      organism="human",
      data_source="sra",
      project=["SRP009615", "SRP012682"],
  )
  counts = bundle.filter(
      resource_type="count_files_gene_or_exon",
      genomic_unit="gene",
  ).stack_count_matrices()

Note:
    Several features require optional dependencies:

    * BiocPy (``biocframe``, ``summarizedexperiment``, ``genomicranges``):
      required for :func:`create_rse`,
      :meth:`~recount3.bundle.R3ResourceBundle.to_ranged_summarized_experiment`,
      and the compute utilities in :mod:`recount3.se`.
      Install with ``pip install biocframe summarizedexperiment genomicranges``.
    * pyBigWig: required for BigWig file access via
      :class:`~recount3._bigwig.BigWigFile`. Install with
      ``pip install pyBigWig``.
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
from recount3.se import (
    build_ranged_summarized_experiment,
    build_summarized_experiment,
    create_rse,
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
    "create_rse",
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
