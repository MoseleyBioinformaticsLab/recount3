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
#
"""Typed Python API for the recount3 RNA-seq data repository.

`recount3 <https://rna.recount.bio/docs/>`_ is a large-scale uniformly
processed RNA-seq resource covering tens of thousands of human and mouse
samples across multiple data sources (SRA, GTEx, TCGA). This package
provides a small, typed Python API for discovering, downloading, and loading
recount3 data files.

The public surface is intentionally flat: core classes and helper functions
are re-exported here for discovery and convenience. The one exception is the
:mod:`recount3.se` submodule, which is exposed as ``r3.se`` rather than
flattened. Its normalization and conversion helpers operate on BiocPy
experiments, so they are called through that prefix, for example
``r3.se.compute_tpm(rse)``.

Typical usage example: high-level (recommended for most cases)::

  import recount3 as r3

  rse = r3.create_rse(
      project="SRP009615",
      organism="human",
      annotation_label="gencode_v26",
  )

Typical usage example: lower-level (multi-project or custom workflows)::

  import recount3 as r3

  bundle = r3.R3ResourceBundle.discover(
      organism="human",
      data_source="sra",
      project=["SRP009615", "SRP001558"],
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
      Install with ``pip install "recount3[biocpy]"``.
    * ``pyBigWig``: required for BigWig file access via
      :class:`~recount3._bigwig.BigWigFile`. Install with
      ``pip install "recount3[bigwig]"``.
"""

from __future__ import annotations

from recount3 import se as se
from recount3._bigwig import BigWigFile
from recount3._descriptions import R3Annotations
from recount3._descriptions import R3BigWig
from recount3._descriptions import R3CollectionMetadata
from recount3._descriptions import R3DataSourceMetadata
from recount3._descriptions import R3DataSources
from recount3._descriptions import R3GeneOrExonCounts
from recount3._descriptions import R3JunctionCounts
from recount3._descriptions import R3ProjectMetadata
from recount3._descriptions import R3ResourceDescription
from recount3.bundle import R3ResourceBundle
from recount3.config import Config
from recount3.config import default_config
from recount3.config import recount3_cache
from recount3.config import recount3_cache_files
from recount3.config import recount3_cache_rm
from recount3.errors import CompatibilityError
from recount3.errors import ConfigurationError
from recount3.errors import DownloadError
from recount3.errors import LoadError
from recount3.errors import MissingRangesError
from recount3.errors import RangesCoverageError
from recount3.errors import RangesError
from recount3.errors import Recount3Error
from recount3.resource import R3Resource
from recount3.se import build_ranged_summarized_experiment
from recount3.se import build_summarized_experiment
from recount3.se import create_rse
from recount3.search import annotation_ext
from recount3.search import annotation_label
from recount3.search import annotation_options
from recount3.search import available_projects
from recount3.search import available_samples
from recount3.search import create_sample_project_lists
from recount3.search import project_homes
from recount3.search import samples_for_project
from recount3.search import search_annotations
from recount3.search import search_bigwig_files
from recount3.search import search_collection
from recount3.search import search_collections
from recount3.search import search_count_files_gene_or_exon
from recount3.search import search_count_files_junctions
from recount3.search import search_data_source_metadata
from recount3.search import search_data_sources
from recount3.search import search_metadata_files
from recount3.search import search_project_all
from recount3.types import CacheMode
from recount3.types import CompatibilityMode
from recount3.types import FieldSpec
from recount3.types import StringOrIterable
from recount3.types import StrPath
from recount3.version import __version__

__all__ = [  # noqa: RUF022 - Exports are grouped by API role.
    "__version__",
    # Submodules
    "se",
    # Config
    "Config",
    "default_config",
    "recount3_cache",
    "recount3_cache_rm",
    "recount3_cache_files",
    # Errors
    "LoadError",
    "RangesError",
    "DownloadError",
    "Recount3Error",
    "CompatibilityError",
    "ConfigurationError",
    "MissingRangesError",
    "RangesCoverageError",
    # Types
    "StrPath",
    "CacheMode",
    "FieldSpec",
    "StringOrIterable",
    "CompatibilityMode",
    # Descriptions
    "R3BigWig",
    "R3Annotations",
    "R3DataSources",
    "R3JunctionCounts",
    "R3ProjectMetadata",
    "R3GeneOrExonCounts",
    "R3CollectionMetadata",
    "R3DataSourceMetadata",
    "R3ResourceDescription",
    # Core
    "BigWigFile",
    "R3Resource",
    "R3ResourceBundle",
    # Builders
    "create_rse",
    "build_summarized_experiment",
    "build_ranged_summarized_experiment",
    # Search
    "project_homes",
    "annotation_ext",
    "annotation_label",
    "available_samples",
    "search_collection",
    "annotation_options",
    "available_projects",
    "search_annotations",
    "search_collections",
    "search_project_all",
    "samples_for_project",
    "search_bigwig_files",
    "search_data_sources",
    "search_metadata_files",
    "create_sample_project_lists",
    "search_data_source_metadata",
    "search_count_files_junctions",
    "search_count_files_gene_or_exon",
]
