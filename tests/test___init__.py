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
#   may be used to endorse or promote products derived from this software without
#   specific prior written permission.
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
