"""Optional BiocPy dependency handling.

Lazy, cached imports for BiocPy-backed classes.
"""

from __future__ import annotations

import importlib
import types
from typing import TYPE_CHECKING

from recount3 import errors

if TYPE_CHECKING:
    import biocframe  # type: ignore[import-not-found]
    import genomicranges  # type: ignore[import-not-found]
    import summarizedexperiment  # type: ignore[import-not-found]


def _format_biocpy_install_instructions(module_name: str) -> str:
    """Formats an error message with BiocPy installation instructions."""
    return (
        f"Optional dependency '{module_name}' is required when "
        "using this feature. Install the BiocPy packages, "
        "for example:\n\n"
        "  pip install biocframe genomicranges summarizedexperiment\n"
    )


def _import_optional_module(module_name: str) -> types.ModuleType:
    """Imports an optional module or raises a CompatibilityError.

    Args:
        module_name: The module name to import.

    Returns:
        The imported module.

    Raises:
        errors.CompatibilityError: If the module is not installed.
    """
    try:
        return importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        raise errors.CompatibilityError(
            _format_biocpy_install_instructions(module_name),
        ) from exc


def get_biocframe_class() -> type["biocframe.BiocFrame"]:
    """Returns the BiocFrame class (biocframe.BiocFrame)."""
    module = _import_optional_module("biocframe")
    return module.BiocFrame


def get_genomicranges_class() -> type["genomicranges.GenomicRanges"]:
    """Returns the GenomicRanges class (genomicranges.GenomicRanges)."""
    module = _import_optional_module("genomicranges")
    return module.GenomicRanges


def get_summarizedexperiment_class() -> type["summarizedexperiment.SummarizedExperiment"]:
    """Returns the SummarizedExperiment class."""
    module = _import_optional_module("summarizedexperiment")
    return module.SummarizedExperiment


def get_ranged_summarizedexperiment_class(
) -> type["summarizedexperiment.RangedSummarizedExperiment"]:
    """Returns the RangedSummarizedExperiment class."""
    module = _import_optional_module("summarizedexperiment")
    return module.RangedSummarizedExperiment
