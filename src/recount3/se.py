"""High-level builders and utilities for SummarizedExperiment objects.

This module provides user-facing helpers for constructing BiocPy
:sup:`SummarizedExperiment` and
:sup:`RangedSummarizedExperiment` objects from recount3 projects, plus
utilities for working with SRA-style sample attributes.

The heavy lifting is implemented as methods on
:class:`recount3.bundle.R3ResourceBundle`:

* :meth:`recount3.bundle.R3ResourceBundle.to_summarized_experiment`
* :meth:`recount3.bundle.R3ResourceBundle.to_ranged_summarized_experiment`

The functions in this module are thin wrappers around those methods (for
example, :func:`create_rse`) and convenience utilities that operate on
BiocPy objects directly (for example, :func:`expand_sra_attributes`).

Typical usage example:

  from recount3 import se

  # High-level one-shot helper, mirroring the R create_rse() API.
  rse = se.create_rse(
      project="SRP009615",
      genomic_unit="gene",
      organism="human",
      data_source="sra",
  )
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import logging
import pandas as pd

from recount3.bundle import R3ResourceBundle
from recount3.search import annotation_ext

# Optional BiocPy imports, used only for expand_sra_attributes().
try:  # pragma: no cover - optional, exercised only when installed.
    from biocframe import BiocFrame  # type: ignore[import-not-found]
    from summarizedexperiment import (  # type: ignore[import-not-found]
        RangedSummarizedExperiment,
        SummarizedExperiment,
    )
except Exception:  # pragma: no cover
    BiocFrame = None
    SummarizedExperiment = None
    RangedSummarizedExperiment = None


def _require_biocpy() -> None:
    """Ensure BiocPy packages are importable for SE/RSE utilities.

    Raises:
      ImportError: If any required BiocPy package is missing.
    """
    if any(
        x is None
        for x in (BiocFrame, SummarizedExperiment, RangedSummarizedExperiment)
    ):
        raise ImportError(
            "This feature requires BiocPy packages. Install with:\n"
            "  pip install summarizedexperiment biocframe genomicranges\n"
            "or via conda (conda-forge)."
        )


def _expand_sra_attributes_df(
    col_df: pd.DataFrame,
    *,
    column: str = "sra.sample_attributes",
    prefix: str = "sra_attribute.",
) -> pd.DataFrame:
    """Return a copy of column metadata with SRA attributes expanded.

    recount3 stores SRA sample attributes in a single string column such
    as ``'sra.sample_attributes'`` using the encoding::

        "age;;67.78|biomaterial_provider;;LIBD|disease;;Control|..."

    This helper mirrors the behavior of the recount3 R function
    ``expand_sra_attributes()``: each entry is split on ``'|'`` into
    individual attributes, then on ``';;'`` into ``key`` and ``value``
    pairs. New columns named ``{prefix}{key}`` (with spaces in ``key``
    replaced by ``'_'``) are appended to the metadata.

    Args:
      col_df: Column metadata DataFrame.
      column: Name of the column carrying the raw SRA attribute strings.
        If this column is missing, ``col_df`` is returned unchanged.
      prefix: Prefix to prepend to attribute names when forming new
        columns.

    Returns:
      A new DataFrame with the same index as ``col_df`` and additional
      columns containing one parsed SRA attribute per column.
    """
    if column not in col_df.columns:
        return col_df.copy()

    ser = col_df[column]
    rows: list[dict[str, Any]] = []

    for value in ser:
        row_map: dict[str, Any] = {}
        if isinstance(value, str) and value:
            parts = value.split("|")
            for entry in parts:
                if not entry or ";;" not in entry:
                    continue
                key, val = entry.split(";;", 1)
                key_clean = key.strip()
                if not key_clean:
                    continue
                col_name = f"{prefix}{key_clean.replace(' ', '_')}"
                row_map[col_name] = val.strip()
        rows.append(row_map)

    if not rows:
        return col_df.copy()

    attrs_df = pd.DataFrame(rows, index=col_df.index)
    merged = pd.concat([col_df, attrs_df], axis=1)
    return merged


def _normalize_genomic_unit(genomic_unit: str) -> str:
    """Return a normalized genomic unit string and validate it.

    Args:
      genomic_unit: Requested feature level.

    Returns:
      Lowercase genomic unit string.

    Raises:
      ValueError: If the genomic unit is not one of ``"gene"``,
        ``"exon"``, or ``"junction"``.
    """
    gu = str(genomic_unit).strip().lower()
    valid = {"gene", "exon", "junction"}
    if gu not in valid:
        raise ValueError(
            f"Invalid genomic_unit {genomic_unit!r}; expected one of "
            f"{sorted(valid)!r}."
        )
    return gu


def _resolve_annotation_extension(
    *,
    organism: str,
    genomic_unit: str,
    annotation: str | None,
    annotation_file_extension: str | None,
) -> str | None:
    """Resolve a single annotation extension for create_rse-style helpers.

    For gene and exon assays, this helper maps a human-readable
    annotation name (for example, ``"gencode_v26"``) or an explicit
    extension (for example, ``"G026"``) to the underlying recount3
    extension. When both ``annotation`` and
    ``annotation_file_extension`` are provided, the explicit extension
    wins.

    Args:
      organism: Organism identifier (``"human"`` or ``"mouse"``).
      genomic_unit: Normalized genomic unit (``"gene"``, ``"exon"``, or
        ``"junction"``).
      annotation: Optional human-readable annotation label.
      annotation_file_extension: Optional explicit annotation extension.

    Returns:
      A single annotation extension string for gene/exon assays, or
      :data:`None` when no annotation is required (for example,
      junction-level assays).
    """
    if genomic_unit not in {"gene", "exon"}:
        return None

    if annotation_file_extension:
        return str(annotation_file_extension).strip()

    if annotation:
        return annotation_ext(organism=organism, annotation=annotation)

    org = organism.strip().lower()
    if org == "human":
        return "G026"
    if org == "mouse":
        return "M023"

    raise ValueError(
        f"Unsupported organism {organism!r}; expected 'human' or 'mouse'."
    )


def build_summarized_experiment(
    bundle: R3ResourceBundle,
    *,
    genomic_unit: str,
    annotation_file_extension: str | None = None,
    assay_name: str = "counts",
    join: str = "inner",
    autoload: bool = True,
) -> SummarizedExperiment:
    """Create a SummarizedExperiment from a resource bundle.

    This is a convenience wrapper around
    :meth:`recount3.bundle.R3ResourceBundle.to_summarized_experiment`
    that preserves the original functional API of the early Python
    recount3 prototypes.

    Args:
      bundle: Resource bundle containing counts and metadata.
      genomic_unit: One of ``"gene"``, ``"exon"``, or ``"junction"``.
      annotation_file_extension: Optional annotation code for gene/exon
        assays (for example, ``"G026"``).
      assay_name: Name for the count assay in the SE.
      join: Join policy across projects (pandas concatenation join).
      autoload: If :data:`True`, load resources transparently.

    Returns:
      A :class:`summarizedexperiment.SummarizedExperiment` instance.
    """
    unit = _normalize_genomic_unit(genomic_unit)
    return bundle.to_summarized_experiment(
        genomic_unit=unit,
        annotation_file_extension=annotation_file_extension,
        assay_name=assay_name,
        join=join,
        autoload=autoload,
    )


def build_ranged_summarized_experiment(
    bundle: R3ResourceBundle,
    *,
    genomic_unit: str,
    annotation_file_extension: str | None = None,
    junction_rr_preferred: bool = True,
    assay_name: str = "counts",
    join: str = "inner",
    autoload: bool = True,
    allow_fallback_to_se: bool = False,
) -> RangedSummarizedExperiment | SummarizedExperiment:
    """Create a RangedSummarizedExperiment when ranges can be resolved.

    This is a convenience wrapper around
    :meth:`recount3.bundle.R3ResourceBundle.to_ranged_summarized_experiment`.

    Args:
      bundle: Resources for counts, metadata, and (optionally)
        annotations.
      genomic_unit: One of ``"gene"``, ``"exon"``, or ``"junction"``.
      annotation_file_extension: Annotation code for gene/exon, if
        desired.
      junction_rr_preferred: Prefer RR for junction coordinates when
        present.
      assay_name: Name for the count assay in the output.
      join: Join policy across projects when stacking.
      autoload: If :data:`True`, load resources transparently.
      allow_fallback_to_se: If :data:`True`, return a plain SE when
        ranges are unavailable.

    Returns:
      A :class:`summarizedexperiment.RangedSummarizedExperiment` object,
      or a plain :class:`summarizedexperiment.SummarizedExperiment` when
      ``allow_fallback_to_se`` is :data:`True` and ranges cannot be
      resolved.
    """
    unit = _normalize_genomic_unit(genomic_unit)
    return bundle.to_ranged_summarized_experiment(
        genomic_unit=unit,
        annotation_file_extension=annotation_file_extension,
        junction_rr_preferred=junction_rr_preferred,
        assay_name=assay_name,
        join=join,
        autoload=autoload,
        allow_fallback_to_se=allow_fallback_to_se,
    )


def create_ranged_summarized_experiment(
    *,
    project: str,
    genomic_unit: str = "gene",
    organism: str = "human",
    data_source: str = "sra",
    annotation: str | None = None,
    annotation_file_extension: str | None = None,
    junction_type: str = "ALL",
    junction_file_extensions: Sequence[str] | None = None,
    include_metadata: bool = True,
    include_bigwig: bool = False,
    join: str = "inner",
    autoload: bool = True,
    allow_fallback_to_se: bool = False,
    strict: bool = True,
) -> RangedSummarizedExperiment | SummarizedExperiment:
    """High-level helper that mirrors recount3's ``create_rse()`` in R.

    This function hides the intermediate bundle construction step by:

      1. Discovering all resources for a project via
         :meth:`recount3.bundle.R3ResourceBundle.discover`.
      2. Stacking expression matrices across projects and samples.
      3. Resolving genomic ranges (GTF for gene/exon; RR for junctions).
      4. Building a BiocPy :class:`RangedSummarizedExperiment` (or SE
         fallback) using the bundle methods.

    Args:
      project: Study or project identifier (for example, ``"SRP009615"``).
      genomic_unit: One of ``{"gene", "exon", "junction"}``.
      organism: Organism identifier (``"human"`` or ``"mouse"``).
      data_source: Data source (``"sra"``, ``"gtex"``, or ``"tcga"``).
      annotation: Optional human-readable annotation label for gene/exon
        assays, such as ``"gencode_v26"`` or ``"gencode_v29"``. Ignored
        for junction-level assays.
      annotation_file_extension: Optional explicit annotation extension
        (for example, ``"G026"``). When provided, this takes precedence
        over ``annotation``. Ignored for junction-level assays.
      junction_type: Junction type; typically ``"ALL"``.
      junction_file_extensions: Iterable of junction artifact extensions
        to include (for example, ``("MM",)`` or ``("MM", "RR")``). If
        :data:`None`, the default ``("MM",)`` is used.
      include_metadata: Whether to include the five project metadata
        tables in the underlying bundle (recommended).
      include_bigwig: Whether to include per-sample BigWig coverage
        resources in the bundle. These can be large.
      join: Join policy across projects when stacking count matrices
        (passed to :func:`build_ranged_summarized_experiment`).
      autoload: If :data:`True`, load resources transparently.
      allow_fallback_to_se: If :data:`True`, construct a plain SE when
        genomic ranges cannot be derived for the requested combination.
      strict: If :data:`True`, propagate validation errors from the
        search layer (for example, missing projects or incompatible
        combinations).

    Returns:
      A :class:`summarizedexperiment.RangedSummarizedExperiment` object,
      or a plain :class:`summarizedexperiment.SummarizedExperiment` when
      ``allow_fallback_to_se`` is :data:`True` and ranges cannot be
      resolved.

    Raises:
      ImportError: If BiocPy packages are not installed.
      ValueError: If inputs are invalid or no counts are found.
      TypeError: If the underlying
        :class:`RangedSummarizedExperiment` constructor rejects all
        compatibility variants.
    """
    unit = _normalize_genomic_unit(genomic_unit)
    ann_ext = _resolve_annotation_extension(
        organism=organism,
        genomic_unit=unit,
        annotation=annotation,
        annotation_file_extension=annotation_file_extension,
    )

    if junction_file_extensions is None:
        junction_exts: tuple[str, ...] = ("MM",)
    else:
        junction_exts = tuple(str(ext).strip() for ext in junction_file_extensions)

    project_bundle = R3ResourceBundle.discover(
        organism=organism,
        data_source=data_source,
        project=project,
        genomic_units=(unit,),
        annotations="default" if ann_ext is None else (ann_ext,),
        junction_exts=junction_exts,
        junction_type=junction_type,
        include_metadata=include_metadata,
        include_bigwig=include_bigwig,
        strict=strict,
        deduplicate=True,
    )

    return project_bundle.to_ranged_summarized_experiment(
        genomic_unit=unit,
        annotation_file_extension=ann_ext,
        junction_rr_preferred=True,
        assay_name="counts",
        join=join,
        autoload=autoload,
        allow_fallback_to_se=allow_fallback_to_se,
    )


def create_rse(
    *,
    project: str,
    genomic_unit: str = "gene",
    organism: str = "human",
    data_source: str = "sra",
    annotation: str | None = None,
    annotation_file_extension: str | None = None,
    junction_type: str = "ALL",
    junction_file_extensions: Sequence[str] | None = None,
    include_metadata: bool = True,
    include_bigwig: bool = False,
    join: str = "inner",
    autoload: bool = True,
    allow_fallback_to_se: bool = False,
    strict: bool = True,
) -> RangedSummarizedExperiment | SummarizedExperiment:
    """Alias for :func:`create_ranged_summarized_experiment`.

    The parameters and behavior are identical; see that function for
    full documentation.
    """
    return create_ranged_summarized_experiment(
        project=project,
        genomic_unit=genomic_unit,
        organism=organism,
        data_source=data_source,
        annotation=annotation,
        annotation_file_extension=annotation_file_extension,
        junction_type=junction_type,
        junction_file_extensions=junction_file_extensions,
        include_metadata=include_metadata,
        include_bigwig=include_bigwig,
        join=join,
        autoload=autoload,
        allow_fallback_to_se=allow_fallback_to_se,
        strict=strict,
    )


def expand_sra_attributes(
    obj: (
        RangedSummarizedExperiment
        | SummarizedExperiment
        | pd.DataFrame
    ),
    *,
    column: str = "sra.sample_attributes",
    prefix: str = "sra_attribute.",
) -> (
    RangedSummarizedExperiment
    | SummarizedExperiment
    | pd.DataFrame
):
    """Expand encoded SRA sample attributes into separate columns.

    This function mirrors the recount3 R helper
    ``expand_sra_attributes()``.

    It understands the SRA encoding used by recount3, where a single
    column (typically ``'sra.sample_attributes'``) contains strings of
    the form::

        "age;;67.78|biomaterial_provider;;LIBD|disease;;Control|..."

    Each ``key;;value`` pair becomes a new column named
    ``{prefix}{key}`` (with spaces in ``key`` replaced by ``'_'``), and
    the parsed values are stored per sample. The original string column
    is preserved.

    The function supports two calling styles:

    * Passing a :class:`pandas.DataFrame` of column metadata, in which
      case a new DataFrame is returned with extra columns.
    * Passing a BiocPy
      :class:`summarizedexperiment.SummarizedExperiment` or
      :class:`summarizedexperiment.RangedSummarizedExperiment`, in which
      case a new object of the same class is returned with updated
      ``column_data``.

    Args:
      obj: Column metadata DataFrame or a BiocPy SE/RSE object.
      column: Name of the column that holds the raw SRA attribute
        strings.
      prefix: Prefix to prepend to the generated attribute column names.

    Returns:
      A new object of the same type as ``obj`` with additional columns
      corresponding to parsed SRA attributes. If the requested column is
      missing, ``obj`` is returned unchanged.

    Raises:
      ImportError: If a SummarizedExperiment/RangedSummarizedExperiment
        is supplied but BiocPy packages are not installed.
      AttributeError: If the BiocPy object does not expose
        ``column_data`` or ``set_column_data`` in the expected API.
      TypeError: If ``obj`` is neither a pandas DataFrame nor a
        supported BiocPy experiment object.
    """
    # DataFrame mode: no BiocPy dependency.
    if isinstance(obj, pd.DataFrame):
        return _expand_sra_attributes_df(
            obj,
            column=column,
            prefix=prefix,
        )

    _require_biocpy()

    if not isinstance(obj, (SummarizedExperiment, RangedSummarizedExperiment)):
        raise TypeError(
            "expand_sra_attributes() expects a pandas DataFrame or a BiocPy "
            "SummarizedExperiment/RangedSummarizedExperiment; got "
            f"{type(obj)!r}."
        )

    col_bf = obj.column_data
    try:
        col_df = col_bf.to_pandas()
    except AttributeError as exc:
        raise AttributeError(
            "Expected column_data to implement to_pandas(); got "
            f"{type(col_bf)!r} instead."
        ) from exc

    if column not in col_df.columns:
        logging.warning(
            "expand_sra_attributes: column %r not present in column_data; "
            "returning object unchanged.",
            column,
        )
        return obj

    expanded_df = _expand_sra_attributes_df(
        col_df,
        column=column,
        prefix=prefix,
    )
    expanded_bf = BiocFrame.from_pandas(expanded_df)
    return obj.set_column_data(expanded_bf)
