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

  rse = se.create_rse(
      project="SRP009615",
      genomic_unit="gene",
      organism="human",
      data_source="sra",
  )
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, TYPE_CHECKING

import logging
import numbers
import numpy as np
import pandas as pd

from recount3.bundle import R3ResourceBundle
from recount3.search import annotation_ext

if TYPE_CHECKING:  # for static type checkers
    from biocframe import BiocFrame
    from summarizedexperiment import (
        RangedSummarizedExperiment,
        SummarizedExperiment,
    )

# Runtime caches for the optional BiocPy classes
_BiocFrame: type[Any] | None = None
_SummarizedExperiment: type[Any] | None = None
_RangedSummarizedExperiment: type[Any] | None = None

def _require_biocpy() -> tuple[type[Any], type[Any], type[Any]]:
    """Ensure BiocPy packages are importable for SE/RSE utilities.

    Returns:
        (BiocFrame, SummarizedExperiment, RangedSummarizedExperiment) classes.

    Raises:
      ImportError: If any required BiocPy package is missing.
    """
    global _BiocFrame, _SummarizedExperiment, _RangedSummarizedExperiment

    if (
        _BiocFrame is None
        or _SummarizedExperiment is None
        or _RangedSummarizedExperiment is None
    ):
        try:  # pragma: no cover - optional, exercised only when installed.
            from biocframe import BiocFrame as BF
            from summarizedexperiment import (
                SummarizedExperiment as SE,
                RangedSummarizedExperiment as RSE,
            )
        except Exception as e:  # pragma: no cover
            raise ImportError(
                "This feature requires BiocPy packages. Install with:\n"
                "  pip install summarizedexperiment biocframe genomicranges\n"
                "or via conda (conda-forge)."
            ) from e

        _BiocFrame = BF
        _SummarizedExperiment = SE
        _RangedSummarizedExperiment = RSE

    # At this point mypy/pyright know these are not None (because of the
    # explicit return type and the guard above).
    return _BiocFrame, _SummarizedExperiment, _RangedSummarizedExperiment


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
    :meth:`recount3.bundle.R3ResourceBundle.to_summarized_experiment`.

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

    BiocFrameCls, SummarizedExperimentCls, RangedSummarizedExperimentCls = _require_biocpy()

    if not isinstance(obj, (SummarizedExperimentCls, RangedSummarizedExperimentCls)):
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
    expanded_bf = BiocFrameCls.from_pandas(expanded_df)
    return obj.set_column_data(expanded_bf)


def _select_shortest_column_match(
    matches: Sequence[tuple[Any, str, str]],
) -> Any:
    """Select the best column match from a list of candidate matches.

    When multiple metadata columns match a canonical name, prefer the shortest
    string representation (it is typically the least-prefixed column).

    Args:
      matches: Candidate matches as tuples (original, string, normalized).

    Returns:
      The original column label (as stored in ``df.columns``).
    """
    best = min(matches, key=lambda item: len(item[1]))
    return best[0]


def compute_read_counts(
    rse: Any,
    round_counts: bool = True,
    avg_mapped_read_length: str = "recount_qc.star.average_mapped_length",
) -> pd.DataFrame:
    """Convert coverage-sum counts into approximate read/fragment counts.

    The "raw_counts" assay used by recount3-style resources represents summed
    per-base coverage over each feature (for example, total coverage across all
    bases in a gene), not the number of reads/fragments overlapping the feature.
    A common approximation for converting coverage-sum values into read/fragment
    counts is to divide by the sample's average mapped read length.

    For each feature i and sample j:

      read_counts[i, j] = raw_counts[i, j] / avg_mapped_read_length[j]

    This operation is performed column-wise: each sample column in the assay is
    divided by a single scalar derived from that sample's metadata.

    Rounding is optional. When enabled, values are rounded to 0 decimals to
    produce integer-like counts, which is convenient for downstream methods that
    assume counts are integer-valued.

    Args:
      rse: A RangedSummarizedExperiment-like object containing a "raw_counts"
        assay and sample metadata in `col_data`.
      round_counts: If True, round the resulting values to 0 decimals.
      avg_mapped_read_length: Name of the metadata column containing average
        mapped read length per sample. The column must be present in `col_data`
        and must contain numeric values.

    Returns:
      A pandas DataFrame of approximate read counts with the same shape as the
      "raw_counts" assay (features x samples). Row and column names are
      preserved when available.

    Raises:
      ValueError: If `rse` is not a RangedSummarizedExperiment, if the
        "raw_counts" assay is missing, if `avg_mapped_read_length` is missing
        from `col_data`, or if the assay and metadata dimensions do not align.
      TypeError: If `round_counts` is not a bool.
    """
    _, _, ranged_summarized_experiment_cls = _require_biocpy()

    if not isinstance(rse, ranged_summarized_experiment_cls):
        raise ValueError(
            "rse must be a RangedSummarizedExperiment to compute read counts."
        )
    if not isinstance(round_counts, bool):
        raise TypeError("round_counts must be a bool.")

    assay_names = getattr(rse, "assay_names", None)
    if not assay_names or "raw_counts" not in assay_names:
        raise ValueError("rse must contain a 'raw_counts' assay.")

    col_data = rse.col_data.to_pandas()
    try:
        avg_len_series = _resolve_metadata_column(col_data, avg_mapped_read_length)
    except ValueError as exc:
        raise ValueError(
            "Required metadata column not found in col_data: "
            f"{avg_mapped_read_length!r}."
        ) from exc

    raw_counts = np.asarray(rse.assay("raw_counts"), dtype=float)
    if raw_counts.ndim != 2:
        raise ValueError(
            "'raw_counts' assay must be a 2D matrix "
            f"(got shape {raw_counts.shape})."
        )

    avg_len_values = avg_len_series.to_numpy()

    for value in avg_len_values:
        if pd.isna(value):
            continue
        if isinstance(value, (str, bytes)):
            raise ValueError(
                f"Metadata column {avg_mapped_read_length!r} must be numeric; "
                f"found {type(value).__name__} value {value!r}."
            )
        if not isinstance(value, numbers.Real):
            raise ValueError(
                f"Metadata column {avg_mapped_read_length!r} must be numeric; "
                f"found {type(value).__name__}."
            )

    avg_len = avg_len_values.astype(float, copy=False)

    if raw_counts.shape[1] != avg_len.shape[0]:
        raise ValueError(
            "Mismatch between number of samples in 'raw_counts' "
            f"({raw_counts.shape[1]}) and length of "
            f"{avg_mapped_read_length!r} ({avg_len.shape[0]})."
        )

    read_counts = raw_counts / avg_len

    if round_counts:
        read_counts = np.rint(read_counts)

    row_names = getattr(rse, "row_names", None)
    col_names = getattr(rse, "col_names", None)

    return pd.DataFrame(
        read_counts,
        index=list(row_names) if row_names is not None else None,
        columns=list(col_names) if col_names is not None else None,
    )


def compute_tpm(rse: RangedSummarizedExperiment) -> pd.DataFrame:
    """Compute Transcripts Per Million (TPM) from raw coverage sums.

    TPM is calculated as:
      1. Approximate Read Counts = Coverage AUC / Avg Read Length
      2. RPK = Read Counts / (Feature Length / 1000)
      3. Scale Factor = Sum(RPK) / 1,000,000
      4. TPM = RPK / Scale Factor

    Args:
      rse: A RangedSummarizedExperiment object containing raw coverage sums.
        Must have feature widths defined in rowRanges.

    Returns:
      A pandas DataFrame of TPM values.

    Raises:
      TypeError: If rse is not a RangedSummarizedExperiment (needs rowRanges).
      ValueError: If feature widths or read lengths are missing.
    """
    _, _, ranged_summarized_experiment_cls = _require_biocpy()

    if not isinstance(rse, ranged_summarized_experiment_cls):
        raise TypeError(
            "compute_tpm requires a RangedSummarizedExperiment with "
            "defined rowRanges to determine feature lengths."
        )

    # 1. Get Read Counts
    reads = compute_read_counts(rse, round_counts=False)

    # 2. Get Feature Lengths (bp)
    # biocutils/genomicranges exposes .width
    try:
        widths = np.array(rse.row_ranges.width)
    except AttributeError as exc:
        raise ValueError(
            "Could not determine feature widths from rowRanges."
        ) from exc

    # 3. Calculate Reads Per Kilobase
    # Numpy broadcasting: (n_features, n_samples) / (n_features, 1)
    kb_widths = widths / 1000.0
    rpk = reads.div(kb_widths, axis=0)

    # 4. Calculate Scaling Factors (per sample)
    # Sum of RPKs per sample / 1,000,000
    scale_factors = rpk.sum(axis=0) / 1_000_000.0

    # 5. Calculate TPM
    tpm = rpk.div(scale_factors, axis=1)

    return tpm


def _coerce_col_data_to_pandas(x: Any) -> pd.DataFrame:
    """Coerce sample metadata into a pandas DataFrame.

    Args:
        x: Either a BiocPy (Ranged)SummarizedExperiment-like object with a
          `.col_data.to_pandas()` method, or a pandas DataFrame already.

    Returns:
        A pandas DataFrame of sample metadata.

    Raises:
        TypeError: If `x` cannot be coerced to a pandas DataFrame.
    """
    if isinstance(x, pd.DataFrame):
        return x

    # If `x` is a (Ranged)SummarizedExperiment, use col_data like the R code.
    if hasattr(x, "col_data") and hasattr(x.col_data, "to_pandas"):
        return x.col_data.to_pandas()

    raise TypeError(
        "Expected a pandas.DataFrame or a SummarizedExperiment-like object with "
        "`col_data.to_pandas()`."
    )


def _resolve_metadata_column(
    metadata: pd.DataFrame,
    column_name: str,
) -> pd.Series:
    """Resolve a metadata column name robustly.

    This mirrors the strictness of the recount3 R implementation (which expects
    exact column names), but also supports the Python-side convention where the
    namespace separator may be `__` instead of `.` for the first separator
    (e.g., `recount_qc.star.average_mapped_length` vs
    `recount_qc__star.average_mapped_length`).

    Args:
        metadata: Sample metadata.
        column_name: Column name to resolve.

    Returns:
        The resolved column as a pandas Series.

    Raises:
        ValueError: If the column cannot be found.
    """
    # Case-insensitive exact match.
    lower_to_actual = {str(c).lower(): c for c in metadata.columns}
    key = column_name.lower()
    if key in lower_to_actual:
        return metadata[lower_to_actual[key]]

    # Try swapping only the *first* namespace separator '.' -> '__'
    # to match conventions used elsewhere in this Python package.
    if "." in column_name:
        namespace, rest = column_name.split(".", 1)
        alt = f"{namespace}__{rest}".lower()
        if alt in lower_to_actual:
            return metadata[lower_to_actual[alt]]

    raise ValueError(
        f"Required metadata column {column_name!r} not found. "
        "If your metadata uses '__' as a namespace separator, pass the "
        "actual column name explicitly."
    )


def _coerce_numeric_column(series: pd.Series, column_name: str) -> pd.Series:
    """Coerce a Series to numeric, erroring on non-numeric non-missing values.

    Args:
        series: Input Series.
        column_name: Name used for error messages.

    Returns:
        Float Series.

    Raises:
        ValueError: If non-missing values cannot be coerced to numeric.
    """
    numeric = pd.to_numeric(series, errors="coerce")
    invalid = series.notna() & numeric.isna()
    if invalid.any():
        examples = series[invalid].head(3).tolist()
        raise ValueError(
            f"Metadata column {column_name!r} contains non-numeric values "
            f"(examples: {examples!r})."
        )
    return numeric.astype(float)


def is_paired_end(
    x: Any,
    avg_mapped_read_length: str = "recount_qc.star.average_mapped_length",
    avg_read_length: str = "recount_seq_qc.avg_len",
) -> pd.Series:
    """Infer paired-end status, matching recount3::is_paired_end().

    In recount3 (R), paired-end status is inferred via:
      ratio <- round(avg_mapped_read_length / avg_read_length, 0)
      ratio must be 1 (single-end) or 2 (paired-end), otherwise NA with warning.
      result <- ratio == 2, with names(result) = external_id.

    Args:
        x: Sample metadata (DataFrame) or a (Ranged)SummarizedExperiment-like
          object with `col_data.to_pandas()`.
        avg_mapped_read_length: Metadata column containing average mapped length.
        avg_read_length: Metadata column containing average read length.

    Returns:
        A pandas Series of dtype "boolean" (True/False/pd.NA), indexed by
        external_id.

    Raises:
        ValueError: If required metadata columns are missing or non-numeric.
    """
    metadata = _coerce_col_data_to_pandas(x)

    external_id = _resolve_metadata_column(metadata, "external_id").astype(str)
    avg_mapped = _coerce_numeric_column(
        _resolve_metadata_column(metadata, avg_mapped_read_length),
        avg_mapped_read_length,
    )
    avg_len = _coerce_numeric_column(
        _resolve_metadata_column(metadata, avg_read_length),
        avg_read_length,
    )

    # R uses round(..., 0). numpy.rint matches "round to nearest, ties to even".
    ratio = np.rint(avg_mapped.to_numpy() / avg_len.to_numpy())
    ratio = pd.Series(ratio, index=external_id, dtype=float)

    invalid = ratio.notna() & ~ratio.isin([1.0, 2.0])
    if invalid.any():
        logging.warning(
            "Unable to determine if samples are paired-end for %d sample(s). "
            "Setting paired_end=NA for those.",
            int(invalid.sum()),
        )
        ratio.loc[invalid] = np.nan

    paired_end = pd.Series(pd.NA, index=external_id, dtype="boolean")
    paired_end.loc[ratio == 2.0] = True
    paired_end.loc[ratio == 1.0] = False
    return paired_end


def compute_scale_factors(
    x: Any,
    by: str = "auc",
    target_size: float = 4e7,
    L: float = 100,
    auc: str = "recount_qc.bc_auc.all_reads_all_bases",
    avg_mapped_read_length: str = "recount_qc.star.average_mapped_length",
    mapped_reads: str = "recount_qc.star.all_mapped_reads",
    paired_end: Sequence[bool] | pd.Series | None = None,
) -> pd.Series:
    """Compute per-sample scaling factors for coverage-sum counts.

    This function produces one scalar scale factor per sample. The intended use
    is to multiply each sample column of a coverage-sum count matrix by the
    corresponding factor to make samples comparable.

    Let C[i, j] be the coverage-sum count for feature i in sample j. If s[j] is
    the scale factor for sample j, scaled counts are computed as:

      scaled[i, j] = C[i, j] * s[j]

    Scale factors are derived from sample metadata. Samples are identified by
    the `external_id` column, and the returned Series is indexed by external_id.

    Two scaling methods are supported:

    1) by="auc"
       Uses a per-sample total coverage metric (AUC) to scale each sample to a
       common target_size:

         s[j] = target_size / auc[j]

       This method preserves relative feature coverage within each sample while
       adjusting overall sample magnitude to be comparable across samples.

    2) by="mapped_reads"
       Uses mapped read counts and read length to normalize samples to a common
       target_size and a common target read length L:

         s[j] = target_size * L * paired_multiplier[j] /
                (mapped_reads[j] * (avg_mapped_read_length[j] ** 2))

       paired_multiplier is:
         - 2 for paired-end samples
         - 1 for single-end samples
         - missing for samples whose paired-end status cannot be inferred

       If `paired_end` is not provided, paired-end status is inferred from
       metadata by comparing average mapped length to average read length:
         ratio = round(avg_mapped_read_length / avg_read_length)
       ratio==2 indicates paired-end, ratio==1 indicates single-end. Other
       ratios are treated as unknown and produce missing paired multipliers.

    Missing values in required metadata propagate to missing scale factors.
    Non-numeric metadata values raise an error.

    Args:
      x: Sample metadata as a pandas DataFrame, or an object with
        `col_data.to_pandas()` that yields sample metadata.
      by: Scaling method: "auc" or "mapped_reads".
      target_size: Target library size used to compute scale factors. Interpreted
        as the number of single-end reads to scale each sample to.
      L: Target read length used only when by="mapped_reads".
      auc: Metadata column name for the per-sample AUC metric.
      avg_mapped_read_length: Metadata column name for average mapped read length
        per sample.
      mapped_reads: Metadata column name for mapped read counts per sample.
      paired_end: Optional paired-end indicator per sample. If provided, it must
        align with the samples in `external_id`. If omitted, paired-end status
        is inferred from metadata.

    Returns:
      A pandas Series of scale factors indexed by `external_id`. The Series name
      is "scale_factor".

    Raises:
      ValueError: If `by` is invalid, required metadata columns are missing, or
        non-numeric metadata values are present.
      TypeError: If `target_size` or `L` are not numeric scalars.
    """
    if by not in ("auc", "mapped_reads"):
        raise ValueError(f"{by=!r} must be either 'auc' or 'mapped_reads'.")

    if not isinstance(target_size, numbers.Real) or isinstance(target_size, bool):
        raise TypeError("target_size must be a numeric scalar.")
    if not isinstance(L, numbers.Real) or isinstance(L, bool):
        raise TypeError("L must be a numeric scalar.")

    metadata = _coerce_col_data_to_pandas(x)

    # Match recount3's stopifnot(): all of these must exist even if by="auc".
    external_id = _resolve_metadata_column(metadata, "external_id").astype(str)

    auc_values = _coerce_numeric_column(
        _resolve_metadata_column(metadata, auc),
        auc,
    )
    avg_mapped_values = _coerce_numeric_column(
        _resolve_metadata_column(metadata, avg_mapped_read_length),
        avg_mapped_read_length,
    )
    mapped_reads_values = _coerce_numeric_column(
        _resolve_metadata_column(metadata, mapped_reads),
        mapped_reads,
    )

    auc_values.index = external_id
    avg_mapped_values.index = external_id
    mapped_reads_values.index = external_id

    if paired_end is None:
        paired_end_series = is_paired_end(
            metadata,
            avg_mapped_read_length=avg_mapped_read_length,
        )
    else:
        paired_end_series = pd.Series(
            paired_end,
            index=external_id,
            dtype="boolean",
        )

    if by == "auc":
        scale_factor = float(target_size) / auc_values
    else:
        pe_multiplier = pd.Series(np.nan, index=external_id, dtype=float)
        pe_multiplier.loc[paired_end_series == True] = 2.0
        pe_multiplier.loc[paired_end_series == False] = 1.0

        denom = mapped_reads_values * (avg_mapped_values**2)
        scale_factor = float(target_size) * float(L) * pe_multiplier / denom

    scale_factor.name = "scale_factor"
    return scale_factor


def transform_counts(
    rse: Any,
    by: str = "auc",
    target_size: float = 4e7,
    L: float = 100,
    round_counts: bool = True,
    **kwargs: Any,
) -> pd.DataFrame:
    """Scale coverage-sum counts to a common library size.

    recount3 "raw_counts" represent summed per-base coverage over each feature,
    not read/fragment counts. This function converts those coverage-sum values
    into scaled counts that are comparable across samples by multiplying each
    sample column by a sample-specific scale factor.

    Scaling is applied independently per sample (per column). For feature i and
    sample j, the returned matrix contains:

      scaled[i, j] = raw_counts[i, j] * scale_factor[j]

    The scale factors are computed from sample metadata (col_data) using one of
    two methods:

    - by="auc":
        scale_factor[j] = target_size / auc[j]

      where auc is a per-sample total coverage metric. This method scales each
      sample to have total coverage approximately equal to target_size.

    - by="mapped_reads":
        scale_factor[j] = (
            target_size * L * paired_multiplier[j]
            / (mapped_reads[j] * (avg_mapped_read_length[j] ** 2))
        )

      where paired_multiplier is 2 for paired-end samples, 1 for single-end
      samples, and missing when paired-end status cannot be inferred. This
      method incorporates mapped reads and read length so that samples with
      different read lengths are normalized onto the same target read length L.

    The returned values remain in the same feature-by-sample shape as the
    input. If round_counts is True, values are rounded to integer-like counts.

    Args:
      rse: A RangedSummarizedExperiment-like object containing a "raw_counts"
        assay and sample metadata in `col_data`.
      by: Scaling method: "auc" or "mapped_reads".
      target_size: Target library size used to compute scale factors. Interpreted
        as the number of single-end reads to scale each sample to.
      L: Target read length used only when by="mapped_reads".
      round_counts: If True, round scaled values to 0 decimals.
      **kwargs: Additional parameters forwarded to compute_scale_factors(). Use
        this to override metadata column names (for example, `auc=...`,
        `mapped_reads=...`, `avg_mapped_read_length=...`) or to provide
        `paired_end=...` when paired-end status should not be inferred.

    Returns:
      A pandas DataFrame of scaled counts with the same dimensions as
      `assay("raw_counts")`. Row and column names are preserved when available.

    Raises:
      ValueError: If `rse` is not a RangedSummarizedExperiment, if the required
        assay or metadata columns are missing, if `by` is invalid, or if the
        assay and metadata dimensions do not align.
      TypeError: If `round_counts` is not a bool, or if numeric parameters are
        not valid scalars.
    """
    _, _, RangedSummarizedExperiment = _require_biocpy()

    if not isinstance(rse, RangedSummarizedExperiment):
        raise ValueError(
            "rse must be a RangedSummarizedExperiment (BiocPy "
            "summarizedexperiment)."
        )

    assay_names = getattr(rse, "assay_names", None)
    if assay_names is None or "raw_counts" not in assay_names:
        raise ValueError("rse must contain a 'raw_counts' assay.")

    if not isinstance(round_counts, bool):
        raise ValueError("round_counts must be a boolean scalar.")

    counts = rse.assay("raw_counts")
    counts_array = np.asarray(counts, dtype=float)

    scale_factor = compute_scale_factors(
        rse,
        by=by,
        target_size=target_size,
        L=L,
        **kwargs,
    )

    if counts_array.ndim != 2:
        raise ValueError(
            f"'raw_counts' assay must be 2D (got shape {counts_array.shape})."
        )
    if counts_array.shape[1] != len(scale_factor):
        raise ValueError(
            "Mismatch between number of samples in 'raw_counts' "
            f"({counts_array.shape[1]}) and number of scale factors "
            f"({len(scale_factor)})."
        )

    scaled = counts_array * scale_factor.to_numpy(dtype=float)

    if round_counts:
        scaled = np.rint(scaled)

    row_names = getattr(rse, "row_names", None)
    col_names = getattr(rse, "col_names", None)

    return pd.DataFrame(
        scaled,
        index=list(row_names) if row_names is not None else None,
        columns=list(col_names) if col_names is not None else None,
    )
