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
import numpy as np
import pandas as pd

from recount3.bundle import R3ResourceBundle
from recount3.search import annotation_ext

if TYPE_CHECKING:  # only for static type checkers
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


def _find_metadata_column(
    df: pd.DataFrame,
    candidates: Sequence[str],
) -> pd.Series | None:
    """Find a column in a DataFrame matching one of the candidate names.

    Matches exact names first, then suffix matches (e.g., 'avg_len' matches
    'recount_seq_qc.avg_len'). This handles the variable prefixing behavior
    of recount3 metadata tables.

    Args:
      df: The DataFrame to search.
      candidates: A sequence of column names to look for.

    Returns:
      The matching Series if found, or None.
    """
    columns = [str(c) for c in df.columns]
    
    # 1. Exact match
    for cand in candidates:
        if cand in columns:
            return df[cand]

    # 2. Suffix match (e.g. 'recount_seq_qc.avg_len' matches 'avg_len')
    for cand in candidates:
        suffix = f".{cand}"
        matches = [col for col in columns if col.endswith(suffix)]
        if matches:
            # If multiple match, pick the shortest (most likely the canonical one)
            # or the first one.
            best = min(matches, key=len)
            return df[best]
            
    return None


def compute_read_counts(
    rse: SummarizedExperiment,
    *,
    round_counts: bool = True,
) -> pd.DataFrame:
    """Compute estimated read counts from raw coverage sums.

    recount3 stores data as "coverage sums" (Area Under Coverage). This function
    converts them to approximate read counts by dividing by the average read
    length per sample.

    Args:
      rse: A SummarizedExperiment or RangedSummarizedExperiment object
        containing the raw coverage counts in the first assay.
      round_counts: If True, round the resulting counts to the nearest integer.

    Returns:
      A pandas DataFrame of estimated read counts.

    Raises:
      ValueError: If the read length metadata cannot be found in col_data.
      ImportError: If BiocPy packages are not installed.
    """
    _require_biocpy()

    # 1. Get raw coverage (AUC)
    # Assumes the first assay is the counts matrix.
    if not rse.assays:
        raise ValueError("RSE object contains no assays.")
    
    # Access first assay regardless of name ("counts", "raw_counts", etc.)
    assay_name = rse.assay_names[0]
    coverage_auc = pd.DataFrame(
        rse.assay(assay_name),
        index=rse.row_names,
        columns=rse.col_names
    )

    # 2. Get average read length
    # Check col_data for standard keys used in recount3 qc tables.
    col_data = rse.col_data.to_pandas()
    candidates = ("avg_len", "avg_read_length", "read_length")
    avg_len = _find_metadata_column(col_data, candidates)

    if avg_len is None:
        raise ValueError(
            "Could not find average read length in col_data. "
            f"Searched for: {candidates}. "
            "Ensure 'recount_seq_qc' metadata is included in the bundle."
        )

    # 3. Compute Read Counts = AUC / AvgReadLength
    # Align avg_len to the columns of the coverage matrix
    avg_len = avg_len.reindex(coverage_auc.columns)
    
    # Handle NaN lengths (div by zero protection)
    if avg_len.isna().any():
        logging.warning(
            "Missing read lengths for %d samples; counts will be NaN.",
            avg_len.isna().sum()
        )

    read_counts = coverage_auc.div(avg_len, axis=1)

    if round_counts:
        read_counts = read_counts.round()

    return read_counts


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
    _require_biocpy()

    if not isinstance(rse, RangedSummarizedExperiment):
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


def transform_counts(
    rse: SummarizedExperiment,
    *,
    by: str = "raw",
    round_counts: bool = True,
) -> SummarizedExperiment:
    """Transform the counts in a SummarizedExperiment.

    Mirrors the behavior of recount3::transform_counts() in R. It returns a
    new SummarizedExperiment with the assay transformed to the requested unit.

    Args:
      rse: Input SummarizedExperiment or RangedSummarizedExperiment with
        raw coverage sums (AUC) in the first assay.
      by: Transformation type.
        * "raw": Return object unchanged (coverage sums).
        * "read_counts": Estimate read counts (AUC / avg_read_len).
        * "tpm": Transcripts Per Million (requires feature lengths).
      round_counts: If True, round read_counts to integers. Ignored for TPM.

    Returns:
      A new SummarizedExperiment containing the transformed assay.
      The assay name is updated to the value of `by`.

    Raises:
      ValueError: If `by` is unknown or metadata is missing.
    """
    _require_biocpy()

    valid_modes = {"raw", "read_counts", "tpm"}
    if by not in valid_modes:
        raise ValueError(f"Invalid transformation 'by={by!r}'. Expected: {valid_modes}")

    if by == "raw":
        return rse

    if by == "read_counts":
        new_matrix = compute_read_counts(rse, round_counts=round_counts)
    elif by == "tpm":
        # TPM usually implies RSE input
        if not isinstance(rse, RangedSummarizedExperiment):
             raise TypeError("Transform 'tpm' requires a RangedSummarizedExperiment.")
        new_matrix = compute_tpm(rse)
    else:
        # Should be unreachable due to check above
        raise ValueError(f"Unknown mode: {by}")

    # Construct new object with transformed assay
    # Prefer to keep the existing metadata, col_data, etc.
    # SE immutable-ish in design; construct a new one.
    
    # BiocPy constructors are flexible.
    # Assume 'new_matrix' is a DataFrame; convert to numpy for the SE.
    new_assay = new_matrix.to_numpy(copy=False)
    
    # Use the specific class constructor to preserve RSE vs SE nature
    cls = type(rse)
    
    # Helper to extract row_ranges only if it exists (RSE)
    kwargs = {
        "assays": {by: new_assay},
        "row_data": rse.row_data,
        "col_data": rse.col_data,
        "metadata": rse.metadata
    }
    
    if isinstance(rse, RangedSummarizedExperiment):
        kwargs["row_ranges"] = rse.row_ranges

    return cls(**kwargs)
