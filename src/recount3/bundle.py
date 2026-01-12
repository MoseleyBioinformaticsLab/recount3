"""Resource bundles, project discovery, and concatenation helpers.

This module defines :class:`R3ResourceBundle`, a general-purpose
container for groups of :class:`recount3.resource.R3Resource` objects.

Bundles support lazy loading, filtering by description fields,
project-aware discovery, and high-level helpers for combining recount3
resources into BiocPy objects such as
:class:`summarizedexperiment.SummarizedExperiment` and
:class:`summarizedexperiment.RangedSummarizedExperiment`.

Typical usage example:

  from recount3.bundle import R3ResourceBundle

  bundle = R3ResourceBundle.discover(
      organism="human",
      data_source="sra",
      project="SRP009615",
  )

  # Stack raw count matrices across samples.
  counts = bundle.only_counts().stack_count_matrices()

  # Build a RangedSummarizedExperiment for gene-level counts.
  rse = bundle.to_ranged_summarized_experiment(genomic_unit="gene")
"""

from __future__ import annotations

import dataclasses
import functools
import logging
from collections.abc import Callable, Iterable, Iterator, Sequence
from typing import Any, Optional, TYPE_CHECKING

import pandas as pd

from recount3 import _bigwig
from recount3 import errors
from recount3 import resource
from recount3 import search
from recount3 import types as r3_types

# Optional BiocPy imports. These are used only when constructing
# SummarizedExperiment-style objects and are not required for core
# bundle functionality.
if TYPE_CHECKING:
    from biocframe import BiocFrame  # type: ignore[import-not-found]
    from genomicranges import GenomicRanges  # type: ignore[import-not-found]
    from summarizedexperiment import (  # type: ignore[import-not-found]
        RangedSummarizedExperiment,
        SummarizedExperiment,
    )

# Runtime caches for optional BiocPy classes.
_BiocFrame: type[Any] | None = None
_GenomicRanges: type[Any] | None = None
_SummarizedExperiment: type[Any] | None = None
_RangedSummarizedExperiment: type[Any] | None = None


def _require_biocpy() -> tuple[type[Any], type[Any], type[Any], type[Any]]:
    """Ensure that the BiocPy packages used here are importable.

    Returns:
      (BiocFrame, GenomicRanges, SummarizedExperiment, RangedSummarizedExperiment)
      classes.

    Raises:
      ImportError: If any required BiocPy package is missing.
    """
    global _BiocFrame, _GenomicRanges, _SummarizedExperiment, _RangedSummarizedExperiment

    if (
        _BiocFrame is None
        or _GenomicRanges is None
        or _SummarizedExperiment is None
        or _RangedSummarizedExperiment is None
    ):
        try:  # pragma: no cover - exercised only when BiocPy is installed.
            from biocframe import BiocFrame as BF
            from genomicranges import GenomicRanges as GR
            from summarizedexperiment import (
                SummarizedExperiment as SE,
                RangedSummarizedExperiment as RSE,
            )
        except Exception as e:  # pragma: no cover
            raise ImportError(
                "BiocPy integration requires the 'summarizedexperiment', "
                "'biocframe', and 'genomicranges' packages. Install them with:\n"
                "  pip install summarizedexperiment biocframe genomicranges\n"
                "or via conda (conda-forge)."
            ) from e

        _BiocFrame = BF
        _GenomicRanges = GR
        _SummarizedExperiment = SE
        _RangedSummarizedExperiment = RSE

    # Help type checkers see these as non-optional
    assert _BiocFrame is not None
    assert _GenomicRanges is not None
    assert _SummarizedExperiment is not None
    assert _RangedSummarizedExperiment is not None

    return (
        _BiocFrame,
        _GenomicRanges,
        _SummarizedExperiment,
        _RangedSummarizedExperiment,
    )


def _ensure_unique_columns(
    df: pd.DataFrame,
    *,
    empty_prefix: str = "col",
) -> pd.DataFrame:
    """Return a copy with unique, non-empty string column names.

    Many recount3 metadata tables share column names (for example,
    ``external_id``) and concatenation can introduce duplicates. BiocFrame
    drops duplicated names, which then breaks downstream validation. This
    helper ensures that:

    * All column names are strings (``None`` becomes an empty string).
    * Empty names are replaced by ``{empty_prefix}``.
    * Duplicates are suffixed as ``name__2``, ``name__3``, and so on.

    Args:
      df: Input :class:`pandas.DataFrame` whose columns may contain
        duplicates or empty names.
      empty_prefix: Base name used when an empty or ``None`` column name
        is encountered.

    Returns:
      A copy of ``df`` with deduplicated, non-empty column names.
    """
    out = df.copy()
    raw_cols = [("" if c is None else str(c)) for c in out.columns]
    counts: dict[str, int] = {}
    new_cols: list[str] = []

    for name in raw_cols:
        base = name or empty_prefix
        n = counts.get(base, 0) + 1
        counts[base] = n
        new_cols.append(base if n == 1 else f"{base}__{n}")

    out.columns = new_cols
    return out


_METADATA_MERGE_KEYS = ("rail_id", "external_id", "study")
_METADATA_KEY_SYNONYMS = {
    # Backwards-compatibility with older recount3 metadata.
    "study_acc": "study",
    "run_acc": "external_id",
    "run_accession": "external_id",
    "run": "external_id",
}
_METADATA_NAMESPACE_SEPARATOR = "__"


def _standardize_metadata_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize metadata column names and key fields.

    This function normalizes column names to lower-case strings, applies a small
    set of known backwards-compatible key renames, and ensures the standard key
    columns exist.

    Args:
      df: Raw metadata table.

    Returns:
      A copy of ``df`` with standardized columns.
    """
    out = df.copy()
    out.columns = [str(col).strip().lower() for col in out.columns]

    rename: dict[str, str] = {}
    for old, new in _METADATA_KEY_SYNONYMS.items():
        if old in out.columns and new not in out.columns:
            rename[old] = new
    if rename:
        out = out.rename(columns=rename)

    for key in _METADATA_MERGE_KEYS:
        if key not in out.columns:
            out[key] = pd.NA

    for key in _METADATA_MERGE_KEYS:
        out[key] = out[key].astype("string")

    return out


def _metadata_origin(res: resource.R3Resource) -> str:
    """Return a stable origin prefix for a metadata resource.

    Args:
      res: A metadata resource.

    Returns:
      A short string used to namespace the resource's non-key columns.
    """
    desc = res.description
    origin = (
        getattr(desc, "table_name", None)
        or getattr(desc, "resource_type", None)
    )
    if not origin:
        origin = "metadata"
    return str(origin).strip().lower()


def _namespace_metadata_columns(
    df: pd.DataFrame,
    *,
    origin: str,
    keys: Sequence[str] = _METADATA_MERGE_KEYS,
    sep: str = _METADATA_NAMESPACE_SEPARATOR,
) -> tuple[pd.DataFrame, dict[str, tuple[str, str]]]:
    """Namespace non-key columns with an origin prefix.

    Args:
      df: Standardized metadata table.
      origin: Prefix string (e.g. ``"recount_qc"``).
      keys: Column names to leave unmodified.
      sep: Separator used between origin and the original column name.

    Returns:
      A tuple (namespaced_df, provenance) where provenance maps new column names
      to (origin, original_name).
    """
    provenance: dict[str, tuple[str, str]] = {}
    rename: dict[str, str] = {}

    for col in df.columns:
        if col in keys:
            continue
        new_name = f"{origin}{sep}{col}"
        rename[col] = new_name
        provenance[new_name] = (origin, col)

    out = df.rename(columns=rename)
    out = _ensure_unique_columns(out, empty_prefix="col")
    return out, provenance


def _outer_merge_metadata_frames(
    frames: Sequence[pd.DataFrame],
    *,
    keys: Sequence[str] = _METADATA_MERGE_KEYS,
) -> pd.DataFrame:
    """Outer-merge metadata frames on a fixed set of key columns."""
    if not frames:
        return pd.DataFrame(columns=list(keys))

    def merge_two(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
        return pd.merge(left, right, on=list(keys), how="outer")

    return functools.reduce(merge_two, frames)


def _choose_alignment_key(
    *,
    sample_ids: Sequence[str],
    merged: pd.DataFrame,
) -> str:
    """Choose the metadata column that best matches the assay sample IDs."""
    candidates = ("external_id", "rail_id")
    sample_set = set(sample_ids)

    best_key = "external_id"
    best_matches = -1
    for key in candidates:
        if key not in merged.columns:
            continue
        matches = merged[key].astype("string").isin(sample_set).sum()
        if matches > best_matches:
            best_matches = int(matches)
            best_key = key
    return best_key


def _collapse_rows_by_key(df: pd.DataFrame, *, key: str) -> pd.DataFrame:
    """Collapse duplicated keys by taking the first non-null value per column."""
    if key not in df.columns:
        return df

    def first_non_null(values: pd.Series) -> Any:
        non_null = values.dropna()
        if non_null.empty:
            return pd.NA
        return non_null.iloc[0]

    # groupby(..., dropna=False) keeps rows with missing keys as a group.
    collapsed = (
        df.groupby(key, dropna=False, sort=False)
        .aggregate(first_non_null)
        .reset_index()
    )
    return collapsed


def _maybe_relabel_counts_columns_to_external_id(
    counts_df: pd.DataFrame,
    col_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Relabel assay columns to external_id when a complete mapping exists."""
    if "external_id" not in col_df.columns:
        return counts_df, col_df

    external_ids = list(col_df["external_id"].astype("string"))
    if len(external_ids) != len(counts_df.columns):
        return counts_df, col_df

    missing_external_id = any(
        val is pd.NA or val is None or str(val).strip() == ""
        for val in external_ids
    )
    if missing_external_id:
        return counts_df, col_df

    external_ids_str = [str(x) for x in external_ids]
    if len(set(external_ids_str)) != len(external_ids_str):
        return counts_df, col_df

    current_ids = [str(c) for c in counts_df.columns]
    if current_ids == external_ids_str:
        return counts_df, col_df

    renamed_counts = counts_df.copy()
    renamed_counts.columns = external_ids_str

    renamed_col = col_df.copy()
    renamed_col.index = external_ids_str
    renamed_col["external_id"] = external_ids_str

    return renamed_counts, renamed_col


def _read_rr_table(res: resource.R3Resource) -> pd.DataFrame:
    """Read an RR-style junction coordinate table into a DataFrame.

    The RR files in recount3 are small TSV(.gz) tables with one row per
    junction containing genomic coordinates and an identifier.

    Args:
      res: A junction resource whose description indicates an RR file.

    Returns:
      A :class:`pandas.DataFrame` containing the parsed RR table.

    Raises:
      ValueError: If the resource cannot be read via its cache path.
    """
    try:
        df = res.load()
        if isinstance(df, pd.DataFrame):
            return df
    except Exception:  # pragma: no cover - fall back to manual parsing.
        pass

    try:
        path = res._cached_path()  # pylint: disable=protected-access
    except Exception as exc:  # pragma: no cover - defensive.
        raise ValueError(f"Cannot resolve cached RR path for: {res.url}") from exc

    import gzip
    import io

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:  # type: ignore[arg-type]
        text = io.TextIOWrapper(fh, encoding="utf-8")
        df = pd.read_csv(text, sep="\t", comment="#")

    return df


def _read_gtf_dataframe(res: resource.R3Resource) -> pd.DataFrame:
    """Read a GTF(.gz) annotation resource into a DataFrame.

    The result contains the standard 9 GTF columns, including the
    ``attributes`` field.

    Args:
      res: An annotation resource describing a GTF or GTF.GZ file.

    Returns:
      A :class:`pandas.DataFrame` with columns::

          seqname, source, feature, start, end, score, strand,
          frame, attributes

    Raises:
      ValueError: If the cached path for the resource cannot be resolved.
    """
    try:
        path = res._cached_path()  # pylint: disable=protected-access
    except Exception as exc:
        raise ValueError(f"Cannot resolve cached GTF path for: {res.url}") from exc

    import gzip
    import io

    cols = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attributes",
    ]
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:  # type: ignore[arg-type]
        text = io.TextIOWrapper(fh, encoding="utf-8")
        df = pd.read_csv(
            text,
            sep="\t",
            comment="#",
            header=None,
            names=cols,
            dtype=str,
        )

    for col in ("start", "end"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def _ranges_from_gtf(
    gtf: pd.DataFrame,
    *,
    feature_kind: str,
) -> pd.DataFrame:
    """Convert a GTF table to feature ranges for genes or exons.

    Args:
      gtf: GTF DataFrame with an ``attributes`` column.
      feature_kind: Either ``"gene"`` or ``"exon"``.

    Returns:
      A DataFrame with columns ``["seqnames", "starts", "ends",
      "strand", "feature_id"]``. The ``feature_id`` column is derived
      from the GTF attributes (for example, ``gene_id`` or ``exon_id``).
    """
    df = gtf.loc[gtf["feature"] == feature_kind].copy()
    if df.empty:
        return pd.DataFrame(
            columns=["seqnames", "starts", "ends", "strand", "feature_id"]
        )

    import re

    def _parse_attr(s: str) -> dict[str, str]:
        items: dict[str, str] = {}
        for part in re.split(r";\s*", s.strip().rstrip(";")):
            if not part:
                continue
            match = re.match(r'([^ ]+)\s+"([^"]+)"', part)
            if match:
                items[match.group(1)] = match.group(2)
        return items

    attrs = df["attributes"].astype(str).map(_parse_attr)

    if feature_kind == "gene":
        feature_ids = [a.get("gene_id") for a in attrs]
    else:
        feature_ids: list[str] = []
        for attr, (_, row) in zip(attrs, df.iterrows()):
            eid = attr.get("exon_id")
            if eid:
                feature_ids.append(eid)
                continue
            # Fallback: coord-based ID; may not match counts index but is
            # better than NaN.
            feature_ids.append(
                f"{row['seqname']}:{row['start']}-{row['end']}:{row['strand']}"
            )

    out = pd.DataFrame(
        {
            "seqnames": df["seqname"].astype(str).values,
            "starts": df["start"].astype(int).values,
            "ends": df["end"].astype(int).values,
            "strand": df["strand"].astype(str).values,
            "feature_id": feature_ids,
        }
    )
    return out


def _to_genomic_ranges(ranges_df: pd.DataFrame) -> Any:
    """Construct a GenomicRanges object from a pandas DataFrame.

    This helper tries multiple constructor variants to support different
    versions of the :mod:`genomicranges` package.

    Args:
      ranges_df: DataFrame with at least ``seqnames``, ``starts``,
        ``ends``, and ``strand`` columns.

    Returns:
      A :class:`genomicranges.GenomicRanges` instance.

    Raises:
      ImportError: If :mod:`genomicranges` cannot be imported.
      TypeError: If no supported constructor is available.
    """
    _, GenomicRangesCls, _, _ = _require_biocpy()

    # Preferred recent API: classmethod from_pandas.
    if hasattr(GenomicRangesCls, "from_pandas"):
        return GenomicRangesCls.from_pandas(ranges_df)

    # Older classmethod name.
    if hasattr(GenomicRangesCls, "fromPandas"):  # pragma: no cover - legacy.
        return GenomicRangesCls.fromPandas(ranges_df)

    # Module-level constructor used in some releases.
    try:  # pragma: no cover - optional compatibility path.
        import genomicranges as genomicranges_module  # type: ignore[import-not-found]
    except Exception as exc:
        raise ImportError(
            "GenomicRanges is installed but no 'from_pandas' style "
            "constructor could be found."
        ) from exc

    if hasattr(genomicranges_module, "from_pandas"):
        return genomicranges_module.from_pandas(ranges_df)  # type: ignore

    raise TypeError(
        "Could not construct a GenomicRanges object from a pandas DataFrame. "
        "No 'from_pandas' or 'fromPandas' constructor was found."
    )


def _construct_se_compat(
    *,
    counts_df: pd.DataFrame,
    row_df: pd.DataFrame,
    col_df: pd.DataFrame,
    assay_name: str,
) -> SummarizedExperiment:
    """Construct a SummarizedExperiment using compatibility-aware variants.

    This function enforces strict shape and name alignment and then tries
    several constructor signatures supported across different versions
    of the BiocPy :mod:`summarizedexperiment` package.

    The following variants are attempted in order:

      1. ``SummarizedExperiment(assays={name: ndarray}, row_data=df,
         column_data=df)``
      2. ``SummarizedExperiment(assays=[ndarray], assay_names=[name],
         row_data=df, column_data=df)``
      3. Same as (1) but using ``col_data=`` instead of ``column_data=``.

    Args:
      counts_df: Feature-by-sample counts matrix.
      row_df: Row metadata. The index is replaced by
        ``counts_df.index``.
      col_df: Column metadata. The index is replaced by
        ``counts_df.columns``.
      assay_name: Name of the assay, for example, ``"counts"``.

    Returns:
      A :class:`summarizedexperiment.SummarizedExperiment` instance.

    Raises:
      ValueError: If shapes between assay and metadata are inconsistent
        or the assay matrix is empty.
      TypeError: If no constructor variant is accepted by the installed
        :mod:`summarizedexperiment` package.
    """
    _, _, SummarizedExperimentCls, _ = _require_biocpy()

    n_features, n_samples = counts_df.shape
    if n_features == 0 or n_samples == 0:
        raise ValueError("Empty assay matrix; cannot build SummarizedExperiment.")

    row_data = row_df.copy()
    col_data = col_df.copy()
    row_data.index = counts_df.index
    col_data.index = counts_df.columns
    row_data = _ensure_unique_columns(row_data, empty_prefix="row")
    col_data = _ensure_unique_columns(col_data, empty_prefix="col")

    if len(row_data) != n_features:
        raise ValueError(
            f"row_data length {len(row_data)} != assay features {n_features}"
        )
    if len(col_data) != n_samples:
        raise ValueError(
            f"column_data length {len(col_data)} != assay samples {n_samples}"
        )

    counts_numeric = counts_df.apply(pd.to_numeric, errors="coerce").fillna(0)
    counts_np = counts_numeric.to_numpy(copy=False)

    errors_seen: list[str] = []

    # Variant 1: dict[str, ndarray] + column_data
    try:
        return SummarizedExperimentCls(
            assays={assay_name: counts_np},
            row_data=row_data,
            column_data=col_data,
        )
    except Exception as exc:  # pylint: disable=broad-except
        errors_seen.append(f"dict[str, ndarray] + column_data → {exc!r}")

    # Variant 2: list[ndarray] + assay_names + column_data
    try:
        return SummarizedExperimentCls(
            assays=[counts_np],
            assay_names=[assay_name],
            row_data=row_data,
            column_data=col_data,
        )
    except Exception as exc:  # pylint: disable=broad-except
        errors_seen.append(
            "list[ndarray] + assay_names + column_data "
            f"→ {exc!r}"
        )

    # Variant 3: swap keyword to col_data
    try:
        return SummarizedExperimentCls(
            assays={assay_name: counts_np},
            row_data=row_data,
            col_data=col_data,
        )
    except Exception as exc:  # pylint: disable=broad-except
        errors_seen.append(f"dict[str, ndarray] + col_data → {exc!r}")

    raise TypeError(
        "Failed to construct SummarizedExperiment; tried variants: "
        + "; ".join(errors_seen)
    )


def _count_compat_keys(res: resource.R3Resource) -> tuple[str, str]:
    """Return the (family, feature_key) pair for a count-file resource.

    The returned values are used to check compatibility before stacking
    count matrices across resources. There are two families:

    * ``\"gene_or_exon\"`` for gene- or exon-level matrices.
    * ``\"junctions\"`` for junction-level matrices.

    Args:
      res: A resource expected to represent a count matrix.

    Returns:
      A tuple ``(family, feature_key)`` where ``family`` groups resource
      types at a coarse level and ``feature_key`` captures the specific
      feature space (such as genomic unit or junction subtype).

    Raises:
      ValueError: If ``res`` is not a recognized count-file type.
    """
    rtype = getattr(res.description, "resource_type", None)
    if rtype == "count_files_gene_or_exon":
        genomic_unit = getattr(res.description, "genomic_unit", None) or ""
        family = "gene_or_exon"
        feature_key = f"{family}:{genomic_unit}"
        return family, feature_key

    if rtype == "count_files_junctions":
        junction_type = getattr(res.description, "junction_type", None) or ""
        junction_ext = (
            getattr(res.description, "junction_file_extension", None) or ""
        )
        family = "junctions"
        feature_key = f"{family}:{junction_type}:{junction_ext}"
        return family, feature_key

    raise ValueError(
        "Resource is not a recognized count-file type for stacking: "
        f"{rtype!r}"
    )


@dataclasses.dataclass(slots=True)
class R3ResourceBundle:
    """Container for a set of :class:`recount3.resource.R3Resource` objects.

    Bundles act as the primary orchestration primitive in this package.
    They keep track of a collection of resources and provide helpers for
    loading, filtering, project-aware workflows, and high-level
    operations such as stacking matrices or building BiocPy objects.

    A bundle may optionally be associated with a single project identity
    via the ``organism``, ``data_source``, and ``project`` attributes. If
    multiple projects are combined into one bundle, these attributes are
    left as :data:`None`.

    Attributes:
      resources: The list of resources contained in the bundle.
      organism: Optional organism identifier (for example, ``"human"`` or
        ``"mouse"``) when the bundle is project-scoped.
      data_source: Optional data source name (for example, ``"sra"``,
        ``"gtex"``, or ``"tcga"``) when the bundle is project-scoped.
      project: Optional study or project identifier (for example,
        ``"SRP009615"``) when the bundle is project-scoped.
    """

    resources: list[resource.R3Resource] = dataclasses.field(
        default_factory=list
    )
    organism: Optional[str] = None
    data_source: Optional[str] = None
    project: Optional[str] = None

    # -------------------------------------------------------------------
    # Construction and basic mutators
    # -------------------------------------------------------------------

    @classmethod
    def discover(
        cls,
        *,
        organism: r3_types.StringOrIterable,
        data_source: r3_types.StringOrIterable,
        project: r3_types.StringOrIterable,
        genomic_units: tuple[str, ...] = ("gene", "exon"),
        annotations: r3_types.StringOrIterable = "default",
        junction_exts: tuple[str, ...] = ("MM",),
        junction_type: str = "ALL",
        include_metadata: bool = True,
        include_bigwig: bool = False,
        strict: bool = True,
        deduplicate: bool = True,
    ) -> R3ResourceBundle:
        """Discover resources for one or more projects and return a bundle.

        This is a generalized version of
        :meth:`recount3.project.R3Project.discover`. It can operate on a
        single ``(organism, data_source, project)`` triple or on the
        Cartesian product of multiple values for each identifier.

        When discovery spans more than one project, the returned bundle
        will contain resources from all projects and the bundle-level
        ``organism``, ``data_source``, and ``project`` attributes will be
        left as :data:`None` to avoid misrepresenting the identity.

        Args:
          organism: Single organism name or iterable of names.
          data_source: Single data source or iterable of data sources.
          project: Single project identifier or iterable of identifiers.
          genomic_units: Gene expression feature levels to include; for
            example, ``("gene", "exon")``.
          annotations: Annotation selection; either ``"default"``,
            ``"all"``, a single annotation code, or an iterable of codes
            or labels understood by
            :func:`recount3.search.annotation_ext`.
          junction_exts: Junction file extensions to include; typically
            ``("MM",)`` for junction counts, with ``"RR"`` and/or
            ``"ID"`` added for coordinates or IDs.
          junction_type: Junction type selector, such as ``"ALL"``.
          include_metadata: Whether to include the 5 project metadata
            tables in the result.
          include_bigwig: Whether to include per-sample BigWig coverage
            resources.
          strict: If :data:`True`, propagate errors for invalid inputs or
            missing projects. If :data:`False`, attempts that fail
            validation are skipped.
          deduplicate: If :data:`True`, remove duplicated resources across
            discovered projects.

        Returns:
          A new :class:`R3ResourceBundle` populated with discovered
          resources. When discovery covers exactly one
          ``(organism, data_source, project)`` triple, the resulting
          bundle's ``organism``, ``data_source``, and ``project``
          attributes are set accordingly.

        Raises:
          ValueError: If all of ``organism``, ``data_source``, or
            ``project`` evaluate to an empty collection after
            normalization.
          recount3.errors.ConfigurationError: If the underlying search
            logic reports configuration problems.
        """
        organism_values = search._as_tuple(organism)
        data_source_values = search._as_tuple(data_source)
        project_values = search._as_tuple(project)

        if not organism_values or not data_source_values or not project_values:
            raise ValueError(
                "Arguments 'organism', 'data_source', and 'project' must "
                "not be empty."
            )

        resources_list: list[resource.R3Resource] = []
        identities: set[tuple[str, str, str]] = set()

        for org in organism_values:
            for src in data_source_values:
                for proj in project_values:
                    discovered = search.search_project_all(
                        organism=org,
                        data_source=src,
                        project=proj,
                        genomic_units=genomic_units,
                        annotations=annotations,
                        junction_file_extension=junction_exts,
                        junction_type=junction_type,
                        include_metadata=include_metadata,
                        include_bigwig=include_bigwig,
                        strict=strict,
                        deduplicate=deduplicate,
                    )
                    if discovered:
                        resources_list.extend(discovered)
                        identities.add((org, src, proj))

        if deduplicate:
            unique: dict[tuple[str, str], resource.R3Resource] = {}
            for res in resources_list:
                key = (res.url or "", res.description.url_path())
                if key not in unique:
                    unique[key] = res
            resources_list = list(unique.values())

        if len(identities) == 1:
            org, src, proj = next(iter(identities))
            return cls(
                resources=resources_list,
                organism=org,
                data_source=src,
                project=proj,
            )

        # Multi-project bundles do not advertise a single project identity.
        return cls(resources=resources_list)

    def add(self, res: resource.R3Resource) -> None:
        """Add a resource to the bundle.

        Args:
          res: The resource to append to :attr:`resources`.
        """
        self.resources.append(res)

    def extend(
        self,
        resources_iter: Iterable[resource.R3Resource],
    ) -> None:
        """Extend the bundle with additional resources.

        Args:
          resources_iter: Iterable of resources to add to the bundle.
        """
        self.resources.extend(resources_iter)

    # -------------------------------------------------------------------
    # Loading and iteration
    # -------------------------------------------------------------------

    def load(
        self,
        *,
        strict: bool = True,
        force: bool = False,
    ) -> R3ResourceBundle:
        """Load all resources and cache their data on each instance.

        This method iterates over :attr:`resources` and calls
        :meth:`recount3.resource.R3Resource.load` on each one.

        Args:
          strict: If :data:`True`, stop at the first exception and
            re-raise it. If :data:`False`, skip resources that fail to
            load.
          force: If :data:`True`, force a reload even when data is
            already cached on the resource.

        Returns:
          This :class:`R3ResourceBundle` instance, to enable chaining.
        """
        for res in self.resources:
            try:
                res.load(force=force)
            except Exception:  # pylint: disable=broad-except
                if strict:
                    raise
        return self

    def iter_loaded(
        self,
        *,
        resource_type: Optional[str] = None,
        autoload: bool = False,
    ) -> Iterator[tuple[resource.R3Resource, Any]]:
        """Yield ``(resource, data)`` pairs for resources with loaded data.

        When ``autoload`` is :data:`True`, resources that have not yet
        been loaded are passed through :meth:`R3Resource.load` before
        yielding.

        Args:
          resource_type: Optional resource-type filter applied to
            ``res.description.resource_type``.
          autoload: If :data:`True`, automatically load resources that
            have not yet been loaded.

        Yields:
          Tuples of ``(resource, loaded_data)`` for each resource that
          matches the optional ``resource_type`` filter and either
          already has cached data or can be loaded successfully.
        """
        for res in self.resources:
            if resource_type is not None:
                rtype = getattr(res.description, "resource_type", None)
                if rtype != resource_type:
                    continue

            if not res.is_loaded():
                if autoload:
                    try:
                        res.load()
                    except Exception:  # pylint: disable=broad-except
                        continue
                else:
                    continue

            obj = res.get_loaded()
            if obj is not None:
                yield res, obj

    def iter_bigwig(
        self,
        *,
        autoload: bool = True,
    ) -> Iterator[tuple[resource.R3Resource, _bigwig.BigWigFile]]:
        """Yield ``(resource, bigwig)`` pairs for BigWig resources.

        Args:
          autoload: If :data:`True`, automatically load BigWig resources
            that have not yet been loaded.

        Yields:
          Pairs of :class:`recount3.resource.R3Resource` and
          :class:`recount3._bigwig.BigWigFile` objects.
        """
        for res, obj in self.iter_loaded(
            resource_type="bigwig_files",
            autoload=autoload,
        ):
            if isinstance(obj, _bigwig.BigWigFile):
                yield res, obj

    def get_loaded(
        self,
        *,
        resource_type: Optional[str] = None,
        autoload: bool = False,
    ) -> list[Any]:
        """Return loaded data objects for resources in the bundle.

        Args:
          resource_type: Optional resource-type filter applied to
            ``res.description.resource_type``.
          autoload: If :data:`True`, automatically load any resources
            that are not yet loaded.

        Returns:
          A list of loaded data objects corresponding to resources in the
          bundle.
        """
        return [
            obj
            for _, obj in self.iter_loaded(
                resource_type=resource_type,
                autoload=autoload,
            )
        ]

    # -------------------------------------------------------------------
    # Filtering and predicates
    # -------------------------------------------------------------------

    def filter(
        self,
        *,
        resource_type: Optional[r3_types.FieldSpec] = None,
        organism: Optional[r3_types.FieldSpec] = None,
        data_source: Optional[r3_types.FieldSpec] = None,
        genomic_unit: Optional[r3_types.FieldSpec] = None,
        project: Optional[r3_types.FieldSpec] = None,
        sample: Optional[r3_types.FieldSpec] = None,
        table_name: Optional[r3_types.FieldSpec] = None,
        junction_type: Optional[r3_types.FieldSpec] = None,
        annotation_file_extension: Optional[r3_types.FieldSpec] = None,
        junction_file_extension: Optional[r3_types.FieldSpec] = None,
        predicate: Optional[Callable[[resource.R3Resource], bool]] = None,
        invert: bool = False,
    ) -> R3ResourceBundle:
        """Return a new bundle containing resources that match criteria.

        Each keyword argument corresponds to an attribute on
        :class:`recount3._descriptions.R3ResourceDescription`. Values are
        interpreted using :func:`recount3.search._match_spec`, allowing
        simple values, iterables of values, or callables.

        Args:
          resource_type: Resource type filter.
          organism: Organism filter.
          data_source: Data source filter.
          genomic_unit: Genomic unit filter.
          project: Project identifier filter.
          sample: Sample identifier filter.
          table_name: Metadata table name filter.
          junction_type: Junction type filter.
          annotation_file_extension: Annotation code filter.
          junction_file_extension: Junction extension filter.
          predicate: Optional callback that receives each resource and
            returns :data:`True` if it should be kept.
          invert: If :data:`True`, invert the final match decision.

        Returns:
          A new :class:`R3ResourceBundle` containing only the resources
          that match all supplied filters and the optional predicate.
        """
        field_specs: dict[str, r3_types.FieldSpec] = {
            "resource_type": resource_type,
            "organism": organism,
            "data_source": data_source,
            "genomic_unit": genomic_unit,
            "project": project,
            "sample": sample,
            "table_name": table_name,
            "junction_type": junction_type,
            "annotation_file_extension": annotation_file_extension,
            "junction_file_extension": junction_file_extension,
        }
        field_specs = {
            key: value for key, value in field_specs.items() if value is not None
        }

        selected: list[resource.R3Resource] = []
        for res in self.resources:
            desc = res.description
            fields_ok = all(
                search._match_spec(
                    getattr(desc, name, None),
                    spec,
                )  # type: ignore[arg-type]
                for name, spec in field_specs.items()
            )

            predicate_ok = True
            if predicate is not None:
                try:
                    predicate_ok = bool(predicate(res))
                except Exception:  # pylint: disable=broad-except
                    predicate_ok = False

            match = fields_ok and predicate_ok
            if invert:
                match = not match

            if match:
                selected.append(res)

        return R3ResourceBundle(
            resources=selected,
            organism=self.organism,
            data_source=self.data_source,
            project=self.project,
        )

    def only_counts(self) -> R3ResourceBundle:
        """Return a bundle restricted to gene/exon or junction count files.

        Returns:
          A new :class:`R3ResourceBundle` containing only resources whose
          ``resource_type`` is ``"count_files_gene_or_exon"`` or
          ``"count_files_junctions"``.
        """
        return self.filter(
            resource_type=("count_files_gene_or_exon", "count_files_junctions")
        )

    def only_metadata(self) -> R3ResourceBundle:
        """Return a bundle restricted to metadata resources.

        Returns:
          A new :class:`R3ResourceBundle` containing only resources whose
          ``resource_type`` is ``"metadata_files"``.
        """
        return self.filter(resource_type="metadata_files")

    def exclude_metadata(self) -> R3ResourceBundle:
        """Return a bundle with metadata resources removed.

        Returns:
          A new :class:`R3ResourceBundle` that excludes resources whose
          ``resource_type`` is ``"metadata_files"``.
        """
        return self.filter(resource_type="metadata_files", invert=True)

    def where(
        self,
        predicate: Callable[[resource.R3Resource], bool],
    ) -> R3ResourceBundle:
        """Predicate-based helper that forwards to :meth:`filter`.

        Args:
          predicate: Function that receives each resource and returns
            :data:`True` if it should be retained in the result.

        Returns:
          A new :class:`R3ResourceBundle` with only resources for which
          ``predicate`` returned :data:`True`.
        """
        return self.filter(predicate=predicate)

    # Convenience aliases that mirror the former R3Project API.

    def counts(self) -> R3ResourceBundle:
        """Return a sub-bundle containing only count-file resources.

        This is a convenience alias for :meth:`only_counts` maintained
        for API continuity with :class:`recount3.project.R3Project`.
        """
        return self.only_counts()

    def metadata(self) -> R3ResourceBundle:
        """Return a sub-bundle containing only metadata resources.

        This is a convenience alias for :meth:`only_metadata` maintained
        for compatibility with :class:`recount3.project.R3Project`.
        """
        return self.only_metadata()

    def bigwigs(self) -> R3ResourceBundle:
        """Return a sub-bundle containing only BigWig resources.

        Returns:
          A new :class:`R3ResourceBundle` containing only resources whose
          type is ``"bigwig_files"``.
        """
        return self.filter(resource_type="bigwig_files")

    # -------------------------------------------------------------------
    # Project identity helpers
    # -------------------------------------------------------------------

    def _resolve_project_identity(
        self,
        organism: Optional[str],
        data_source: Optional[str],
        project: Optional[str],
    ) -> tuple[str, str, str]:
        """Resolve project identity from explicit arguments or attributes.

        Args:
          organism: Optional organism override.
          data_source: Optional data source override.
          project: Optional project override.

        Returns:
          A tuple ``(organism, data_source, project)``.

        Raises:
          ValueError: If any component of the identity is missing or if
            overrides conflict with the bundle's stored identity.
        """
        resolved_organism = organism or self.organism
        resolved_data_source = data_source or self.data_source
        resolved_project = project or self.project

        if (
            not resolved_organism
            or not resolved_data_source
            or not resolved_project
        ):
            raise ValueError(
                "Project identity is incomplete. Provide explicit values for "
                "'organism', 'data_source', and 'project', or construct the "
                "bundle via 'R3ResourceBundle.discover' so the identity is "
                "recorded."
            )

        if (
            self.organism is not None
            and organism is not None
            and organism != self.organism
        ):
            raise ValueError(
                "Explicit 'organism' does not match the bundle's stored "
                f"organism: {organism!r} != {self.organism!r}."
            )

        if (
            self.data_source is not None
            and data_source is not None
            and data_source != self.data_source
        ):
            raise ValueError(
                "Explicit 'data_source' does not match the bundle's stored "
                f"data_source: {data_source!r} != {self.data_source!r}."
            )

        if (
            self.project is not None
            and project is not None
            and project != self.project
        ):
            raise ValueError(
                "Explicit 'project' does not match the bundle's stored "
                f"project: {project!r} != {self.project!r}."
            )

        return resolved_organism, resolved_data_source, resolved_project

    def samples(
        self,
        *,
        organism: Optional[str] = None,
        data_source: Optional[str] = None,
        project: Optional[str] = None,
    ) -> list[str]:
        """Return the list of sample identifiers associated with a project.

        The default behavior uses the bundle's stored project identity, as
        recorded by :meth:`discover`. Explicit keyword arguments can be
        provided to override or define the identity when the bundle was
        not created by :meth:`discover`.

        Args:
          organism: Optional organism identifier override.
          data_source: Optional data source override.
          project: Optional project identifier override.

        Returns:
          A sorted list of sample identifiers for the resolved project.

        Raises:
          ValueError: If the project cannot be resolved or validated.
        """
        org, src, proj = self._resolve_project_identity(
            organism=organism,
            data_source=data_source,
            project=project,
        )
        return search.samples_for_project(
            organism=org,
            data_source=src,
            project=proj,
        )

    # -------------------------------------------------------------------
    # Matrix stacking and BiocPy integration
    # -------------------------------------------------------------------

    def stack_count_matrices(
        self,
        *,
        join: str = "inner",
        axis: int = 1,
        verify_integrity: bool = False,
        autoload: bool = True,
        compat: r3_types.CompatibilityMode = "family",
    ) -> pd.DataFrame:
        """Concatenate count matrices (gene/exon or junction) as DataFrames.

        Args:
          join: Join policy passed to :func:`pandas.concat`.
          axis: Concatenation axis passed to :func:`pandas.concat`.
          verify_integrity: If :data:`True`, raise when labels are not
            unique along the concatenation axis.
          autoload: If :data:`True`, automatically load resources prior to
            concatenation.
          compat: Compatibility mode. ``"family"`` enforces that all
            inputs come from the same high-level family (gene/exon or
            junction), while ``"feature"`` enforces an identical feature
            space (for example, same genomic unit and junction subtype).

        Returns:
          A :class:`pandas.DataFrame` containing the concatenated count
          matrices.

        Raises:
          recount3.errors.CompatibilityError: If incompatible count
            resources are mixed in a way that violates ``compat``.
          TypeError: If a loaded object is not a
            :class:`pandas.DataFrame`.
          ValueError: If no applicable resources are present or if no
            loaded count matrices are found.
        """
        wanted = {
            "count_files_gene_or_exon",
            "count_files_junctions",
        }
        count_resources = [
            res
            for res in self.resources
            if getattr(res.description, "resource_type", None) in wanted
        ]
        if not count_resources:
            raise ValueError(
                "No count-file resources available to stack in this bundle."
            )

        families: set[str] = set()
        features: set[str] = set()
        family_counts: dict[str, int] = {}

        for res in count_resources:
            try:
                family, feature_key = _count_compat_keys(res)
            except ValueError:
                continue
            families.add(family)
            features.add(feature_key)
            family_counts[family] = family_counts.get(family, 0) + 1

        if compat == "family":
            if len(families) > 1:
                details = ", ".join(
                    f"{name}={count}"
                    for name, count in sorted(family_counts.items())
                )
                raise errors.CompatibilityError(
                    "Incompatible count families selected for stacking. "
                    f"Found families: {sorted(families)} ({details}). "
                    "Stack gene/exon with gene/exon, and junctions with "
                    "junctions. Hint: filter first, for example, "
                    'bundle.filter(resource_type="count_files_gene_or_exon") '
                    'or bundle.filter(resource_type="count_files_junctions").'
                )
        elif compat == "feature":
            if len(features) > 1:
                examples = ", ".join(sorted(features))
                raise errors.CompatibilityError(
                    "Feature-level incompatibility detected. All inputs must "
                    "share the same feature key (for example, gene vs exon; "
                    "junction subtype). Distinct feature keys observed: "
                    f"{examples}. Hint: filter by 'genomic_unit' for "
                    "gene/exon or by 'junction_type' / "
                    "'junction_file_extension' for junctions."
                )
        else:
            raise ValueError(f"Unknown compat mode: {compat!r}")

        data_frames: list[pd.DataFrame] = []
        for res, obj in self.iter_loaded(autoload=autoload):
            rtype = getattr(res.description, "resource_type", None)
            if rtype not in wanted:
                continue
            if not isinstance(obj, pd.DataFrame):
                raise TypeError(
                    f"Loaded object for resource {res.url!r} is not a "
                    "pandas.DataFrame instance."
                )
            data_frames.append(obj)

        if not data_frames:
            raise ValueError(
                "No loaded count matrices found. Try 'autoload=True' or call "
                "'bundle.load()' before stacking."
            )

        return pd.concat(
            data_frames,
            axis=axis,  # type: ignore[arg-type]
            join=join,  # type: ignore[arg-type]
            verify_integrity=verify_integrity,
        )

    def _stack_counts_for(
        self,
        *,
        genomic_unit: str,
        join: str = "inner",
        autoload: bool = True,
    ) -> pd.DataFrame:
        """Return a wide counts DataFrame for the requested feature family.

        This is a bundle-scoped helper used by the SummarizedExperiment
        builders. It enforces appropriate compatibility within the
        gene/exon family or the junctions family.

        Args:
          genomic_unit: One of ``"gene"``, ``"exon"``, or ``"junction"``.
          join: Join mode for concatenation (see :func:`pandas.concat`).
          autoload: If :data:`True`, load resources on demand.

        Returns:
          A wide matrix of shape ``(features, samples)`` as a
          :class:`pandas.DataFrame`.

        Raises:
          ValueError: If no compatible count resources exist or stacking
            fails.
        """
        load_errors: list[tuple[str, Exception]] = []

        if genomic_unit in {"gene", "exon"}:
            sel = (
                self.filter(resource_type="count_files_gene_or_exon")
                .filter(genomic_unit=genomic_unit)
            )

            if autoload:
                for res in sel.resources:
                    try:
                        res.load()
                    except Exception as exc:  # pylint: disable=broad-except
                        load_errors.append((res.url, exc))

            try:
                return sel.stack_count_matrices(
                    join=join,
                    axis=1,
                    autoload=False,
                    compat="feature",
                )
            except ValueError as exc:
                if load_errors:
                    first_url, first_exc = load_errors[0]
                    raise ValueError(
                        "Failed to load any gene/exon count matrix. "
                        f"First error: {first_exc!r} (while loading "
                        f"{first_url})."
                    ) from exc
                raise

        sel = self.filter(resource_type="count_files_junctions")
        if autoload:
            for res in sel.resources:
                try:
                    res.load()
                except Exception as exc:  # pylint: disable=broad-except
                    load_errors.append((res.url, exc))

        try:
            return sel.stack_count_matrices(
                join=join,
                axis=1,
                autoload=False,
                compat="family",
            )
        except ValueError as exc:
            if load_errors:
                first_url, first_exc = load_errors[0]
                raise ValueError(
                    "Failed to load any junction count matrix. "
                    f"First error: {first_exc!r} (while loading {first_url})."
                ) from exc
            raise

    def _normalize_sample_metadata(
        self,
        *,
        sample_ids: Sequence[str],
        autoload: bool = True,
    ) -> pd.DataFrame:
        """Merge available metadata tables and align rows to samples.

        1. Load each per-project metadata table.
        2. Standardize column names and key fields.
        3. Namespace non-key columns with the metadata table name.
        4. Outer-merge all tables on ``rail_id``, ``external_id``, and
            ``study``.
        5. Align the merged rows to ``sample_ids``.

        Args:
            sample_ids: Column names of the counts matrix to align to.
            autoload: If True, load metadata resources on demand.

        Returns:
            A DataFrame indexed by ``sample_ids`` containing merged, namespaced
            metadata.

        The returned DataFrame columns are guaranteed to be unique.

        Raises:
            ValueError: If metadata cannot be aligned to the requested samples.
        """
        sample_ids = [str(s) for s in sample_ids]
        frames: list[pd.DataFrame] = []
        provenance: dict[str, tuple[str, str]] = {}

        for res, obj in self.only_metadata().iter_loaded(autoload=autoload):
            if not isinstance(obj, pd.DataFrame):
                continue
            origin = _metadata_origin(res)
            standardized = _standardize_metadata_frame(obj)
            namespaced, prov = _namespace_metadata_columns(
                standardized,
                origin=origin,
            )
            frames.append(namespaced)
            provenance.update(prov)

        if not frames:
            out = pd.DataFrame(
                index=sample_ids,
                data={"external_id": sample_ids},
            )
            out.index.name = None
            return _ensure_unique_columns(out, empty_prefix="col")

        merged = _outer_merge_metadata_frames(frames)

        align_key = _choose_alignment_key(sample_ids=sample_ids, merged=merged)
        merged = _collapse_rows_by_key(merged, key=align_key)

        if align_key not in merged.columns:
            raise ValueError(
                f"Cannot align metadata: missing alignment key {align_key!r}."
            )

        aligned = merged.set_index(align_key, drop=False).reindex(sample_ids)
        aligned.index.name = None

        # Always expose external_id as the canonical user-facing sample ID.
        if "external_id" not in aligned.columns:
            aligned["external_id"] = sample_ids
        elif aligned["external_id"].isna().any():
            aligned["external_id"] = aligned["external_id"].fillna(
                pd.Series(sample_ids, index=aligned.index)
            )

        aligned = _ensure_unique_columns(aligned, empty_prefix="col")
        aligned.attrs["recount3_metadata_provenance"] = provenance
        return aligned

    def to_summarized_experiment(
        self,
        *,
        genomic_unit: str,
        annotation_file_extension: Optional[str] = None,
        assay_name: str = "counts",
        join: str = "inner",
        autoload: bool = True,
    ) -> SummarizedExperiment:
        """Build a BiocPy :class:`SummarizedExperiment` from this bundle.

        This method stacks compatible count matrices, merges available
        sample metadata, and constructs a BiocPy
        :class:`SummarizedExperiment` using a compatibility-aware
        constructor that supports multiple versions of the
        :mod:`summarizedexperiment` package. 

        Args:
          genomic_unit: Genomic unit to summarize, such as ``"gene"``,
            ``"exon"``, or ``"junction"``.
          annotation_file_extension: Optional annotation code for gene or
            exon summarizations (for example, ``"G026"``). When provided
            and ``genomic_unit`` is gene or exon, only count resources
            with matching annotation are used.
          assay_name: Name assigned to the count assay within the
            :class:`SummarizedExperiment`.
          join: Join policy used when concatenating counts across
            resources.
          autoload: If :data:`True`, load resources when needed.

        Returns:
          A BiocPy :class:`SummarizedExperiment` instance.

        Raises:
          ImportError: If BiocPy packages are not installed.
          ValueError: If no counts are found or shapes are inconsistent.
          TypeError: If the underlying
            :class:`SummarizedExperiment` constructor rejects all
            compatibility variants.
        """
        _require_biocpy()

        working = self
        if genomic_unit in {"gene", "exon"} and annotation_file_extension:
            working = (
                self.filter(resource_type="count_files_gene_or_exon")
                .filter(
                    genomic_unit=genomic_unit,
                    annotation_file_extension=annotation_file_extension,
                )
            )

        counts_df = working._stack_counts_for(
            genomic_unit=genomic_unit,
            join=join,
            autoload=autoload,
        )

        counts_df.index = [str(i) for i in counts_df.index]
        counts_df.columns = [str(c) for c in counts_df.columns]

        col_df = self._normalize_sample_metadata(
            sample_ids=list(counts_df.columns),
            autoload=autoload,
        )
        counts_df, col_df = _maybe_relabel_counts_columns_to_external_id(
            counts_df,
            col_df,
        )

        row_df = pd.DataFrame(
            {"feature_id": [str(idx) for idx in counts_df.index]}
        )

        return _construct_se_compat(
            counts_df=counts_df,
            row_df=row_df,
            col_df=col_df,
            assay_name=assay_name,
        )

    def to_ranged_summarized_experiment(
        self,
        *,
        genomic_unit: str,
        annotation_file_extension: Optional[str] = None,
        junction_rr_preferred: bool = True,
        assay_name: str = "counts",
        join: str = "inner",
        autoload: bool = True,
        allow_fallback_to_se: bool = False,
    ) -> RangedSummarizedExperiment | SummarizedExperiment:
        """Build a BiocPy :class:`RangedSummarizedExperiment` when possible.

        For ``"gene"`` and ``"exon"`` genomic units, row ranges are
        derived from a matching GTF(.gz) annotation resource. For
        ``"junction"``, this method prefers an RR table (junction
        coordinates) when available.

        When row ranges cannot be resolved and
        ``allow_fallback_to_se`` is :data:`True`, a plain
        :class:`SummarizedExperiment` is returned instead.

        Args:
          genomic_unit: One of ``"gene"``, ``"exon"``, or ``"junction"``.
          annotation_file_extension: Annotation code for gene/exon
            assays, if desired.
          junction_rr_preferred: If :data:`True`, prefer RR junction
            files for coordinate definitions when they are available.
          assay_name: Name assigned to the assay.
          join: Join policy across projects when stacking.
          autoload: If :data:`True`, load resources transparently.
          allow_fallback_to_se: If :data:`True`, construct a plain
            :class:`SummarizedExperiment` when genomic ranges cannot be
            derived for the requested combination.

        Returns:
          A :class:`RangedSummarizedExperiment` instance, or a plain
          :class:`SummarizedExperiment` when ``allow_fallback_to_se`` is
          :data:`True` and ranges are unavailable.

        Raises:
          ImportError: If BiocPy packages are not installed.
          ValueError: If counts are missing or ranges cannot be
            determined and ``allow_fallback_to_se`` is :data:`False`.
          TypeError: If the
            :class:`RangedSummarizedExperiment` constructor rejects all
            compatibility variants.
        """
        BiocFrameCls, GenomicRangesCls, SummarizedExperimentCls, RangedSummarizedExperimentCls = _require_biocpy()

        working = self
        counts_df = working._stack_counts_for(
            genomic_unit=genomic_unit,
            join=join,
            autoload=autoload,
        )

        counts_df.index = [str(i) for i in counts_df.index]
        counts_df.columns = [str(c) for c in counts_df.columns]

        col_df = self._normalize_sample_metadata(
            sample_ids=list(counts_df.columns),
            autoload=autoload,
        )
        counts_df, col_df = _maybe_relabel_counts_columns_to_external_id(
            counts_df,
            col_df,
        )

        row_data_df = pd.DataFrame(
            {"feature_id": [str(idx) for idx in counts_df.index]}
        )
        feature_ids: list[str] = list(row_data_df["feature_id"])

        ranges_df: Optional[pd.DataFrame] = None

        if genomic_unit in {"gene", "exon"}:
            ann_bundle = self.filter(resource_type="annotations")
            if annotation_file_extension:
                ann_bundle = ann_bundle.filter(
                    annotation_file_extension=annotation_file_extension
                )

            gtf_res = next(iter(ann_bundle.resources), None)
            if gtf_res is not None:
                try:
                    if autoload:
                        gtf_res.download(path=None, cache_mode="enable")
                    gtf = _read_gtf_dataframe(gtf_res)
                    feature_kind = "gene" if genomic_unit == "gene" else "exon"
                    ranges = _ranges_from_gtf(
                        gtf,
                        feature_kind=feature_kind,
                    )
                    idxed = ranges.set_index("feature_id").reindex(feature_ids)
                    ranges_df = idxed[
                        ["seqnames", "starts", "ends", "strand"]
                    ].reset_index(drop=True)

                    enrich_cols = [
                        col
                        for col in ranges.columns
                        if col
                        not in {
                            "seqnames",
                            "starts",
                            "ends",
                            "strand",
                            "feature_id",
                        }
                    ]
                    if enrich_cols:
                        anno_extra = ranges.set_index("feature_id")[
                            enrich_cols
                        ].reindex(feature_ids)
                        row_data_df = pd.concat(
                            [
                                row_data_df.reset_index(drop=True),
                                anno_extra.reset_index(drop=True),
                            ],
                            axis=1,
                        )
                except Exception as exc:  # pylint: disable=broad-except
                    logging.warning(
                        "Falling back: failed to parse GTF for %s ranges "
                        "(reason: %r).",
                        genomic_unit,
                        exc,
                    )

        elif genomic_unit == "junction" and junction_rr_preferred:
            rr_res = next(
                (
                    res
                    for res in self.filter(
                        resource_type="count_files_junctions"
                    ).resources
                    if getattr(
                        res.description,
                        "junction_file_extension",
                        "",
                    )
                    == "RR"
                ),
                None,
            )
            if rr_res is not None:
                try:
                    rr = _read_rr_table(rr_res)
                    candidates = ("junction_id", "jxn_id", "jid", "id", "name")
                    id_col = next(
                        (col for col in candidates if col in rr.columns),
                        None,
                    )
                    if id_col is not None:
                        rr = rr.set_index(rr[id_col].astype(str))
                        aligned = rr.reindex(feature_ids)
                        std = aligned.rename(
                            columns={
                                "chrom": "seqnames",
                                "chr": "seqnames",
                                "start": "starts",
                                "end": "ends",
                            }
                        )
                        if not {
                            "seqnames",
                            "starts",
                            "ends",
                        }.issubset(std.columns):
                            raise ValueError(
                                "RR file is missing coordinate columns."
                            )
                        if "strand" not in std.columns:
                            std["strand"] = "*"
                        ranges_df = std[
                            ["seqnames", "starts", "ends", "strand"]
                        ].reset_index(drop=True)
                except Exception as exc:  # pylint: disable=broad-except
                    logging.warning(
                        "Falling back: failed to use RR for junction ranges "
                        "(reason: %r).",
                        exc,
                    )

        if ranges_df is None:
            if allow_fallback_to_se:
                return self.to_summarized_experiment(
                    genomic_unit=genomic_unit,
                    annotation_file_extension=annotation_file_extension,
                    assay_name=assay_name,
                    join=join,
                    autoload=autoload,
                )
            raise ValueError(
                "Could not derive genomic ranges for the requested object. "
                "Pass allow_fallback_to_se=True to receive a plain "
                "SummarizedExperiment."
            )

        row_data = row_data_df.copy()
        row_data.index = counts_df.index

        col_data = col_df.copy()
        col_data.index = counts_df.columns

        row_data = _ensure_unique_columns(row_data, empty_prefix="row")
        col_data = _ensure_unique_columns(col_data, empty_prefix="col")

        counts_numeric = counts_df.apply(pd.to_numeric, errors="coerce").fillna(0)

        gr = _to_genomic_ranges(ranges_df)

        errors_seen: list[str] = []

        # Variant 1: dict[str, ndarray] + column_data
        try:
            return RangedSummarizedExperimentCls(
                assays={assay_name: counts_numeric.to_numpy(copy=False)},
                row_data=row_data,
                row_ranges=gr,
                column_data=col_data,
            )
        except Exception as exc:
            errors_seen.append(
                "dict[str, ndarray] + column_data → "
                f"{exc!r}"
            )

        # Variant 2: list[ndarray] + assay_names + column_data
        try:
            return RangedSummarizedExperimentCls(
                assays=[counts_numeric.to_numpy(copy=False)],
                assay_names=[assay_name],
                row_data=row_data,
                row_ranges=gr,
                column_data=col_data,
            )
        except Exception as exc:
            errors_seen.append(
                "list[ndarray] + assay_names + column_data "
                f"→ {exc!r}"
            )

        # Variant 3: swap keyword to col_data
        try:
            return RangedSummarizedExperimentCls(
                assays={assay_name: counts_numeric.to_numpy(copy=False)},
                row_data=row_data,
                row_ranges=gr,
                col_data=col_data,
            )
        except Exception as exc:
            errors_seen.append(
                "dict[str, ndarray] + col_data → "
                f"{exc!r}"
            )

        raise TypeError(
            "Failed to construct RangedSummarizedExperiment; tried variants: "
            + "; ".join(errors_seen)
        )

    # -------------------------------------------------------------------
    # One-shot materialization
    # -------------------------------------------------------------------

    def download(
        self,
        *,
        dest: str = ".",
        overwrite: bool = False,
        cache: r3_types.CacheMode = "enable",
    ) -> None:
        """Download all resources in the bundle to a local destination.

        This method is a convenience wrapper around
        :meth:`recount3.resource.R3Resource.download` for each contained
        resource. For more advanced workflows (event logs, streaming to a
        ZIP archive, and so on) prefer the command-line interface,
        ``recount3 download``.

        Args:
          dest: Destination directory or ``.zip`` path. When a directory
            is provided, each resource is materialized as a separate file
            under that directory. When a path ending in ``.zip`` is
            provided, resources are written into that archive.
          overwrite: If :data:`True`, allow overwriting existing files in
            directory mode.
          cache: Cache behavior: ``"enable"``, ``"disable"``, or
            ``"update"`` as defined by :class:`recount3.types.CacheMode`.

        Raises:
          ValueError: Propagated from underlying resource download
            failures, for example when an unsupported cache mode is
            selected.
        """
        for res in self.resources:
            res.download(path=dest, cache_mode=cache, overwrite=overwrite)
