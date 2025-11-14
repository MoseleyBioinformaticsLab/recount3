"""High-level builders for SummarizedExperiment objects.

This module converts one or more recount3 resources into either a
:class:`summarizedexperiment.SummarizedExperiment` (SE) or a
:class:`summarizedexperiment.RangedSummarizedExperiment` (RSE).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import gzip
import io
import logging
import re

import pandas as pd

from .bundle import R3ResourceBundle
from .resource import R3Resource
from .project import R3Project
from .search import annotation_ext


# Attempt optional imports (pattern matches _bigwig.py).
try:  # pragma: no cover - optional, exercised only if installed
    from biocframe import BiocFrame  # type: ignore[import-not-found]
    from genomicranges import GenomicRanges  # type: ignore[import-not-found]
    from summarizedexperiment import (  # type: ignore[import-not-found]
        SummarizedExperiment,
        RangedSummarizedExperiment,
    )
except Exception:  # pragma: no cover
    BiocFrame = None
    GenomicRanges = None
    SummarizedExperiment = None
    RangedSummarizedExperiment = None


# ----------------------------- helpers -------------------------------------

def _ensure_unique_columns(df: pd.DataFrame, *, empty_prefix: str = "col") -> pd.DataFrame:
    """Return a copy with unique, non-empty string column names.
 
    Many recount3 metadata tables share common column names (e.g., ``external_id``),
    and concatenation can introduce duplicates. BiocFrame drops duplicated names
    which then breaks downstream validation. This helper ensures that:
 
    * All column names are strings (``None`` becomes an empty string).
    * Empty names are replaced by ``{empty_prefix}``.
    * Duplicates are suffixed as ``name__2``, ``name__3``, ...
 
    Args:
      df: Input DataFrame whose columns may contain duplicates or empties.
      empty_prefix: Base name used when an empty/None column name is encountered.
 
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

def _require_biocpy() -> None:
    """Ensure BiocPy packages are importable.

    Raises:
      ImportError: If any required BiocPy package is missing.
    """
    if any(x is None for x in (BiocFrame, GenomicRanges, SummarizedExperiment)):
        raise ImportError(
            "This feature requires BiocPy packages. Install with:\n"
            "  pip install summarizedexperiment biocframe genomicranges\n"
            "or via conda (conda-forge)."
        )


def _stack_counts_for(  #TODO: Should be a method of Bundle?
    bundle: R3ResourceBundle,
    *,
    genomic_unit: str,
    join: str = "inner",
    autoload: bool = True,
) -> pd.DataFrame:
    """Return a wide counts DataFrame for the requested feature family.

    Args:
      bundle: Resource bundle containing counts and metadata.
      genomic_unit: One of 'gene', 'exon', or 'junction'.
      join: Join mode for concatenation (see pandas.concat).
      autoload: If True, load resources on-demand.

    Returns:
      A wide matrix of shape (features, samples) as a pandas DataFrame.

    Raises:
      ValueError: If no compatible counts exist.
    """
    load_errors: list[tuple[str, Exception]] = []

    if genomic_unit in {"gene", "exon"}:
        sel = bundle.filter(resource_type="count_files_gene_or_exon").filter(
            genomic_unit=genomic_unit
        )
        # Eager-load for clearer errors because iter_loaded() suppresses them.
        if autoload:
            for res in sel.resources:
                try:
                    res.load()
                except Exception as exc:  # pylint: disable=broad-except
                    load_errors.append((res.url, exc))
        # Require feature-level compatibility within the family (gene vs exon).
        try:
            df = sel.stack_count_matrices(
                join=join, axis=1, autoload=False, compat="feature"
            )
            return df
        except ValueError as exc:
            if load_errors:
                first_url, first_exc = load_errors[0]
                raise ValueError(
                    "Failed to load any gene/exon count matrix. "
                    f"First error: {first_exc!r} (while loading {first_url})."
                ) from exc
            raise

    sel = bundle.filter(resource_type="count_files_junctions")
    if autoload:
        for res in sel.resources:
            try:
                res.load()
            except Exception as exc:  # pylint: disable=broad-except
                load_errors.append((res.url, exc))
    # Family-level compatibility is sufficient for junctions.
    try:
        df = sel.stack_count_matrices(
            join=join, axis=1, autoload=False, compat="family"
        )
        return df
    except ValueError as exc:
        if load_errors:
            first_url, first_exc = load_errors[0]
            raise ValueError(
                "Failed to load any junction count matrix. "
                f"First error: {first_exc!r} (while loading {first_url})."
            ) from exc
        raise


def _normalize_sample_metadata(
    bundle: R3ResourceBundle,
    *,
    sample_ids: Sequence[str],
    autoload: bool = True,
) -> pd.DataFrame:
    """Merge available metadata tables and align rows to samples.

    We align by trying multiple common sample keys present in recount3 tables:
    'external_id', 'rail_id', 'run', 'run_accession', 'sample', 'sample_id'.

    Args:
      bundle: Bundle with metadata resources.
      sample_ids: Column names of the counts matrix to align to.
      autoload: If True, load metadata resources as needed.

    Returns:
      A pandas DataFrame indexed by `sample_ids`. Missing values are NaN.
    """
    keys = (
        "external_id",
        "rail_id",
        "run",
        "run_accession",
        "sample",
        "sample_id",
    )

    parts: list[pd.DataFrame] = []
    for res, obj in bundle.only_metadata().iter_loaded(autoload=autoload):
        if not isinstance(obj, pd.DataFrame):
            continue
        df = obj.copy()
        # All-lowercase is common across recount3 helpers; normalize politely.
        df.columns = [str(c).strip() for c in df.columns]
        key = next((k for k in keys if k in df.columns), None)
        if not key:
            continue
        # Align by the discovered key.
        df = df.set_index(key)
        aligned = df.reindex(sample_ids)
        parts.append(aligned)

    if not parts:
        # Provide a minimal column_data when nothing is available.
        return pd.DataFrame(index=list(sample_ids), data={"external_id": sample_ids})

    out = pd.concat(parts, axis=1, join="outer")
    # Ensure index order/labels match exactly.
    out = out.reindex(sample_ids)
    out.index.name = None
    return out

def _expand_sra_attributes_df(
    col_df: pd.DataFrame,
    *,
    column: str = "sra.sample_attributes",
    prefix: str = "sra_attribute.",
) -> pd.DataFrame:
    """Return a copy of column metadata with SRA attributes expanded.

    recount3 stores SRA sample attributes in a single string column such as
    'sra.sample_attributes' using the encoding::

        "age;;67.78|biomaterial_provider;;LIBD|disease;;Control|..."

    This helper mirrors the behavior of the recount3 R function
    ``expand_sra_attributes()``: each entry is split on '|' into individual
    attributes, then on ';;' into ``key`` and ``value`` pairs. New columns
    named ``{prefix}{key}`` (with spaces in ``key`` replaced by '_') are
    appended to the metadata.

    Args:
      col_df: Column metadata DataFrame.
      column: Name of the column carrying the raw SRA attribute strings. If
        this column is missing, ``col_df`` is returned unchanged.
      prefix: Prefix to prepend to attribute names when forming new columns.

    Returns:
      A new DataFrame with the same index as ``col_df`` and additional
      columns containing one parsed SRA attribute per column.
    """
    if column not in col_df.columns:
        return col_df.copy()

    ser = col_df[column]
    rows: list[dict[str, object]] = []

    for value in ser:
        row_map: dict[str, object] = {}
        if isinstance(value, str) and value:
            parts = value.split("|")
            for entry in parts:
                if not entry:
                    continue
                if ";;" not in entry:
                    continue
                key, val = entry.split(";;", 1)
                key_clean = key.strip()
                if not key_clean:
                    continue
                col_name = f"{prefix}{key_clean.replace(' ', '_')}"
                # Last occurrence wins if the same key appears multiple times.
                row_map[col_name] = val.strip()
        rows.append(row_map)

    if not rows:
        return col_df.copy()

    attrs_df = pd.DataFrame(rows, index=col_df.index)
    # attrs_df may be empty (no valid attributes); concat handles that.
    merged = pd.concat([col_df, attrs_df], axis=1)
    return merged


@dataclass(slots=True)
class _JxnRange:
    """Lightweight holder for a junction range."""

    seqnames: str
    starts: int
    ends: int
    strand: str
    jid: str


def _read_rr_table(resource: R3Resource) -> pd.DataFrame:
    """Read an 'RR' junction coordinate table into a DataFrame.

    The RR file is a small TSV(.gz) table with one row per junction containing
    genomic coordinates and an identifier (often a 'junction_id' or 'jid').

    Args:
      resource: An RR-type junction resource.

    Returns:
      A pandas DataFrame containing RR columns.

    Raises:
      ValueError: If the resource is not RR-like or cannot be parsed.
    """
    # We route through R3Resource.load if available; if not recognized
    # (rare for RR), stream-parse the gzipped file.
    try:
        df = resource.load()
        if isinstance(df, pd.DataFrame):
            return df
    except Exception:  # Fall back to manual parsing below.
        pass

    # Manual path (transparent for .gz and plain text).
    # Use a cached path if possible; we rely on the internal helper for now.
    try:
        path = resource._cached_path()  # pylint: disable=protected-access
    except Exception as exc:
        raise ValueError(f"Cannot resolve cached RR path for: {resource.url}") from exc

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:  # type: ignore[arg-type]
        text = io.TextIOWrapper(fh, encoding="utf-8")
        df = pd.read_csv(text, sep="\t", comment="#")
    return df


def _ranges_from_gtf(
    gtf: pd.DataFrame, *, feature_kind: str
) -> pd.DataFrame:
    """Convert GTF rows for a feature type ('gene' or 'exon') into ranges.

    Args:
      gtf: A DataFrame with standard GTF columns plus an 'attributes' field.
      feature_kind: Either 'gene' or 'exon'.

    Returns:
      DataFrame with columns ['seqnames', 'starts', 'ends', 'strand', 'feature_id'].
    """
    # Filter to the requested feature lines and parse attributes for IDs.
    df = gtf.loc[gtf["feature"] == feature_kind].copy()
    if df.empty:
        return pd.DataFrame(
            columns=["seqnames", "starts", "ends", "strand", "feature_id"]
        )

    def _parse_attr(s: str) -> dict[str, str]:
        items = {}
        for part in re.split(r";\s*", s.strip().rstrip(";")):
            if not part:
                continue
            # GTF attributes: key "value"
            m = re.match(r'([^ ]+)\s+"([^"]+)"', part)
            if m:
                items[m.group(1)] = m.group(2)
        return items

    attrs = df["attributes"].astype(str).map(_parse_attr)
    if feature_kind == "gene":
        fid = [a.get("gene_id") for a in attrs]
    else:
        # Prefer exon_id, fall back to constructed ID if absent.
        fids = []
        for a, (_, row) in zip(attrs, df.iterrows()):
            eid = a.get("exon_id")
            if eid:
                fids.append(eid)
                continue
            # Fallback: chr:start-end:strand (unlikely to match counts index
            # in all cases but better than NaN).
            fids.append(
                f"{row['seqname']}:{row['start']}-{row['end']}:{row['strand']}"
            )
        fid = fids

    out = pd.DataFrame(
        {
            "seqnames": df["seqname"].astype(str).values,
            "starts": df["start"].astype(int).values,
            "ends": df["end"].astype(int).values,
            "strand": df["strand"].astype(str).values,
            "feature_id": fid,
        }
    )
    return out


def _read_gtf_dataframe(res: R3Resource) -> pd.DataFrame:
    """Read a GTF(.gz) into a lightly parsed DataFrame.

    Args:
      res: An annotations resource (.gtf or .gtf.gz).

    Returns:
      DataFrame with the standard 9 GTF columns, including 'attributes'.

    Raises:
      ValueError: If the file cannot be read.
    """
    try:
        path = res._cached_path()  # pylint: disable=protected-access
    except Exception as exc:
        raise ValueError(f"Cannot resolve cached GTF path for: {res.url}") from exc

    opener = gzip.open if str(path).endswith(".gz") else open
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
    with opener(path, "rb") as fh:  # type: ignore[arg-type]
        text = io.TextIOWrapper(fh, encoding="utf-8")
        df = pd.read_csv(
            text, sep="\t", comment="#", header=None, names=cols, dtype=str
        )
    # Normalize numeric columns.
    for c in ("start", "end"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


# --------------------------- public API -------------------------------------


def build_summarized_experiment(
    bundle: R3ResourceBundle,
    *,
    genomic_unit: str,
    annotation_file_extension: str | None = None,
    assay_name: str = "counts",
    join: str = "inner",
    autoload: bool = True,
) -> "SummarizedExperiment":
    """Create a SummarizedExperiment from a bundle.

    This stacks compatible count matrices, aligns sample metadata, and builds
    a BiocPy SummarizedExperiment using a compatibility-first constructor
    (handles multiple versions of the library).

    Args:
      bundle: Resources (counts + metadata) to assemble.
      genomic_unit: One of {'gene', 'exon', 'junction'}.
      annotation_file_extension: Optional annotation key for gene/exon.
      assay_name: Name for the count assay in the SE.
      join: Join policy across projects (pandas concat join).
      autoload: If True, load resources transparently.

    Returns:
      A :class:`summarizedexperiment.SummarizedExperiment`.

    Raises:
      ImportError: If BiocPy packages are not installed.
      ValueError: If no counts found or shapes do not align.
      TypeError: If all constructor variants are rejected by the library.
    """
    _require_biocpy()

    work = bundle
    if genomic_unit in {"gene", "exon"} and annotation_file_extension:
        work = work.filter(resource_type="count_files_gene_or_exon").filter(
            genomic_unit=genomic_unit,
            annotation_file_extension=annotation_file_extension,
        )
    counts_df = _stack_counts_for(
        work, genomic_unit=genomic_unit, join=join, autoload=autoload
    )

    # Align metadata to column order.
    sample_ids = list(map(str, counts_df.columns))
    col_df = _normalize_sample_metadata(bundle, sample_ids=sample_ids, autoload=autoload)

    # Minimal row_data: expose feature ids.
    row_df = pd.DataFrame({"feature_id": list(map(str, counts_df.index))})

    return _construct_se_compat(
        counts_df=counts_df, row_df=row_df, col_df=col_df, assay_name=assay_name
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
) -> "RangedSummarizedExperiment | SummarizedExperiment":
    """Create a RangedSummarizedExperiment when ranges can be resolved.

    For 'gene' and 'exon', derive row ranges from a matching GTF(.gz).
    For 'junction', prefer an RR table (junction coordinates). If ranges
    cannot be determined and `allow_fallback_to_se` is True, return a plain SE.

    Args:
      bundle: Resources for counts, metadata, and (optionally) annotations.
      genomic_unit: One of {'gene', 'exon', 'junction'}.
      annotation_file_extension: Annotation key for gene/exon, if desired.
      junction_rr_preferred: Prefer RR for junction coordinates when present.
      assay_name: Name for the count assay in the output.
      join: Join policy across projects when stacking.
      autoload: If True, load resources transparently.
      allow_fallback_to_se: Return a plain SE if ranges are unavailable.

    Returns:
      A :class:`summarizedexperiment.RangedSummarizedExperiment` or a plain SE.

    Raises:
      ImportError: If BiocPy packages are not installed.
      ValueError: If counts are missing or ranges cannot be determined and
        `allow_fallback_to_se` is False.
      TypeError: If all constructor variants are rejected by the library.
    """
    _require_biocpy()

    counts_df = _stack_counts_for(
        bundle, genomic_unit=genomic_unit, join=join, autoload=autoload
    )
    # Build metadata frames aligned to counts.
    sample_ids = list(map(str, counts_df.columns))
    col_df = _normalize_sample_metadata(
        bundle, sample_ids=sample_ids, autoload=autoload
    )
    row_data_df = pd.DataFrame({"feature_id": list(map(str, counts_df.index))})
    feature_ids: list[str] = list(row_data_df["feature_id"])

    # Resolve ranges.
    ranges_df: pd.DataFrame | None = None

    if genomic_unit in {"gene", "exon"}:
        ann = bundle.filter(resource_type="annotations")
        if annotation_file_extension:
            ann = ann.filter(annotation_file_extension=annotation_file_extension)
        gtf_res = next((r for r in ann.resources), None)
        if gtf_res:
            try:
                gtf = _read_gtf_dataframe(gtf_res)
                feat = "gene" if genomic_unit == "gene" else "exon"
                rd = _ranges_from_gtf(gtf, feature_kind=feat)
                idxed = rd.set_index("feature_id").reindex(feature_ids)
                ranges_df = idxed[["seqnames", "starts", "ends", "strand"]].reset_index(
                    drop=True
                )
                # Optionally enrich row_data with extra annotation columns.
                enrich_cols = [
                    c
                    for c in rd.columns
                    if c
                    not in {"seqnames", "starts", "ends", "strand", "feature_id"}
                ]
                if enrich_cols:
                    anno_extra = rd.set_index("feature_id")[enrich_cols].reindex(
                        feature_ids
                    )
                    row_data_df = pd.concat(
                        [
                            row_data_df.reset_index(drop=True),
                            anno_extra.reset_index(drop=True),
                        ],
                        axis=1,
                    )
            except Exception as exc:  # pylint: disable=broad-except
                logging.warning(
                    "Falling back: failed to parse GTF for %s ranges (reason: %r).",
                    genomic_unit,
                    exc,
                )

    elif genomic_unit == "junction" and junction_rr_preferred:
        rr_res = next(
            (
                r
                for r in bundle.filter(resource_type="count_files_junctions").resources
                if getattr(r.description, "junction_file_extension", "") == "RR"
            ),
            None,
        )
        if rr_res:
            try:
                rr = _read_rr_table(rr_res)
                candidates = ("junction_id", "jxn_id", "jid", "id", "name")
                id_col = next((c for c in candidates if c in rr.columns), None)
                if id_col:
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
                    if not {"seqnames", "starts", "ends"}.issubset(std.columns):
                        raise ValueError("RR file missing coordinate columns.")
                    if "strand" not in std.columns:
                        std["strand"] = "*"
                    ranges_df = std[
                        ["seqnames", "starts", "ends", "strand"]
                    ].reset_index(drop=True)
            except Exception as exc:  # pylint: disable=broad-except
                logging.warning(
                    "Falling back: failed to use RR for junction ranges (reason: %r).",
                    exc,
                )

    if ranges_df is None:
        if allow_fallback_to_se:
            return _construct_se_compat(
                counts_df=counts_df,
                row_df=row_data_df,
                col_df=col_df,
                assay_name=assay_name,
            )
        raise ValueError(
            "Could not derive genomic ranges for the requested object. "
            "Pass allow_fallback_to_se=True to receive a plain SE."
        )

    # Align indices and construct RSE with the same compatibility path.
    rd = row_data_df.copy()
    rd.index = counts_df.index
    cd = col_df.copy()
    cd.index = counts_df.columns
    counts_numeric = counts_df.apply(pd.to_numeric, errors="coerce").fillna(0)
    gr = GenomicRanges.from_pandas(ranges_df)

    # Try the common variants first; surface first error if none succeed.
    errors: list[str] = []

    # Variant 1: dict[str, ndarray] + column_data
    try:
        return RangedSummarizedExperiment(
            assays={assay_name: counts_numeric.to_numpy(copy=False)},
            row_data=rd,            # pandas; library will coerce if needed
            row_ranges=gr,
            column_data=cd,         # column_data (most common)
        )
    except Exception as exc:  # pylint: disable=broad-except
        errors.append(f"dict[str, ndarray] + column_data → {exc!r}")

    # Variant 2: list[ndarray] + assay_names + column_data
    try:
        return RangedSummarizedExperiment(
            assays=[counts_numeric.to_numpy(copy=False)],
            assay_names=[assay_name],
            row_data=rd,
            row_ranges=gr,
            column_data=cd,
        )
    except Exception as exc:  # pylint: disable=broad-except
        errors.append(f"list[ndarray]+assay_names + column_data → {exc!r}")

    # Variant 3: swap keyword to col_data
    try:
        return RangedSummarizedExperiment(
            assays={assay_name: counts_numeric.to_numpy(copy=False)},
            row_data=rd,
            row_ranges=gr,
            col_data=cd,
        )
    except Exception as exc:  # pylint: disable=broad-except
        errors.append(f"dict[str, ndarray] + col_data → {exc!r}")

    raise TypeError(
        "Failed to construct RangedSummarizedExperiment; tried variants: "
        + "; ".join(errors)
    )

def _normalize_genomic_unit(genomic_unit: str) -> str:
    """Return a normalized genomic unit string and validate it.

    Args:
      genomic_unit: Requested feature level.

    Returns:
      Lowercase genomic unit string.

    Raises:
      ValueError: If the genomic unit is not one of 'gene', 'exon',
        or 'junction'.
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
    """Resolve a single annotation extension for create_rse helpers.

    For gene/exon assays, this helper maps a human-readable annotation
    name (for example, "gencode_v26") or an explicit extension (for
    example, "G026") to the underlying recount3 extension. When both
    ``annotation`` and ``annotation_file_extension`` are provided, the
    explicit extension wins.

    Args:
      organism: Organism identifier ("human" or "mouse").
      genomic_unit: Normalized genomic unit ("gene", "exon", or
        "junction").
      annotation: Optional human-readable annotation label.
      annotation_file_extension: Optional explicit annotation extension.

    Returns:
      A single annotation extension string for gene/exon assays, or
      None when no annotation is required (for example, junction-level
      assays).
    """
    if genomic_unit not in {"gene", "exon"}:
        return None

    if annotation_file_extension:
        return str(annotation_file_extension).strip()

    if annotation:
        return annotation_ext(organism=organism, annotation=annotation)

    # Fall back to the search-layer defaults used by search_project_all().
    org = organism.strip().lower()
    if org == "human":
        return "G026"
    if org == "mouse":
        return "M023"

    raise ValueError(
        f"Unsupported organism {organism!r}; expected 'human' or 'mouse'."
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
) -> "RangedSummarizedExperiment | SummarizedExperiment":
    """High-level helper that mirrors recount3's create_rse() in R.

    This function hides the intermediate bundle construction step by:
      1) Discovering all resources for a project via :class:`R3Project`.
      2) Stacking expression matrices across projects/samples.
      3) Resolving genomic ranges (GTF for gene/exon; RR for junctions).
      4) Building a BiocPy RangedSummarizedExperiment (or SE fallback).

    Args:
      project: Study or project identifier (for example, "SRP009615").
      genomic_unit: One of {"gene", "exon", "junction"}.
      organism: Organism identifier ("human" or "mouse").
      data_source: Data source ("sra", "gtex", or "tcga").
      annotation: Optional human-readable annotation label for gene/exon
        assays, such as "gencode_v26" or "gencode_v29". Ignored for
        junction-level assays.
      annotation_file_extension: Optional explicit annotation extension
        (for example, "G026"). When provided, this takes precedence over
        ``annotation``. Ignored for junction-level assays.
      junction_type: Junction type; typically "ALL".
      junction_file_extensions: Iterable of junction artifact extensions
        to include (for example, ("MM",) or ("MM", "RR")). If None, the
        default ("MM",) is used.
      include_metadata: Whether to include the five project metadata
        tables in the underlying bundle (recommended).
      include_bigwig: Whether to include per-sample BigWig coverage
        resources in the bundle. These can be very large.
      join: Join policy across projects when stacking count matrices
        (passed to :func:`build_ranged_summarized_experiment`).
      autoload: If True, load resources transparently.
      allow_fallback_to_se: If True, construct a plain SE when genomic
        ranges cannot be derived for the requested combination.
      strict: If True, propagate validation errors from the search layer
        (for example, missing projects or incompatible combinations).

    Returns:
      A :class:`summarizedexperiment.RangedSummarizedExperiment` object,
      or a plain :class:`summarizedexperiment.SummarizedExperiment` when
      ``allow_fallback_to_se`` is True and ranges cannot be resolved.

    Raises:
      ImportError: If BiocPy packages are not installed.
      ValueError: If inputs are invalid or no counts are found.
      TypeError: If the underlying RangedSummarizedExperiment
        constructor rejects all compatibility variants.
    """
    _require_biocpy()

    unit = _normalize_genomic_unit(genomic_unit)
    ann_ext = _resolve_annotation_extension(
        organism=organism,
        genomic_unit=unit,
        annotation=annotation,
        annotation_file_extension=annotation_file_extension,
    )

    jexts: tuple[str, ...]
    if junction_file_extensions is None:
        jexts = ("MM",)
    else:
        jexts = tuple(str(e).strip() for e in junction_file_extensions)

    project_bundle = R3Project.discover(
        organism=organism,
        data_source=data_source,
        project=project,
        genomic_units=(unit,),
        annotations=("default" if ann_ext is None else (ann_ext,)),
        junction_exts=jexts,
        junction_type=junction_type,
        include_metadata=include_metadata,
        include_bigwig=include_bigwig,
        strict=strict,
        deduplicate=True,
    )

    return build_ranged_summarized_experiment(
        project_bundle,
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
) -> "RangedSummarizedExperiment | SummarizedExperiment":
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
        "RangedSummarizedExperiment"
        | "SummarizedExperiment"
        | pd.DataFrame
    ),
    *,
    column: str = "sra.sample_attributes",
    prefix: str = "sra_attribute.",
) -> (
    "RangedSummarizedExperiment"
    | "SummarizedExperiment"
    | pd.DataFrame
):
    """Expand encoded SRA sample attributes into separate columns.

    This function mirrors the recount3 R helper ``expand_sra_attributes()``.

    It understands the SRA encoding used by recount3, where a single column
    (typically ``'sra.sample_attributes'``) contains strings of the form::

        "age;;67.78|biomaterial_provider;;LIBD|disease;;Control|..."

    Each ``key;;value`` pair becomes a new column named
    ``{prefix}{key}`` (with spaces in ``key`` replaced by '_'), and the
    parsed values are stored per sample. The original string column is
    preserved.

    The function supports two calling styles:

    * Passing a :class:`pandas.DataFrame` of column metadata, in which case a
      new DataFrame is returned with extra columns.
    * Passing a BiocPy
      :class:`summarizedexperiment.SummarizedExperiment` or
      :class:`summarizedexperiment.RangedSummarizedExperiment`, in which case
      a new object of the same class is returned with updated
      ``column_data``.

    Args:
      obj: Column metadata DataFrame or a BiocPy SE/RSE object.
      column: Name of the column that holds the raw SRA attribute strings.
      prefix: Prefix to prepend to the generated attribute column names.

    Returns:
      A new object of the same type as ``obj`` with additional columns
      corresponding to parsed SRA attributes. If the requested column is
      missing, ``obj`` is returned unchanged.

    Raises:
      ImportError: If a SummarizedExperiment/RangedSummarizedExperiment is
        supplied but BiocPy packages are not installed.
      AttributeError: If the BiocPy object does not expose ``column_data`` or
        ``set_column_data`` in the expected API.
      TypeError: If ``obj`` is neither a pandas DataFrame nor a supported
        BiocPy experiment object.
    """
    # DataFrame mode: no BiocPy dependency.
    if isinstance(obj, pd.DataFrame):
        return _expand_sra_attributes_df(
            obj,
            column=column,
            prefix=prefix,
        )

    # SummarizedExperiment / RangedSummarizedExperiment mode.
    _require_biocpy()

    # After _require_biocpy() these are guaranteed to be real classes.
    if not isinstance(obj, (SummarizedExperiment, RangedSummarizedExperiment)):
        raise TypeError(
            "expand_sra_attributes() expects a pandas DataFrame or a BiocPy "
            f"SummarizedExperiment/RangedSummarizedExperiment; got {type(obj)!r}."
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


# --------------------------- constructors -----------------------------------


def _construct_se_compat(
    *,
    counts_df: pd.DataFrame,
    row_df: pd.DataFrame,
    col_df: pd.DataFrame,
    assay_name: str,
) -> "SummarizedExperiment":
    """Construct an SE using several compatible constructor variants.

    This enforces strict shape/name alignment and then tries multiple shapes
    supported by different SummarizedExperiment versions:

      1) assays={name: ndarray}, row_data=DataFrame, column_data=DataFrame
      2) assays=[ndarray], assay_names=[name], row_data=DataFrame, column_data=DataFrame
      3) Same as (1) but with `col_data=` instead of `column_data=`

    Args:
      counts_df: Feature x sample counts matrix (DataFrame).
      row_df: Row data (will be aligned and indexed to counts_df.index).
      col_df: Column data (will be aligned and indexed to counts_df.columns).
      assay_name: Assay slot name (e.g., "counts").

    Returns:
      A :class:`summarizedexperiment.SummarizedExperiment`.

    Raises:
      ValueError: If shapes do not match.
      TypeError: If no constructor variant is accepted.
    """
    _require_biocpy()

    # Align and validate shapes.
    n_features, n_samples = counts_df.shape
    if n_features == 0 or n_samples == 0:
        raise ValueError("Empty assay matrix; cannot build SE.")
    rd = row_df.copy()
    cd = col_df.copy()
    rd.index = counts_df.index
    cd.index = counts_df.columns
    rd = _ensure_unique_columns(rd, empty_prefix="row")
    cd = _ensure_unique_columns(cd, empty_prefix="col")
    if len(rd) != n_features:
        raise ValueError(
            f"row_data length {len(rd)} != assay features {n_features}"
        )
    if len(cd) != n_samples:
        raise ValueError(
            f"column_data length {len(cd)} != assay samples {n_samples}"
        )

    counts_numeric = counts_df.apply(pd.to_numeric, errors="coerce").fillna(0)
    counts_np = counts_numeric.to_numpy(copy=False)

    errors: list[str] = []

    # Variant 1: dict[str, ndarray] + column_data
    try:
        return SummarizedExperiment(
            assays={assay_name: counts_np},
            row_data=rd,
            column_data=cd,
        )
    except Exception as exc:  # pylint: disable=broad-except
        errors.append(f"dict[str, ndarray] + column_data → {exc!r}")

    # Variant 2: list[ndarray] + assay_names + column_data
    try:
        return SummarizedExperiment(
            assays=[counts_np],
            assay_names=[assay_name],
            row_data=rd,
            column_data=cd,
        )
    except Exception as exc:  # pylint: disable=broad-except
        errors.append(f"list[ndarray]+assay_names + column_data → {exc!r}")

    # Variant 3: swap keyword to col_data
    try:
        return SummarizedExperiment(
            assays={assay_name: counts_np},
            row_data=rd,
            col_data=cd,
        )
    except Exception as exc:  # pylint: disable=broad-except
        errors.append(f"dict[str, ndarray] + col_data → {exc!r}")

    raise TypeError(
        "Failed to construct SummarizedExperiment; tried variants: "
        + "; ".join(errors)
    )
