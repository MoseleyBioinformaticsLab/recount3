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
#
"""Resource bundles, project discovery, and concatenation helpers.

This module defines :class:`R3ResourceBundle`, a general-purpose
container for groups of :class:`~recount3.resource.R3Resource` objects.

Bundles support lazy loading, filtering by description fields,
project-aware discovery, and high-level helpers for combining recount3
resources into BiocPy objects such as
:class:`~summarizedexperiment.SummarizedExperiment` and
:class:`~summarizedexperiment.RangedSummarizedExperiment`.

When discovery covers exactly one ``(organism, data_source, project)``
triple, the bundle's ``organism``, ``data_source``, and ``project``
attributes are set accordingly. For multi-project bundles these attributes
are ``None`` to avoid misrepresenting the identity.

Filtering with FieldSpec
  :meth:`~R3ResourceBundle.filter` accepts a
  :data:`~recount3.types.FieldSpec` for each description field. Three
  forms are accepted:

  * A string: exact match (e.g. ``genomic_unit="gene"``).
  * An iterable of strings: keep if the field is any of the given values.
  * A callable: called with the field value; truthy return keeps the resource.
  * ``None`` (default): no filtering on that field.

Typical usage example::

  from recount3 import R3ResourceBundle

  bundle = R3ResourceBundle.discover(
      organism="human",
      data_source="sra",
      project="SRP009615",
  )

  # Filter to gene-level count resources and stack into a DataFrame:
  counts = bundle.filter(
      resource_type="count_files_gene_or_exon",
      genomic_unit="gene",
  ).stack_count_matrices()

  # Filter with a callable predicate:
  meta_only = bundle.filter(
      resource_type=lambda t: t == "metadata_files"
  )

Note:
    The :meth:`~R3ResourceBundle.to_summarized_experiment` and
    :meth:`~R3ResourceBundle.to_ranged_summarized_experiment` methods
    require the BiocPy package ``summarizedexperiment``, which might be
    difficult to install on Windows. Install with::

        pip install summarizedexperiment
"""

from __future__ import annotations

import dataclasses
import functools
import gzip
import io
import logging
import re
from collections import Counter
from collections.abc import Callable, Iterable, Iterator, Sequence
from typing import Any, Optional, Mapping, TYPE_CHECKING

import pandas as pd
import numpy as np

from recount3 import _bigwig
from recount3 import _utils
from recount3 import errors
from recount3 import resource
from recount3 import search
from recount3 import types as r3_types

if TYPE_CHECKING:  # pragma: no cover
    import summarizedexperiment  # type: ignore[import-not-found]
    import genomicranges  # type: ignore[import-not-found]

_EXON_ID_ATTR_RE = re.compile(r'\bexon_id\s+"([^"]+)"')
_RECOUNT_EXON_ID_ATTR_RE = re.compile(r'\brecount_exon_id\s+"([^"]+)"')


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


def _default_assay_name(genomic_unit: str, assay_name: str) -> str:
    """Return the appropriate assay name for the given genomic unit.

    The R reference package uses ``"raw_counts"`` for gene and exon assays
    (which represent base-pair coverage sums) and ``"counts"`` for junction
    assays (which are actual read counts). This helper applies the same
    convention when the caller uses the default assay name.

    Args:
      genomic_unit: One of ``"gene"``, ``"exon"``, or ``"junction"``.
      assay_name: The caller-provided assay name.

    Returns:
      The resolved assay name.
    """
    if assay_name == "raw_counts" and genomic_unit == "junction":
        return "counts"
    return assay_name


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
    origin = getattr(desc, "table_name", None) or getattr(
        desc, "resource_type", None
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
    except Exception:  # pylint: disable=broad-exception-caught
        pass

    try:
        path = res._cached_path()  # pylint: disable=protected-access
    except Exception as exc:
        raise ValueError(
            f"Cannot resolve cached RR path for: {res.url}"
        ) from exc

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:  # type: ignore[arg-type]
        text = io.TextIOWrapper(fh, encoding="utf-8")
        df = pd.read_csv(text, sep="\t", comment="#")

    return df


_GTF_ATTR_PAIR_RE = re.compile(
    r'\s*(?P<key>[^\s;]+)\s+"?(?P<value>[^";]+)"?\s*(?:;|$)'
)


def _parse_gtf_attributes(attrs: pd.Series) -> pd.DataFrame:
    """Parse a GTF attributes column into a wide (column-per-key) DataFrame.

    Args:
      attrs: Series containing raw GTF attribute strings.

    Returns:
      DataFrame indexed like attrs, with one column per attribute key.
    """
    if attrs.empty:
        return pd.DataFrame(index=attrs.index)

    extracted = (
        attrs.fillna("").astype("string").str.extractall(_GTF_ATTR_PAIR_RE)
    )
    if extracted.empty:
        return pd.DataFrame(index=attrs.index)

    extracted = extracted.reset_index(level=1, drop=True)
    extracted.index.name = "_row"

    wide = (
        extracted.reset_index()
        .pivot_table(
            index="_row",
            columns="key",
            values="value",
            aggfunc="first",  # pyright: ignore[reportArgumentType]
        )
        .reindex(attrs.index)
    )
    wide.columns = [str(c) for c in wide.columns]
    return wide


def _coerce_gtf_phase_column(phase_column: pd.Series) -> pd.Series:
    """Coerce a GTF frame/phase column to a nullable integer Series.

    GTF "frame" (a.k.a. phase) is typically one of {"0", "1", "2"} or "."
    for features where it is not applicable. This returns a pandas nullable
    integer (Int64) Series where "." and missing values become <NA>.

    Args:
      phase_column: A Series containing the 8th GTF column ("frame"/phase).

    Returns:
      A pandas Series of dtype Int64 with values in {0, 1, 2} or <NA>.
    """
    # Normalize to string to handle mixed types (object, int, float, NA).
    frame_str = phase_column.astype("string").str.strip()

    # GTF uses "." to indicate missing.
    frame_str = frame_str.replace(
        {".": pd.NA, "": pd.NA}  # pyright: ignore[reportArgumentType]
    )  # pyright: ignore[reportArgumentType]

    phase = pd.to_numeric(frame_str, errors="coerce").astype("Int64")

    # Validate: phase should be 0/1/2 when present.
    invalid = phase.notna() & ~phase.isin([0, 1, 2])
    if invalid.any():
        bad_values = sorted(
            {str(v) for v in frame_str[invalid].dropna().unique()}
        )
        logging.warning(
            "GTF phase/frame column contains unexpected values; "
            "coercing them to <NA>. Values=%s",
            bad_values[:10],
        )
        phase = phase.mask(invalid, pd.NA)

    return phase


def _coerce_gtf_bp_length(
    score: pd.Series,
    *,
    starts: pd.Series,
    ends: pd.Series,
) -> pd.Series:
    """Derive a bp_length column for GTF features.

    Args:
      score: GTF "score" column (6th column). Often "." in standard GTFs.
      starts: Numeric start positions (1-based, inclusive).
      ends: Numeric end positions (1-based, inclusive).

    Returns:
      A pandas Series of dtype Int64 giving bp_length (width).
    """
    starts_num = pd.to_numeric(starts, errors="coerce").astype("Int64")
    ends_num = pd.to_numeric(ends, errors="coerce").astype("Int64")
    width = (ends_num - starts_num + 1).astype("Int64")

    score_str = (
        score.astype("string").str.strip().replace({".": pd.NA, "": pd.NA})  # pyright: ignore[reportArgumentType]
    )  # pyright: ignore[reportArgumentType]
    score_num = pd.to_numeric(score_str, errors="coerce")

    comparable = score_num.notna() & width.notna()
    if not comparable.any():
        return width

    score_rounded = score_num.round()
    is_integer_like = (score_num - score_rounded).abs() <= 1e-6

    score_int = score_rounded.astype("Int64")
    matches_width = comparable & is_integer_like & (score_int == width)

    match_rate = float(matches_width.sum()) / float(comparable.sum())

    # Tolerate small noise.
    if match_rate >= 0.95:
        return score_int.where(score_int.notna(), width)

    return width


_ENSEMBL_VERSION_SUFFIX_RE = r"\.\d+$"


def _strip_ensembl_version(values: pd.Series) -> pd.Series:
    """Strip trailing Ensembl version suffix (e.g. ENSG... .12 -> ENSG...).

    Args:
      values: Series of identifiers.

    Returns:
      Series with trailing '.<digits>' removed when present.
    """
    return (
        values.astype("string")
        .str.strip()
        .str.replace(_ENSEMBL_VERSION_SUFFIX_RE, "", regex=True)
    )


def _align_ranges_to_features(
    ranges: pd.DataFrame,
    *,
    feature_ids: Sequence[str],
) -> pd.DataFrame:
    """Align annotation ranges to a list of feature IDs.

    Tries exact feature_id matching first. If any features are missing,
    tries a secondary match after stripping Ensembl version suffixes
    ('.<digits>') on both sides.

    Args:
      ranges: DataFrame containing a 'feature_id' column and coordinate columns.
      feature_ids: Feature IDs from the counts matrix (in desired order).

    Returns:
      DataFrame reindexed to feature_ids order, containing the same columns as
      `ranges` except for 'feature_id' (which is used as the index).

    Raises:
      ValueError: If required coordinate columns are missing from `ranges`, or if
        version-stripped matching is ambiguous (conflicting coordinates).
    """
    required = {"feature_id", "seqnames", "starts", "ends", "strand"}
    missing_cols = required - set(ranges.columns)
    if missing_cols:
        raise ValueError(
            "ranges is missing required columns: " f"{sorted(missing_cols)}"
        )

    feature_index = pd.Index([str(x) for x in feature_ids])

    exact = ranges.set_index("feature_id").reindex(feature_index)

    coord_cols = ["seqnames", "starts", "ends", "strand"]
    missing_any = exact[coord_cols].isna().any(axis=1)
    if not missing_any.any():
        return exact

    # Build a version-stripped index for the annotation ranges.
    ranges_uv = ranges.copy()
    ranges_uv["_feature_id_unversioned"] = _strip_ensembl_version(
        ranges_uv["feature_id"]
    )

    dup_mask = ranges_uv["_feature_id_unversioned"].duplicated(keep=False)
    if dup_mask.any():
        dup = ranges_uv.loc[dup_mask, ["_feature_id_unversioned"] + coord_cols]
        nunique = dup.groupby("_feature_id_unversioned", sort=False)[
            coord_cols
        ].nunique(dropna=False)
        conflicting = nunique.max(axis=1) > 1
        if conflicting.any():
            bad = list(nunique[conflicting].index[:10])
            raise ValueError(
                "After stripping Ensembl versions, multiple annotation rows map to "
                "the same ID but with conflicting coordinates. Example IDs: "
                f"{bad!r}"
            )
        ranges_uv = ranges_uv.drop_duplicates(
            "_feature_id_unversioned", keep="first"
        )

    ranges_uv = ranges_uv.set_index("_feature_id_unversioned")

    feat_uv = _strip_ensembl_version(pd.Series(feature_index, dtype="string"))
    fallback = ranges_uv.reindex(pd.Index(feat_uv.astype(str)))
    fallback.index = feature_index

    # Fill missing exact matches with version-stripped matches.
    filled = exact.combine_first(fallback)

    if filled[coord_cols].isna().any(axis=1).any():
        # Left to caller to decide whether to error
        return filled

    logging.info(
        "Aligned annotation ranges using Ensembl version-stripped fallback for "
        "%d features.",
        int(missing_any.sum()),
    )
    return filled


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
        raise ValueError(
            f"Cannot resolve cached GTF path for: {res.url}"
        ) from exc

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

    attrs = df["attributes"].astype("string")
    attrs_wide = _parse_gtf_attributes(attrs)

    if feature_kind == "gene" and "gene_id" in attrs_wide.columns:
        feature_ids = attrs_wide["gene_id"].astype("string")
    elif feature_kind != "gene" and "recount_exon_id" in attrs_wide.columns:
        feature_ids = attrs_wide["recount_exon_id"].astype("string")
    elif feature_kind != "gene" and "exon_id" in attrs_wide.columns:
        feature_ids = attrs_wide["exon_id"].astype("string")
    else:
        feature_ids = attrs.str.extract(
            _RECOUNT_EXON_ID_ATTR_RE,
            expand=False,
        )
        missing = feature_ids.isna()
        if missing.any():  # pragma: no branch
            exon_ids = attrs[missing].str.extract(
                _EXON_ID_ATTR_RE,
                expand=False,
            )
            feature_ids = feature_ids.fillna(exon_ids)

    missing = feature_ids.isna()
    if missing.any():
        # Last-resort fallback (does not guarantee compatibility with counts).
        coords = (
            df.loc[missing, "seqname"].astype(str)
            + "|"
            + df.loc[missing, "start"].astype(str)
            + "|"
            + df.loc[missing, "end"].astype(str)
            + "|"
            + df.loc[missing, "strand"].astype(str)
        )
        feature_ids.loc[missing] = coords.values

    out = pd.DataFrame(
        {
            "seqnames": df["seqname"].astype("string"),
            "starts": df["start"].astype(int),
            "ends": df["end"].astype(int),
            "strand": df["strand"].astype("string"),
            "feature_id": feature_ids.astype("string"),
            "source": df["source"].astype("string"),
            "type": df["feature"].astype("string"),
            "bp_length": _coerce_gtf_bp_length(
                df["score"], starts=df["start"], ends=df["end"]
            ),
            "phase": _coerce_gtf_phase_column(df["frame"]),
        },
        index=df.index,
    )

    if not attrs_wide.empty:
        attrs_extra = attrs_wide.drop(columns=out.columns, errors="ignore")
        out = out.join(attrs_extra)

    if "level" in out.columns:
        out["level"] = pd.to_numeric(out["level"], errors="coerce").astype(
            "Int64"
        )

    if out["feature_id"].duplicated().any():
        dup_mask = out["feature_id"].duplicated(keep=False)
        dup = out.loc[
            dup_mask,
            ["feature_id", "seqnames", "starts", "ends", "strand"],
        ]
        nunique = dup.groupby("feature_id")[
            ["seqnames", "starts", "ends", "strand"]
        ].nunique()
        conflicting = nunique.max(axis=1) > 1
        if conflicting.any():
            example_ids = list(nunique[conflicting].index[:5])
            raise ValueError(
                "GTF contains duplicate feature_id values with conflicting "
                f"genomic ranges (example feature_id values: {example_ids})."
            )

        dropped = int(dup["feature_id"].nunique())
        logging.info(
            "Dropping duplicate GTF rows for %d feature_id values with "
            "identical ranges.",
            dropped,
        )
        out = out.drop_duplicates(subset=["feature_id"], keep="first")

    return out


def _peek_gtf_feature_counts(
    res: resource.R3Resource, *, max_lines: int = 50000
) -> Counter[str]:
    """Quickly scan the GTF and count feature types without loading the whole file."""
    try:
        path = res._cached_path()  # pylint: disable=protected-access
    except Exception:  # pylint: disable=broad-exception-caught
        res.download(path=None, cache_mode="enable")
        path = res._cached_path()  # pylint: disable=protected-access

    opener = gzip.open if str(path).endswith(".gz") else open
    counts: Counter[str] = Counter()
    seen = 0
    with opener(path, "rb") as fh:  # type: ignore[arg-type]
        for line in io.TextIOWrapper(fh, encoding="utf-8"):
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                counts[parts[2]] += 1
                seen += 1
            if seen >= max_lines:
                break
    return counts


def _select_gtf_resource_for_unit(
    bundle: R3ResourceBundle,
    *,
    genomic_unit: str,
    annotation_extension: Optional[str],
) -> Optional[resource.R3Resource]:
    """Pick the most appropriate annotation resource for gene/exon ranges."""
    feature_kind = "gene" if genomic_unit == "gene" else "exon"

    ann = bundle.filter(resource_type="annotations")
    if annotation_extension:
        ann = ann.filter(annotation_extension=annotation_extension)

    candidates = list(ann.resources)
    if not candidates:
        return None

    def score(res: resource.R3Resource) -> int:
        desc = res.description
        path = (desc.url_path() if hasattr(desc, "url_path") else "") or ""
        url = res.url or ""
        text = f"{path} {url}".lower()

        s = 0
        gu = getattr(desc, "genomic_unit", None)
        if isinstance(gu, str) and gu.lower() == genomic_unit:
            s += 200

        if genomic_unit == "gene" and ("gene" in text or "genes" in text):
            s += 120
        if genomic_unit == "exon" and ("exon" in text or "exons" in text):
            s += 120

        if ".gtf" in text:
            s += 20
        return s

    candidates = sorted(candidates, key=score, reverse=True)

    for res in candidates:
        try:
            feats = _peek_gtf_feature_counts(res, max_lines=50000)
            if feats.get(feature_kind, 0) > 0:
                logging.info(
                    "Selected annotation resource for %s: %s (peek features=%s)",
                    genomic_unit,
                    res.description.url_path(),
                    dict(feats.most_common(5)),
                )
                return res
        except Exception as exc:  # pylint: disable=broad-exception-caught
            logging.warning(
                "Failed to peek GTF features for %s: %r", res.url, exc
            )

    return candidates[0]


def _to_genomic_ranges(ranges_df: pd.DataFrame) -> genomicranges.GenomicRanges:
    """Construct a GenomicRanges object from a pandas DataFrame.

    Args:
      ranges_df: DataFrame with at least ``seqnames``, ``starts``,
        ``ends``, and ``strand`` columns.

    Returns:
      A :class:`genomicranges.GenomicRanges` instance.

    Raises:
      ImportError: If :mod:`genomicranges` cannot be imported.
    """
    genomic_ranges_cls = _utils.get_genomicranges_class()
    return genomic_ranges_cls.from_pandas(ranges_df)


def _construct_summarized_experiment(
    *,
    counts_df: pd.DataFrame,
    row_df: pd.DataFrame,
    col_df: pd.DataFrame,
    assay_name: str,
    metadata: Optional[Mapping[str, Any]] = None,
) -> summarizedexperiment.SummarizedExperiment:
    """Construct a SummarizedExperiment using the latest BiocPy API.

    This function targets the modern `summarizedexperiment.SummarizedExperiment`
    constructor that accepts:

      - assays: dict[str, 2D matrix]
      - row_data: feature metadata (coerced to BiocFrame by BiocPy)
      - column_data: sample metadata (coerced to BiocFrame by BiocPy)
      - row_names / column_names: optional but provided explicitly here
      - metadata: optional experiment-level metadata

    It enforces:
      - non-empty assay matrix
      - strict shape agreement between assay and row/column data
      - numeric assay values (raises on non-numeric, fills missing with 0)

    Args:
      counts_df: Feature-by-sample counts matrix.
      row_df: Feature metadata DataFrame (will be reindexed to counts_df.index).
      col_df: Sample metadata DataFrame (will be reindexed to counts_df.columns).
      assay_name: Name of the assay (e.g. "raw_counts").
      metadata: Optional experiment-level metadata.

    Returns:
      A SummarizedExperiment instance.

    Raises:
      ValueError: If shapes are inconsistent or assay contains invalid values.
    """
    summarized_experiment_cls = _utils.get_summarizedexperiment_class()

    if counts_df.ndim != 2:
        raise ValueError("counts_df must be a 2D DataFrame.")

    n_features, n_samples = counts_df.shape
    if n_features == 0 or n_samples == 0:
        raise ValueError(
            "Empty assay matrix; cannot build SummarizedExperiment."
        )

    row_names = [str(x) for x in counts_df.index]
    col_names = [str(x) for x in counts_df.columns]

    row_data = row_df.copy()
    col_data = col_df.copy()
    row_data.index = row_names
    col_data.index = col_names

    row_data = _ensure_unique_columns(row_data, empty_prefix="row")
    col_data = _ensure_unique_columns(col_data, empty_prefix="col")

    if len(row_data) != n_features:
        raise ValueError(
            f"row_data length {len(row_data)} != assay features {n_features}."
        )
    if len(col_data) != n_samples:
        raise ValueError(
            f"column_data length {len(col_data)} != assay samples {n_samples}."
        )

    raw = counts_df.copy()
    numeric = raw.apply(pd.to_numeric, errors="coerce")

    invalid = raw.notna() & numeric.isna()
    if invalid.any().any():
        bad_positions = list(zip(*np.where(invalid.to_numpy())))
        examples = []
        for r, c in bad_positions[:5]:
            examples.append(
                {
                    "row": row_names[r],
                    "col": col_names[c],
                    "value": raw.iat[r, c],
                }
            )
        raise ValueError(
            "counts_df contains non-numeric values that cannot be coerced. "
            f"Examples: {examples!r}"
        )

    counts_np = numeric.fillna(0).to_numpy(dtype=float, copy=False)

    biocframe_cls = _utils.get_biocframe_class()
    return summarized_experiment_cls(
        assays={assay_name: counts_np},
        row_data=biocframe_cls.from_pandas(row_data),
        column_data=biocframe_cls.from_pandas(col_data),
        row_names=row_names,
        column_names=col_names,
        metadata=dict(metadata) if metadata is not None else None,
    )


def _construct_ranged_summarized_experiment(
    *,
    counts_df: pd.DataFrame,
    row_df: pd.DataFrame,
    col_df: pd.DataFrame,
    ranges_df: pd.DataFrame,
    assay_name: str,
    metadata: Optional[Mapping[str, Any]] = None,
) -> summarizedexperiment.RangedSummarizedExperiment:
    """Construct a RangedSummarizedExperiment using the latest BiocPy API.

    This targets the modern `summarizedexperiment.RangedSummarizedExperiment`
    constructor with kwargs:

      - assays: dict[str, 2D matrix]
      - row_ranges: GenomicRanges
      - row_data / column_data: pandas DataFrames (coerced to BiocFrame by BiocPy)
      - row_names / column_names: provided explicitly for stability
      - metadata: optional experiment metadata

    Enforces:
      - non-empty 2D counts matrix
      - strict shape agreement between assay, metadata, and ranges
      - numeric assay values (raises on non-numeric non-missing values)

    Args:
      counts_df: Feature-by-sample counts matrix.
      row_df: Feature metadata (reindexed to counts_df.index).
      col_df: Sample metadata (reindexed to counts_df.columns).
      ranges_df: DataFrame with at least ["seqnames", "starts", "ends", "strand"].
        Must have the same number/order of rows as counts_df.
      assay_name: Name of the assay (e.g. "raw_counts").
      metadata: Optional experiment-level metadata.

    Returns:
      A RangedSummarizedExperiment instance.

    Raises:
      ValueError: If shapes/columns are inconsistent or counts are not numeric.
    """
    ranged_summarized_experiment_cls = (
        _utils.get_ranged_summarizedexperiment_class()
    )

    if counts_df.ndim != 2:
        raise ValueError("counts_df must be a 2D DataFrame.")

    n_features, n_samples = counts_df.shape
    if n_features == 0 or n_samples == 0:
        raise ValueError(
            "Empty assay matrix; cannot build RangedSummarizedExperiment."
        )

    required_range_cols = {"seqnames", "starts", "ends", "strand"}
    missing = required_range_cols - set(ranges_df.columns)
    if missing:
        raise ValueError(
            f"ranges_df is missing required columns: {sorted(missing)!r}."
        )

    if len(ranges_df) != n_features:
        raise ValueError(
            f"ranges_df length {len(ranges_df)} != assay features {n_features}."
        )

    row_names = [str(x) for x in counts_df.index]
    col_names = [str(x) for x in counts_df.columns]

    row_data = row_df.copy()
    col_data = col_df.copy()
    row_data.index = row_names
    col_data.index = col_names

    row_data = _ensure_unique_columns(row_data, empty_prefix="row")
    col_data = _ensure_unique_columns(col_data, empty_prefix="col")

    if len(row_data) != n_features:
        raise ValueError(
            f"row_data length {len(row_data)} != assay features {n_features}."
        )
    if len(col_data) != n_samples:
        raise ValueError(
            f"column_data length {len(col_data)} != assay samples {n_samples}."
        )

    ranges_aligned = ranges_df.copy()
    ranges_aligned.index = row_names
    ranges_aligned["starts"] = pd.to_numeric(
        ranges_aligned["starts"], errors="coerce"
    )
    ranges_aligned["ends"] = pd.to_numeric(
        ranges_aligned["ends"], errors="coerce"
    )
    if (
        ranges_aligned[["seqnames", "starts", "ends", "strand"]]
        .isna()
        .any()
        .any()
    ):
        raise ValueError(
            "ranges_df contains missing values in coordinate columns."
        )

    raw = counts_df.copy()
    numeric = raw.apply(pd.to_numeric, errors="coerce")
    invalid = raw.notna() & numeric.isna()
    if invalid.any().any():
        bad_positions = list(zip(*np.where(invalid.to_numpy())))
        examples = []
        for r, c in bad_positions[:5]:
            examples.append(
                {
                    "row": row_names[r],
                    "col": col_names[c],
                    "value": raw.iat[r, c],
                }
            )
        raise ValueError(
            "counts_df contains non-numeric values that cannot be coerced. "
            f"Examples: {examples!r}"
        )

    counts_np = numeric.fillna(0).to_numpy(dtype=float, copy=False)

    gr = _to_genomic_ranges(ranges_aligned)

    biocframe_cls = _utils.get_biocframe_class()
    return ranged_summarized_experiment_cls(
        assays={assay_name: counts_np},
        row_ranges=gr,
        row_data=biocframe_cls.from_pandas(row_data),
        column_data=biocframe_cls.from_pandas(col_data),
        row_names=row_names,
        column_names=col_names,
        metadata=dict(metadata) if metadata is not None else None,
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
            getattr(res.description, "junction_extension", None) or ""
        )
        family = "junctions"
        feature_key = f"{family}:{junction_type}:{junction_ext}"
        return family, feature_key

    raise ValueError(
        "Resource is not a recognized count-file type for stacking: "
        f"{rtype!r}"
    )


def _make_unique_names(
    names: Sequence[str],
    *,
    suffix: str = "__dup",
) -> list[str]:
    """Return a stable, order-preserving unique-ified version of `names`.

    Example:
      ["a", "b", "a"] -> ["a", "b", "a__dup2"]
    """
    seen: Counter[str] = Counter()
    out: list[str] = []
    for name in names:
        seen[name] += 1
        if seen[name] == 1:
            out.append(name)
        else:
            out.append(f"{name}{suffix}{seen[name]}")
    return out


def _dedupe_ranges_on_feature_id(ranges: pd.DataFrame) -> pd.DataFrame:
    """Ensure `ranges.feature_id` is unique, erroring if duplicates disagree."""
    if "feature_id" not in ranges.columns:
        raise ValueError("ranges is missing required column 'feature_id'.")

    dup_mask = ranges["feature_id"].duplicated(keep=False)
    if not dup_mask.any():
        return ranges

    key_cols = ["seqnames", "starts", "ends", "strand"]
    have_cols = [c for c in key_cols if c in ranges.columns]
    if len(have_cols) != len(key_cols):
        raise ValueError(
            "ranges is missing one or more required columns: "
            f"{set(key_cols) - set(have_cols)}"
        )

    dup = ranges.loc[dup_mask, ["feature_id"] + key_cols]
    nunique = dup.groupby("feature_id", sort=False)[key_cols].nunique(
        dropna=False
    )
    inconsistent = nunique.max(axis=1) > 1
    if inconsistent.any():
        bad_ids = list(inconsistent[inconsistent].index[:10])
        raise ValueError(
            "Duplicate feature_id values map to different coordinates; "
            f"example feature_ids={bad_ids!r}"
        )

    logging.warning(
        "Annotation ranges contain duplicate feature_id values; keeping the "
        "first occurrence for alignment."
    )
    return ranges.drop_duplicates(subset=["feature_id"], keep="first")


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

        Examples:
            Discover all default resources for a single project::

                bundle = R3ResourceBundle.discover(
                    organism="human",
                    data_source="sra",
                    project="SRP009615",
                )

            Discover gene counts only across two projects::

                bundle = R3ResourceBundle.discover(
                    organism="human",
                    data_source="sra",
                    project=["SRP009615", "SRP012682"],
                    genomic_units=("gene",),
                )

            Include BigWig coverage files alongside counts::

                bundle = R3ResourceBundle.discover(
                    organism="human",
                    data_source="sra",
                    project="SRP009615",
                    include_bigwig=True,
                )
        """
        organism_values = search._normalize_to_tuple(organism)
        data_source_values = search._normalize_to_tuple(data_source)
        project_values = search._normalize_to_tuple(project)

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
                        junction_extension=junction_exts,
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
        annotation_extension: Optional[r3_types.FieldSpec] = None,
        junction_extension: Optional[r3_types.FieldSpec] = None,
        predicate: Optional[Callable[[resource.R3Resource], bool]] = None,
        invert: bool = False,
    ) -> R3ResourceBundle:
        """Return a new bundle containing resources that match criteria.

        Each keyword argument corresponds to an attribute on
        :class:`recount3._descriptions.R3ResourceDescription`. Values are
        interpreted using :func:`recount3.search.match_spec`, allowing
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
          annotation_extension: Annotation code filter.
          junction_extension: Junction extension filter.
          predicate: Optional callback that receives each resource and
            returns :data:`True` if it should be kept.
          invert: If :data:`True`, invert the final match decision.

        Returns:
          A new :class:`R3ResourceBundle` containing only the resources
          that match all supplied filters and the optional predicate.

        Examples:
            Keep only gene-level resources::

                gene_bundle = bundle.filter(genomic_unit="gene")

            Keep gene or exon resources (iterable form)::

                ge_bundle = bundle.filter(genomic_unit=["gene", "exon"])

            Keep resources whose type contains "count" (callable form)::

                counts = bundle.filter(
                    resource_type=lambda t: "count" in (t or "")
                )

            Invert a filter to exclude metadata tables::

                no_meta = bundle.filter(
                    resource_type="metadata_files", invert=True
                )
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
            "annotation_extension": annotation_extension,
            "junction_extension": junction_extension,
        }
        field_specs = {
            key: value
            for key, value in field_specs.items()
            if value is not None
        }

        selected: list[resource.R3Resource] = []
        for res in self.resources:
            desc = res.description
            fields_ok = all(
                search.match_spec(
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
        join_policy: str = "inner",
        axis: int = 1,
        verify_integrity: bool = False,
        autoload: bool = True,
        compat: r3_types.CompatibilityMode = "family",
    ) -> pd.DataFrame:
        """Concatenate count matrices (gene/exon or junction) as DataFrames.

        Args:
          join_policy: Join policy passed to :func:`pandas.concat`.
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

        Examples:
            Stack gene counts across all projects in the bundle::

                df = bundle.filter(
                    resource_type="count_files_gene_or_exon",
                    genomic_unit="gene",
                ).stack_count_matrices()

            Require identical feature sets (faster; fails if annotations differ)::

                df = bundle.filter(
                    resource_type="count_files_gene_or_exon"
                ).stack_count_matrices(compat="feature")
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
                    "'junction_extension' for junctions."
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
            join=join_policy,  # type: ignore[arg-type]
            verify_integrity=verify_integrity,
        )

    def _stack_counts_for(
        self,
        *,
        genomic_unit: str,
        join_policy: str = "inner",
        autoload: bool = True,
    ) -> pd.DataFrame:
        """Return a wide counts DataFrame for the requested feature family.

        This is a bundle-scoped helper used by the SummarizedExperiment
        builders. It enforces appropriate compatibility within the
        gene/exon family or the junctions family.

        Args:
          genomic_unit: One of ``"gene"``, ``"exon"``, or ``"junction"``.
          join_policy: Join mode for concatenation (see :func:`pandas.concat`).
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
            sel = self.filter(resource_type="count_files_gene_or_exon").filter(
                genomic_unit=genomic_unit
            )

            if autoload:
                for res in sel.resources:
                    try:
                        res.load()
                    except Exception as exc:  # pylint: disable=broad-except
                        url = (
                            res.url or f"<no-url:{res.description.url_path()}>"
                        )
                        load_errors.append((url, exc))

            try:
                return sel.stack_count_matrices(
                    join_policy=join_policy,
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

        sel = self.filter(resource_type="count_files_junctions").filter(
            junction_extension="MM"
        )
        if autoload:
            for res in sel.resources:
                try:
                    res.load()
                except Exception as exc:  # pylint: disable=broad-except
                    url = res.url or f"<no-url:{res.description.url_path()}>"
                    load_errors.append((url, exc))

        try:
            return sel.stack_count_matrices(
                join_policy=join_policy,
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

    def _add_bigwig_urls(
        self,
        col_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Add a ``BigWigURL`` column to sample metadata.

        For each sample (identified by ``external_id`` in ``col_df``), this
        method constructs the duffel URL for its BigWig coverage file by
        inspecting the bundle's count resources to infer ``organism``,
        ``data_source``, and ``project``.

        This matches the R reference behavior in ``create_rse_manual()``
        where ``metadata$BigWigURL`` is populated via ``locate_url()``.

        Args:
          col_df: Sample metadata DataFrame with an ``external_id`` column.

        Returns:
          A copy of ``col_df`` with a ``BigWigURL`` column appended. If
          the necessary resource attributes cannot be inferred, the column
          contains ``pd.NA``.
        """
        from recount3._descriptions import R3BigWig
        from recount3.config import default_config

        out = col_df.copy()

        if "external_id" not in out.columns:
            out["BigWigURL"] = pd.NA
            return out

        count_types = {"count_files_gene_or_exon", "count_files_junctions"}
        organism = data_source = project = None
        for res in self.resources:
            desc = res.description
            rt = getattr(desc, "resource_type", None)
            if rt in count_types:
                organism = getattr(desc, "organism", None)
                data_source = getattr(desc, "data_source", None)
                project = getattr(desc, "project", None)
                if organism and data_source and project:
                    break

        if not (organism and data_source and project):
            out["BigWigURL"] = pd.NA
            return out

        file_source_col = None
        for col in out.columns:
            if col.lower().endswith("file_source"):
                file_source_col = col
                break
        if file_source_col is not None:
            first_source = out[file_source_col].dropna().iloc[0] if not out[file_source_col].dropna().empty else None
            if first_source and isinstance(first_source, str):
                resolved = first_source.strip().rstrip("/").rsplit("/", 1)[-1]
                if resolved:
                    data_source = resolved

        cfg = default_config()
        base_url = cfg.base_url.rstrip("/")

        urls: list[str | None] = []
        for ext_id in out["external_id"]:
            ext_id_str = str(ext_id).strip()
            if not ext_id_str or ext_id_str == "<NA>":
                urls.append(None)
                continue
            try:
                bw_desc = R3BigWig(
                    resource_type="bigwig_files",
                    organism=organism,
                    data_source=data_source,
                    project=project,
                    sample=ext_id_str,
                )
                urls.append(f"{base_url}/{bw_desc.url_path()}")
            except Exception:  # pylint: disable=broad-except
                urls.append(None)

        out["BigWigURL"] = urls
        return out

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
        annotation_extension: Optional[str] = None,
        assay_name: str = "raw_counts",
        join_policy: str = "inner",
        autoload: bool = True,
    ) -> summarizedexperiment.SummarizedExperiment:
        """Build a BiocPy :class:`SummarizedExperiment` from this bundle.

        This method stacks compatible count matrices, merges available
        sample metadata, and constructs a BiocPy
        :class:`SummarizedExperiment` using a compatibility-aware
        constructor that supports multiple versions of the
        :mod:`summarizedexperiment` package.

        Args:
          genomic_unit: Genomic unit to summarize, such as ``"gene"``,
            ``"exon"``, or ``"junction"``.
          annotation_extension: Optional annotation code for gene or
            exon summarizations (for example, ``"G026"``). When provided
            and ``genomic_unit`` is gene or exon, only count resources
            with matching annotation are used.
          assay_name: Name assigned to the coverage-sum assay within the
            :class:`SummarizedExperiment`. (default: ``"raw_counts"``).
          join_policy: Join policy used when concatenating counts across
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
        working = self
        if genomic_unit in {"gene", "exon"} and annotation_extension:
            working = self.filter(
                resource_type="count_files_gene_or_exon"
            ).filter(
                genomic_unit=genomic_unit,
                annotation_extension=annotation_extension,
            )

        counts_df = working._stack_counts_for(  # pylint: disable=protected-access
            genomic_unit=genomic_unit,
            join_policy=join_policy,
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

        col_df = self._add_bigwig_urls(col_df)

        original_feature_ids = list(counts_df.index)
        dup_count = int(pd.Index(original_feature_ids).duplicated().sum())

        row_names = original_feature_ids
        if dup_count:
            logging.warning(
                "Counts contain %d duplicate feature IDs; generating unique row_names "
                "and preserving originals in row_data['feature_id'].",
                dup_count,
            )
            row_names = _make_unique_names(original_feature_ids)

        counts_df = counts_df.copy()
        counts_df.index = row_names

        row_df = pd.DataFrame(
            {"feature_id": original_feature_ids},
            index=row_names,
        )

        return _construct_summarized_experiment(
            counts_df=counts_df,
            row_df=row_df,
            col_df=col_df,
            assay_name=_default_assay_name(genomic_unit, assay_name),
        )

    def to_ranged_summarized_experiment(
        self,
        *,
        genomic_unit: str,
        annotation_extension: Optional[str] = None,
        prefer_rr_junction_coordinates: bool = True,
        assay_name: str = "raw_counts",
        join_policy: str = "inner",
        autoload: bool = True,
        allow_fallback_to_se: bool = False,
    ) -> (
        summarizedexperiment.RangedSummarizedExperiment
        | summarizedexperiment.SummarizedExperiment
    ):
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
          annotation_extension: Annotation code for gene/exon
            assays, if desired.
          prefer_rr_junction_coordinates: If :data:`True`, prefer RR junction
            files for coordinate definitions when they are available.
          assay_name: Name assigned to the coverage-sum assay within the
            :class:`SummarizedExperiment` (default: ``"raw_counts"``).
          join_policy: Join policy across projects when stacking.
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
        last_ranges_error: Exception | None = None

        working = self
        counts_df = working._stack_counts_for(  # pylint: disable=protected-access
            genomic_unit=genomic_unit,
            join_policy=join_policy,
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

        col_df = self._add_bigwig_urls(col_df)

        original_feature_ids = list(counts_df.index)
        dup_count = int(pd.Index(original_feature_ids).duplicated().sum())

        row_names = original_feature_ids
        if dup_count:
            logging.warning(
                "Counts contain %d duplicate feature IDs; generating unique row_names "
                "and preserving originals in row_data['feature_id'].",
                dup_count,
            )
            row_names = _make_unique_names(original_feature_ids)

        counts_df = counts_df.copy()
        counts_df.index = row_names

        row_data_df = pd.DataFrame(
            {"feature_id": original_feature_ids},
            index=row_names,
        )

        feature_ids = original_feature_ids

        ranges_df: Optional[pd.DataFrame] = None

        if genomic_unit in {"gene", "exon"}:
            gtf_res = _select_gtf_resource_for_unit(
                self,
                genomic_unit=genomic_unit,
                annotation_extension=annotation_extension,
            )
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
                    ranges = _dedupe_ranges_on_feature_id(ranges)
                    idxed = _align_ranges_to_features(
                        ranges, feature_ids=feature_ids
                    )
                    missing_mask = (
                        idxed["seqnames"].isna()
                        | idxed["starts"].isna()
                        | idxed["ends"].isna()
                        | idxed["strand"].isna()
                    )
                    if missing_mask.any():
                        missing_ids = [
                            str(feature_ids[i])
                            for i, m in enumerate(missing_mask)
                            if m
                        ][:10]
                        raise ValueError(
                            "Annotation does not contain ranges for some feature IDs present in "
                            "the counts matrix. Example missing feature_ids: "
                            f"{missing_ids}. This usually indicates an annotation mismatch; try "
                            "setting annotation_extension to match the counts resource."
                        )

                    ranges_df = idxed[
                        ["seqnames", "starts", "ends", "strand"]
                    ].copy()
                    ranges_df.index = row_names
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

                    base_cols = ["seqnames", "starts", "ends", "strand"]
                    ranges_df = idxed[base_cols + enrich_cols].copy()
                    ranges_df.index = row_names

                    if enrich_cols:  # pragma: no branch
                        row_data_df = row_data_df.join(ranges_df[enrich_cols])
                except Exception as exc:  # pylint: disable=broad-except
                    logging.warning(
                        "Falling back: failed to parse GTF for %s ranges "
                        "(reason: %r).",
                        genomic_unit,
                        exc,
                    )
                    last_ranges_error = exc

        elif genomic_unit == "junction" and prefer_rr_junction_coordinates:
            try:
                n_features, _ = counts_df.shape

                rr_candidates = (
                    self.filter(resource_type="count_files_junctions")
                    .filter(junction_extension="RR")
                    .resources
                )
                rr_res = rr_candidates[0] if rr_candidates else None
                if rr_res is None:
                    raise ValueError(
                        "No RR junction coordinate resource found in bundle."
                    )

                if autoload:
                    rr_res.download(path=None, cache_mode="enable")

                rr = _read_rr_table(rr_res)

                std = rr.rename(
                    columns={
                        "chromosome": "seqnames",
                        "chrom": "seqnames",
                        "chr": "seqnames",
                        "start": "starts",
                        "end": "ends",
                    }
                ).copy()

                required = {"seqnames", "starts", "ends"}
                missing = required.difference(std.columns)
                if missing:
                    raise ValueError(
                        f"RR file missing required columns: {sorted(missing)}"
                    )

                if "strand" not in std.columns:
                    std["strand"] = "*"

                std["seqnames"] = std["seqnames"].astype(str)
                std["strand"] = (
                    std["strand"].astype(str).replace({"?": "*"}).fillna("*")
                )
                std["starts"] = pd.to_numeric(
                    std["starts"], errors="raise"
                ).astype(int)
                std["ends"] = pd.to_numeric(std["ends"], errors="raise").astype(
                    int
                )

                if len(std) != n_features:
                    raise ValueError(
                        f"RR row count {len(std)} != MM feature count {n_features}; cannot build junction ranges."
                    )

                id_candidates = ("junction_id", "jxn_id", "jid", "id", "name")
                id_col = next(
                    (c for c in id_candidates if c in std.columns), None
                )

                if id_col is not None:
                    row_names_from_rr = std[id_col].astype(str).tolist()
                else:
                    row_names_from_rr = [
                        f"{seq}:{start}-{end}:{strand}"
                        for seq, start, end, strand in zip(
                            std["seqnames"],
                            std["starts"],
                            std["ends"],
                            std["strand"],
                        )
                    ]

                if pd.Index(row_names_from_rr).duplicated().any():
                    row_names_from_rr = _make_unique_names(row_names_from_rr)

                counts_df = counts_df.copy()
                counts_df.index = row_names_from_rr

                ranges_df = std[["seqnames", "starts", "ends", "strand"]].copy()
                ranges_df.index = row_names_from_rr

                row_data_df = pd.DataFrame(
                    {
                        "feature_id": row_names_from_rr,
                        "mm_row": original_feature_ids,
                    },
                    index=row_names_from_rr,
                )

            except Exception as exc:  # pylint: disable=broad-except
                logging.warning(
                    "Falling back: failed to derive junction ranges from RR (reason: %r).",
                    exc,
                )
                last_ranges_error = exc
                ranges_df = None

        if ranges_df is None:
            if allow_fallback_to_se:
                return self.to_summarized_experiment(
                    genomic_unit=genomic_unit,
                    annotation_extension=annotation_extension,
                    assay_name=assay_name,
                    join_policy=join_policy,
                    autoload=autoload,
                )
            message = (
                "Could not derive genomic ranges for the requested object. "
                "Pass allow_fallback_to_se=True to receive a plain "
                "SummarizedExperiment."
            )
            if last_ranges_error is not None:
                raise ValueError(message) from last_ranges_error
            raise ValueError(message)

        return _construct_ranged_summarized_experiment(
            counts_df=counts_df,
            row_df=row_data_df,
            col_df=col_df,
            ranges_df=ranges_df,
            assay_name=_default_assay_name(genomic_unit, assay_name),
        )

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
