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
  :data:`~recount3.types.FieldSpec` for each description field. The
  following forms are accepted:

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

        pip install "recount3[biocpy]"
"""

from __future__ import annotations

from datetime import datetime, timezone
import concurrent.futures
import dataclasses
import functools
import gzip
import io
import logging
import re
from collections import Counter
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from typing import Any, Optional, TYPE_CHECKING

import pandas as pd
import numpy as np
from scipy import sparse

from recount3 import _bigwig
from recount3 import _utils
from recount3 import errors
from recount3 import resource
from recount3 import search
from recount3 import types as r3_types

if TYPE_CHECKING:  # pragma: no cover
    import biocframe
    from numpy.typing import NDArray
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
    ``external_id``) and concatenation can introduce duplicates.
    :class:`~biocframe.BiocFrame` drops duplicated names, which then breaks
    downstream validation. This
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
    out = df.copy(deep=False)
    raw_cols = [("" if c is None else str(c)) for c in out.columns]
    counts: dict[str, int] = {}
    reserved = {name or empty_prefix for name in raw_cols}
    new_cols: list[str] = []

    for name in raw_cols:
        base = name or empty_prefix
        n = counts.get(base, 0) + 1
        candidate = base if n == 1 else f"{base}__{n}"
        while n > 1 and candidate in reserved:
            n += 1
            candidate = f"{base}__{n}"
        counts[base] = n
        reserved.add(candidate)
        new_cols.append(candidate)

    out.columns = new_cols
    return out


def _default_assay_name(genomic_unit: str, assay_name: str) -> str:
    """Return the appropriate assay name for the given genomic unit.

    Args:
      genomic_unit: One of ``"gene"``, ``"exon"``, or ``"junction"``.
      assay_name: The caller-provided assay name.

    Returns:
      ``"counts"`` for junctions when the caller supplies ``"raw_counts"``;
      otherwise the caller-provided name is returned unchanged.
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

    The merge keys are rendered as text so that tables can be joined on them
    regardless of how each file's columns were typed. That text form is
    produced by :func:`recount3._utils.canonical_identifier_series` rather
    than by a plain cast, because a numeric identifier parsed as ``int64``
    in one table and ``float64`` in another would otherwise stringify to
    ``"123488"`` and ``"123488.0"`` and fail to join.

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
        out[key] = _utils.canonical_identifier_series(out[key])

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
    """Collapse duplicated keys, keeping the first non-null per column."""
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
    """Relabel count columns and sample metadata with verified external IDs.

    Args:
        counts_df: Feature-by-sample counts whose columns match the metadata
            row order.
        col_df: Sample metadata with an optional ``external_id`` column.

    Returns:
        Count and metadata frames with matching external-ID labels. The inputs
        are returned unchanged if IDs are absent, incomplete, duplicated, or
        already match. Relabeling shares count storage without modifying it.
    """
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

    renamed_counts = counts_df.copy(deep=False)
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
        path = res._cached_path()
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
    r'(?P<key>[^\s;"]+)\s+'
    r'(?:"(?P<quoted>[^"]*)"|(?P<bare>[^;"]*?))'
    r"\s*(?:;|$)"
)


def _parse_gtf_attributes(attrs: pd.Series) -> pd.DataFrame:
    """Parse a GTF attributes column into a wide (column-per-key) DataFrame.

    Values are read with GTF quoting rules, matching
    ``rtracklayer::import.gff``: a double-quoted value is taken whole, so
    ``gene_name "has;semicolon"`` keeps its semicolon and an empty
    ``gene_name ""`` yields an empty string rather than a missing one.
    Unquoted values such as ``level 2`` run to the delimiting semicolon.

    Reading quoted values as a unit is what keeps the scan aligned with the
    field: stopping at a separator inside quotes would leave the remainder
    of that value to be rescanned, and its trailing words would be parsed
    as further ``key value`` pairs and become spurious columns.

    Repeated keys retain their last value. Missing attribute strings produce
    rows with missing values; the input index is preserved.

    Args:
      attrs: Series containing raw GTF attribute strings.

    Returns:
      DataFrame indexed like attrs, with one column per attribute key.
    """
    rows = []
    for value in attrs.fillna(""):
        row = {}
        for match in _GTF_ATTR_PAIR_RE.finditer(str(value)):
            key = match.group("key")
            quoted = match.group("quoted")
            row[key] = quoted if quoted is not None else match.group("bare")
        rows.append(row)
    return pd.DataFrame(rows, index=attrs.index)


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
    """Parse covered feature lengths from the GTF score column.

    Args:
        score: GTF score values encoding covered lengths in bases. Dots and
            empty strings denote missing lengths.
        starts: One-based start coordinates, retained for call compatibility;
            these do not participate in the conversion.
        ends: Inclusive end coordinates, retained for call compatibility;
            these do not participate in the conversion.

    Returns:
        Numeric lengths with missing scores preserved. Covered gene lengths
        can differ from genomic spans, so coordinates are not used as a fallback.

    Raises:
        ValueError: If a nonmissing score cannot be converted to a number.
    """
    return pd.to_numeric(score.replace({".": pd.NA, "": pd.NA}), errors="raise")


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

    Duplicate IDs are matched by occurrence when multiplicities agree. An empty
    feature request produces a valid empty range table. Output order always
    follows the requested features.

    Args:
      ranges: DataFrame containing a 'feature_id' column and coordinate columns.
      feature_ids: Feature IDs from the counts matrix (in desired order).

    Returns:
      DataFrame reindexed to feature_ids order, containing the same columns as
      `ranges` except for 'feature_id' (which is used as the index).
      Unmatched unique IDs retain missing coordinates for the caller to handle.

    Raises:
      ValueError: If required coordinate columns are missing from `ranges`,
        or if version-stripped matching is ambiguous (conflicting
        coordinates).
      recount3.errors.RangesCoverageError: If duplicate occurrences cannot
        be matched completely and unambiguously.
    """
    required = {"feature_id", "seqnames", "starts", "ends", "strand"}
    missing_cols = required - set(ranges.columns)
    if missing_cols:
        raise ValueError(
            "ranges is missing required columns: " f"{sorted(missing_cols)}"
        )

    feature_index = pd.Index([str(x) for x in feature_ids])

    ids = pd.Index(ranges["feature_id"].astype(str))
    if ids.equals(feature_index):
        exact = ranges.drop(columns="feature_id").copy(deep=False)
        exact.index = feature_index
        return exact
    if ids.has_duplicates or feature_index.has_duplicates:
        available = pd.Series(ids).value_counts()
        requested = pd.Series(feature_index).value_counts()
        for name, count in requested.items():
            if (
                name in available
                and available[name] != count
                and (count > 1 or available[name] > 1)
            ):
                raise errors.RangesCoverageError(
                    f"Ambiguous duplicate feature occurrences for {name!r}."
                )
        source_index = pd.MultiIndex.from_arrays(
            [ids, pd.Series(ids).groupby(pd.Series(ids), sort=False).cumcount()]
        )
        target_index = pd.MultiIndex.from_arrays(
            [
                feature_index,
                pd.Series(feature_index)
                .groupby(pd.Series(feature_index), sort=False)
                .cumcount(),
            ]
        )
        exact = (
            ranges.drop(columns="feature_id")
            .set_axis(source_index)
            .reindex(target_index)
        )
        exact.index = feature_index
        if exact[["seqnames", "starts", "ends", "strand"]].isna().any().any():
            raise errors.RangesCoverageError(
                "Cannot match duplicate feature occurrences to annotation."
            )
        return exact
    exact = ranges.set_index("feature_id").reindex(feature_index)

    coord_cols = ["seqnames", "starts", "ends", "strand"]
    missing_any = exact[coord_cols].isna().any(axis=1)
    if not missing_any.any():
        return exact

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
                "After stripping Ensembl versions, multiple annotation "
                "rows map to the same ID but with conflicting "
                "coordinates. Example IDs: "
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
    filled = exact.combine_first(fallback).reindex(feature_index)

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
        path = res._cached_path()
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
    """Extract feature coordinates and annotation metadata from a GTF table.

    Args:
        gtf: Parsed GTF columns including ``feature``, ``seqname``, ``start``,
            ``end``, ``strand``, ``source``, ``score``, ``frame``, and
            ``attributes``.
        feature_kind: Feature type to retain: ``"gene"`` or ``"exon"``.

    Returns:
        A table with one-based inclusive coordinates, ``feature_id``, source,
        feature type, covered ``bp_length``, nullable phase, and parsed GTF
        attributes. Gene IDs come from ``gene_id``; exon IDs prefer
        ``recount_exon_id`` over ``exon_id``. Missing IDs fall back to coordinate
        strings. Duplicate IDs with identical ranges retain separate annotation
        rows. No matching features produces an empty coordinate/ID table.

    Raises:
        KeyError: If a required GTF column is missing.
        ValueError: If coordinates or covered lengths cannot be parsed, or
            duplicate feature IDs describe conflicting genomic ranges.
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

    return out


def _classify_ranges_failure(
    exc: Optional[BaseException], *, source: str = "annotation"
) -> str:
    """Describe why genomic ranges could not be derived.

    Four outcomes are worth telling apart, because the user action differs
    for each: nothing to derive ranges from was in the bundle; the file
    could not be retrieved (worth retrying, or a mirror problem); it was
    retrieved but could not be parsed (damaged or unexpected content); or
    it parsed cleanly but does not describe every counted feature (a
    mismatch, fixed by selecting the matching ``annotation_extension``).

    Args:
      exc: The exception raised while deriving ranges, or :data:`None`
        when no attempt was made at all.
      source: What ranges were to be read from, for the returned phrase --
        ``"annotation"`` for gene and exon units, the RR coordinate file
        for junctions.

    Returns:
      A short human-readable phrase naming the outcome.
    """
    if exc is None or isinstance(exc, errors.MissingRangesError):
        return f"no {source} providing genomic ranges was in the bundle"
    if isinstance(exc, errors.RangesCoverageError):
        return f"the {source} does not cover every counted feature"
    if isinstance(exc, (gzip.BadGzipFile, EOFError)):
        return f"the {source} could not be parsed"
    if isinstance(exc, (errors.DownloadError, OSError)):
        return f"the {source} could not be retrieved"
    return f"the {source} could not be parsed"


_RR_SOURCE = "RR coordinate file"


def _peek_gtf_feature_counts(
    res: resource.R3Resource,
    *,
    max_lines: int = 50000,
    autoload: bool = True,
) -> Counter[str]:
    """Scan the GTF and count feature types without loading it fully.

    Args:
      res: The annotation resource to inspect.
      max_lines: Stop after this many feature rows.
      autoload: If :data:`True`, download the annotation when it is not
        already cached. If :data:`False`, an uncached annotation raises
        :exc:`FileNotFoundError` instead of being fetched.

    Returns:
      A :class:`collections.Counter` mapping GTF feature type to the number
      of rows seen for it within the scanned prefix of the file.
    """
    path = res.ensure_cached(download=autoload)

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
    autoload: bool = True,
) -> Optional[resource.R3Resource]:
    """Pick the most appropriate annotation resource for gene/exon ranges.

    An unambiguous matching descriptor is selected without file I/O.
    Otherwise candidates are ranked by description and URL and confirmed by
    peeking at the file itself. Confirmation needs the annotation on disk,
    so with ``autoload`` enabled a candidate that is not cached yet is
    downloaded first; with ``autoload`` disabled, uncached candidates are
    skipped and the ranking alone decides.

    Args:
      bundle: The bundle to search for annotation resources.
      genomic_unit: Either ``"gene"`` or ``"exon"``.
      annotation_extension: Restrict candidates to this annotation code
        when given.
      autoload: If :data:`True`, download candidate annotations as needed
        to inspect them. If :data:`False`, inspect only what is cached.

    Returns:
      The selected resource, or :data:`None` when the bundle holds no
      annotation resources at all.
    """
    feature_kind = "gene" if genomic_unit == "gene" else "exon"

    ann = bundle.filter(resource_type="annotations")
    if annotation_extension:
        ann = ann.filter(annotation_extension=annotation_extension)

    candidates = list(ann.resources)
    exact = [
        r
        for r in candidates
        if _description_value(r, "genomic_unit") == genomic_unit
    ]
    if len(exact) == 1:
        return exact[0]
    if len(exact) > 1:
        raise errors.CompatibilityError(
            "Multiple matching annotations; select one annotation resource."
        )
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
            feats = _peek_gtf_feature_counts(
                res, max_lines=50000, autoload=autoload
            )
        except FileNotFoundError:
            # Expected whenever autoload is off and the candidate has not
            # been downloaded yet; the ranking above still applies.
            logging.debug(
                "Not inspecting annotation %s: it is not cached and "
                "autoload is disabled.",
                res.url,
            )
            continue
        except Exception as exc:  # pylint: disable=broad-exception-caught
            logging.warning(
                "Skipping annotation candidate %s: could not inspect its "
                "GTF features (%r).",
                res.url,
                exc,
            )
            continue

        if feats.get(feature_kind, 0) > 0:
            logging.info(
                "Selected annotation resource for %s: %s (peek features=%s)",
                genomic_unit,
                res.description.url_path(),
                dict(feats.most_common(5)),
            )
            return res

    return candidates[0]


def _to_genomic_ranges(ranges_df: pd.DataFrame) -> genomicranges.GenomicRanges:
    """Construct a :class:`~genomicranges.GenomicRanges` from a DataFrame.

    Args:
      ranges_df: DataFrame with at least ``seqnames``, ``starts``,
        ``ends``, and ``strand`` columns using one-based inclusive coordinates.
        The index supplies range names; other columns become range metadata.

    Returns:
      A :class:`genomicranges.GenomicRanges` instance.

    Raises:
      ImportError: If required BiocPy packages cannot be imported.
      KeyError: If a required coordinate column is absent.
      ValueError: If the range constructor rejects the coordinates.
    """
    genomic_ranges_cls = _utils.get_genomicranges_class()
    coords = ["seqnames", "starts", "ends", "strand"]
    gr = genomic_ranges_cls.from_pandas(
        ranges_df[coords].reset_index(drop=True)
    )
    extras = {
        col: ranges_df[col].to_numpy() for col in ranges_df if col not in coords
    }
    if extras:
        gr = gr.set_mcols(_utils.get_biocframe_class()(extras))
    return gr.set_names([str(x) for x in ranges_df.index])


def _numeric_assay(
    counts_df: pd.DataFrame, *, copy: bool = True
) -> NDArray[Any] | sparse.csc_matrix:
    """Convert count data to a validated dense array or CSC sparse matrix.

    Args:
        counts_df: Feature-by-sample counts. Sparse columns must all use zero
            fill values. Dense numeric strings are converted when necessary.
        copy: Whether to copy assay storage. With ``False``, dense output may
            share storage with the input; sparse conversion may still allocate.

    Returns:
        A two-dimensional numeric array, or a CSC matrix when all input columns
        are sparse. Numeric dtypes are preserved when conversion is unnecessary.

    Raises:
        ValueError: If counts are nonnumeric, missing, nonfinite, negative,
            not two-dimensional, or use a nonzero sparse fill value.
    """
    if not isinstance(counts_df, pd.DataFrame) or counts_df.ndim != 2:
        raise ValueError("counts_df must be a 2D DataFrame.")
    is_sparse = len(counts_df.columns) > 0 and all(
        isinstance(dtype, pd.SparseDtype) for dtype in counts_df.dtypes
    )
    if is_sparse:
        if any(dtype.fill_value != 0 for dtype in counts_df.dtypes):
            raise ValueError("Sparse counts must have zero fill values.")
        matrix = sparse.csc_matrix(counts_df.sparse.to_coo(), copy=copy)
        values = matrix.data
    else:
        if not all(pd.api.types.is_numeric_dtype(d) for d in counts_df.dtypes):
            try:
                counts_df = counts_df.apply(pd.to_numeric, errors="raise")
            except (ValueError, TypeError) as exc:
                raise ValueError(
                    "counts_df contains non-numeric values."
                ) from exc
        matrix = counts_df.to_numpy(copy=copy)
        if matrix.dtype.kind == "O":
            matrix = counts_df.to_numpy(dtype=float, na_value=np.nan, copy=copy)
        values = matrix
    if not np.isfinite(values).all():
        raise ValueError("counts_df contains missing or non-finite counts.")
    if (values < 0).any():
        raise ValueError("counts_df contains negative counts.")
    return matrix


def _validate_count_frame(frame: pd.DataFrame) -> None:
    """Validate count columns without materializing a complete second assay.

    Sparse columns are checked through their stored values, preserving implicit
    zeros. Validation does not modify the input frame.

    Args:
        frame: Feature-by-sample count table with unique column labels.

    Raises:
        ValueError: If a column contains nonnumeric, missing, nonfinite, or
            negative counts, or has a nonzero sparse fill value.
    """
    for name in frame.columns:
        column = frame[name]
        if isinstance(column.dtype, pd.SparseDtype):
            if column.dtype.fill_value != 0:
                raise ValueError("Sparse counts must have zero fill values.")
            values = column.array.sp_values
        else:
            try:
                numeric = pd.to_numeric(column, errors="raise")
                values = numeric.to_numpy()
                if values.dtype.kind == "O":
                    values = numeric.to_numpy(dtype=float, na_value=np.nan)
            except (ValueError, TypeError) as exc:
                raise ValueError(
                    "counts_df contains non-numeric values."
                ) from exc
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(
                "Counts must be finite, non-missing and nonnegative."
            )


def _experiment_frame(
    frame: pd.DataFrame, names: list[str], label: str
) -> biocframe.BiocFrame:
    """Align metadata with an assay dimension and construct a BiocFrame.

    Args:
        frame: Feature or sample metadata. A RangeIndex denotes positional
            alignment; other indices must identify the requested rows uniquely.
        names: Assay identifiers in the required order.
        label: Metadata dimension name used in errors and empty column labels.

    Returns:
        Metadata with the requested row names, unique nonempty column names,
        and copied column arrays independent of the input.

    Raises:
        ValueError: If metadata length or identifiers do not match the assay.
        ImportError: If BiocFrame is unavailable.
    """
    if len(frame) != len(names):
        raise ValueError(
            f"{label} length {len(frame)} != assay dimension {len(names)}."
        )
    frame = frame.copy(deep=False)
    if not isinstance(frame.index, pd.RangeIndex):
        frame.index = frame.index.map(str)
        if list(frame.index) != names:
            if not frame.index.is_unique or set(frame.index) != set(names):
                raise ValueError(
                    f"{label} identifiers do not match assay identifiers."
                )
            frame = frame.reindex(names)
    frame.index = names
    frame = _ensure_unique_columns(frame, empty_prefix=label)
    return _utils.get_biocframe_class()(
        {col: frame[col].to_numpy(copy=True) for col in frame.columns},
        row_names=names,
        number_of_rows=len(names),
    )


def _construct_summarized_experiment(
    *,
    counts_df: pd.DataFrame,
    row_df: pd.DataFrame,
    col_df: pd.DataFrame,
    assay_name: str,
    metadata: Mapping[str, Any] | None = None,
) -> summarizedexperiment.SummarizedExperiment:
    """Construct a validated experiment with independently owned assay storage.

    Args:
        counts_df: Feature-by-sample counts with feature and sample identifiers
            on the index and columns. Sparse inputs remain sparse.
        row_df: Feature metadata aligned by identifier, or position for a
            RangeIndex.
        col_df: Sample metadata aligned by identifier, or position for a
            RangeIndex.
        assay_name: Name assigned to the count assay.
        metadata: Experiment-level annotations, copied into a new dictionary.

    Returns:
        A SummarizedExperiment with synchronized names, metadata, and a dense
        or CSC sparse count assay.

    Raises:
        ValueError: If counts are invalid or metadata cannot be aligned.
        TypeError: If the experiment constructor rejects an input type.
        ImportError: If the required BiocPy packages are unavailable.
    """
    matrix = _numeric_assay(counts_df)
    rows, cols = list(counts_df.index.map(str)), list(
        counts_df.columns.map(str)
    )
    return _utils.get_summarizedexperiment_class()(
        assays={assay_name: matrix},
        row_data=_experiment_frame(row_df, rows, "row_data"),
        column_data=_experiment_frame(col_df, cols, "column_data"),
        row_names=rows,
        column_names=cols,
        metadata=dict(metadata or {}),
    )


def _construct_ranged_summarized_experiment(
    *,
    counts_df: pd.DataFrame,
    row_df: pd.DataFrame,
    col_df: pd.DataFrame,
    ranges_df: pd.DataFrame,
    assay_name: str,
    metadata: Mapping[str, Any] | None = None,
) -> summarizedexperiment.RangedSummarizedExperiment:
    """Construct a validated experiment with genomic ranges and owned counts.

    Args:
        counts_df: Feature-by-sample counts with named rows and columns.
            Sparse inputs remain sparse.
        row_df: Feature metadata aligned by identifier, or position for a
            RangeIndex.
        col_df: Sample metadata aligned by identifier, or position for a
            RangeIndex.
        ranges_df: One row per feature with ``seqnames``, ``starts``, ``ends``,
            and ``strand`` columns. Coordinates are one-based and inclusive;
            additional columns become range metadata. A RangeIndex denotes
            positional alignment; other indices must match feature identifiers.
        assay_name: Name assigned to the count assay.
        metadata: Experiment-level annotations, copied into a new dictionary.

    Returns:
        A RangedSummarizedExperiment with aligned feature ranges, metadata,
        names, and a dense or CSC sparse count assay.

    Raises:
        ValueError: If counts, coordinates, dimensions, or identifiers are
            invalid or required range columns are absent.
        TypeError: If the experiment constructor rejects an input type.
        ImportError: If the required BiocPy packages are unavailable.
    """
    matrix = _numeric_assay(counts_df)
    rows, cols = list(counts_df.index.map(str)), list(
        counts_df.columns.map(str)
    )
    required = {"seqnames", "starts", "ends", "strand"}
    if required - set(ranges_df):
        raise ValueError(
            f"ranges_df is missing required columns: {sorted(required-set(ranges_df))}."
        )
    if len(ranges_df) != len(rows):
        raise ValueError(
            f"ranges_df length {len(ranges_df)} != assay features {len(rows)}."
        )
    ranges = ranges_df.copy(deep=False)
    if (
        not isinstance(ranges.index, pd.RangeIndex)
        and list(ranges.index.map(str)) != rows
    ):
        if not ranges.index.is_unique or set(ranges.index.map(str)) != set(
            rows
        ):
            raise ValueError(
                "ranges_df identifiers do not match assay identifiers."
            )
        ranges = ranges.reindex(rows)
    ranges.index = rows
    _validate_coordinates(ranges)
    return _utils.get_ranged_summarizedexperiment_class()(
        assays={assay_name: matrix},
        row_ranges=_to_genomic_ranges(ranges),
        row_data=_experiment_frame(row_df, rows, "row_data"),
        column_data=_experiment_frame(col_df, cols, "column_data"),
        row_names=rows,
        column_names=cols,
        metadata=dict(metadata or {}),
    )


def _validate_coordinates(frame: pd.DataFrame) -> None:
    """Validate one-based inclusive coordinates and normalized strand values.

    Args:
        frame: Range table containing ``seqnames``, ``starts``, ``ends``, and
            ``strand``. Starts and ends must be finite integers with
            ``1 <= starts <= ends``; strands must be ``+``, ``-``, or ``*``.

    Raises:
        KeyError: If a required coordinate column is absent.
        ValueError: If coordinates or strands are missing or invalid.
    """
    if frame[["seqnames", "starts", "ends", "strand"]].isna().any().any():
        raise ValueError(
            "ranges_df contains missing values in coordinate columns."
        )
    for col in ("starts", "ends"):
        values = pd.to_numeric(frame[col], errors="raise").to_numpy(dtype=float)
        if not np.isfinite(values).all() or (values != np.floor(values)).any():
            raise ValueError("Genomic coordinates must be finite integers.")
    if (pd.to_numeric(frame["starts"]) < 1).any() or (
        pd.to_numeric(frame["ends"]) < pd.to_numeric(frame["starts"])
    ).any():
        raise ValueError(
            "Genomic coordinates must be positive inclusive intervals."
        )
    if not frame["strand"].isin(["+", "-", "*"]).all():
        raise ValueError("Invalid genomic strand.")


def _description_value(res: resource.R3Resource, key: str) -> str | None:
    """Read a string-valued resource description field.

    Args:
        res: Resource whose description supplies the field.
        key: Description attribute name.

    Returns:
        The field value if it is a string, otherwise ``None``.
    """
    value = getattr(res.description, key, None)
    return value if isinstance(value, str) else None


def _project_key(
    res: resource.R3Resource,
) -> tuple[str | None, str | None, str | None]:
    """Identify the organism, source, and project associated with a resource.

    Args:
        res: Resource to identify.

    Returns:
        An ``(organism, data_source, project)`` tuple. Missing or nonstring
        fields are represented by ``None``.
    """
    return tuple(
        _description_value(res, key)
        for key in ("organism", "data_source", "project")
    )


def _merge_count_frames(
    frames: Sequence[pd.DataFrame], join_policy: str
) -> pd.DataFrame:
    """Combine count tables by feature identity while preserving sample order.

    Args:
        frames: Nonempty sequence of validated feature-by-sample count tables.
            Repeated feature IDs are allowed only when all feature indices
            have exactly the same ordering.
        join_policy: ``"inner"`` intersects feature IDs; ``"outer"`` takes
            their union and fills newly introduced feature rows with zeros.

    Returns:
        Counts with columns concatenated in input order. Existing missing
        values are preserved. A single input returns a shallow copy.

    Raises:
        ValueError: If the join policy is unsupported.
        recount3.errors.CompatibilityError: If duplicate feature IDs prevent
            unambiguous alignment across different feature indices.
    """
    if join_policy not in ("inner", "outer"):
        raise ValueError("join_policy must be 'inner' or 'outer'.")
    if len(frames) == 1:
        return frames[0].copy(deep=False)
    if all(frame.index.equals(frames[0].index) for frame in frames[1:]):
        return pd.concat(frames, axis=1)
    if any(not frame.index.is_unique for frame in frames):
        raise errors.CompatibilityError(
            "Cannot align duplicate feature IDs across projects."
        )
    index = frames[0].index
    for frame in frames[1:]:
        index = (
            index.intersection(frame.index, sort=False)
            if join_policy == "inner"
            else index.union(frame.index, sort=False)
        )
    return pd.concat(
        [frame.reindex(index, fill_value=0) for frame in frames], axis=1
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
    match rtype:
        case "count_files_gene_or_exon":
            genomic_unit = getattr(res.description, "genomic_unit", None) or ""
            family = "gene_or_exon"
            feature_key = f"{family}:{genomic_unit}"
            return family, feature_key
        case "count_files_junctions":
            junction_type = (
                getattr(res.description, "junction_type", None) or ""
            )
            junction_ext = (
                getattr(res.description, "junction_extension", None) or ""
            )
            family = "junctions"
            feature_key = f"{family}:{junction_type}:{junction_ext}"
            return family, feature_key
        case _:
            raise ValueError(
                "Resource is not a recognized count-file type for stacking: "
                f"{rtype!r}"
            )


def _make_unique_names(
    names: Sequence[str],
    *,
    suffix: str = "__dup",
) -> list[str]:
    """Make feature names unique while preserving their order.

    Args:
        names: Original feature names, including possible duplicates.
        suffix: Text inserted before the occurrence number on duplicates.

    Returns:
        Names with later occurrences suffixed starting at two. Generated names
        avoid collisions with both original and previously generated names.

    Examples:
        >>> _make_unique_names(["a", "b", "a"])
        ['a', 'b', 'a__dup2']
    """
    seen: Counter[str] = Counter()
    reserved = set(names)
    out: list[str] = []
    for name in names:
        seen[name] += 1
        if seen[name] == 1:
            out.append(name)
        else:
            candidate = f"{name}{suffix}{seen[name]}"
            while candidate in reserved:
                seen[name] += 1
                candidate = f"{name}{suffix}{seen[name]}"
            out.append(candidate)
            reserved.add(candidate)
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
    _range_cache: Any = dataclasses.field(
        default=None, repr=False, compare=False
    )

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

        It can operate on a single ``(organism, data_source, project)``
        triple or on the Cartesian product of multiple values for each
        identifier.

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
                    project=["SRP009615", "SRP001558"],
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

        This is a convenience alias for :meth:`only_counts`.
        """
        return self.only_counts()

    def metadata(self) -> R3ResourceBundle:
        """Return a sub-bundle containing only metadata resources.

        This is a convenience alias for :meth:`only_metadata`.
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

            Require an identical feature space; fails if gene and exon are
            mixed (the annotation build is not constrained)::

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

        match compat:
            case "family":
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
                        "bundle.filter("
                        'resource_type="count_files_gene_or_exon") '
                        "or bundle.filter("
                        'resource_type="count_files_junctions").'
                    )
            case "feature":
                if len(features) > 1:
                    examples = ", ".join(sorted(features))
                    raise errors.CompatibilityError(
                        "Feature-level incompatibility detected. All inputs "
                        "must share the same feature key (for example, gene "
                        "vs exon; junction subtype). Distinct feature keys "
                        f"observed: {examples}. Hint: filter by 'genomic_unit' "
                        "for gene/exon or by 'junction_type' / "
                        "'junction_extension' for junctions."
                    )
            case _:
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

        This is a bundle-scoped helper used by the
        :class:`~summarizedexperiment.SummarizedExperiment` builders. It
        enforces appropriate compatibility within the gene/exon family or
        the junctions family.

        Args:
          genomic_unit: One of ``"gene"``, ``"exon"``, or ``"junction"``.
          join_policy: ``"inner"`` intersects feature rows; ``"outer"`` unions them
            and fills newly introduced rows with zero.
          autoload: If :data:`True`, load resources on demand.

        Returns:
          A wide matrix of shape ``(features, samples)`` as a
          :class:`pandas.DataFrame`.

        Raises:
          ValueError: If no compatible count resources exist or stacking
            fails, including a required count resource failing to load.
          recount3.errors.CompatibilityError: If feature spaces are incompatible.
          recount3.errors.RangesError: If multi-project junction coordinates are
            missing or ambiguous.
        """
        if join_policy not in ("inner", "outer"):
            raise ValueError("join_policy must be 'inner' or 'outer'.")
        unit = _utils._normalize_genomic_unit(genomic_unit)
        if unit == "junction":
            selected = self.filter(
                resource_type="count_files_junctions", junction_extension="MM"
            )
        else:
            selected = self.filter(
                resource_type="count_files_gene_or_exon", genomic_unit=unit
            )
        if not selected.resources:
            raise ValueError(
                "No count-file resources available for the requested genomic unit."
            )
        feature_key = (
            "junction_type" if unit == "junction" else "annotation_extension"
        )
        compatibility = {
            (
                _description_value(res, "organism"),
                _description_value(res, feature_key),
            )
            for res in selected.resources
        }
        if len(compatibility) > 1:
            raise errors.CompatibilityError(
                "Selected counts have incompatible organisms, annotations or junction formats."
            )
        frames = []
        for res in selected.resources:
            if not res.is_loaded():
                if not autoload:
                    raise ValueError(
                        f"Count resource is not loaded: {res.url}. "
                        "Load it first or use autoload=True."
                    )
                try:
                    res.load()
                except Exception as exc:
                    raise ValueError(
                        f"Failed to load requested count matrix {res.url}: {exc}"
                    ) from exc
            frame = res.get_loaded()
            if not isinstance(frame, pd.DataFrame):
                raise TypeError(f"Loaded counts are not a DataFrame: {res.url}")
            if frame.columns.has_duplicates:
                raise ValueError(f"Duplicate sample identifiers in {res.url}")
            for names in (frame.index, frame.columns):
                if names.isna().any() or any(not str(x).strip() for x in names):
                    raise ValueError(
                        f"Missing feature or sample identifiers in {res.url}"
                    )
            _validate_count_frame(frame)
            frame = frame.copy(deep=False)
            frame.index = frame.index.map(str)
            frame.columns = frame.columns.map(str)
            if unit == "junction" and len(selected.resources) > 1:
                group = R3ResourceBundle(
                    [
                        r
                        for r in self.resources
                        if _project_key(r) == _project_key(res)
                    ]
                )
                ranges = group._junction_ranges(frame, autoload=autoload)
                frame.index = ranges.index
            frames.append(frame)
        return _merge_count_frames(frames, join_policy)

    def _junction_ranges(
        self, counts: pd.DataFrame, *, autoload: bool
    ) -> pd.DataFrame:
        """Resolve the coordinate sidecar for a single junction count matrix.

        Loaded sidecars are reused; newly parsed sidecars are cached on their
        resource. Missing and unknown strands are normalized to ``*``.

        Args:
            counts: Count table whose row order and length match the matrix's
                coordinate sidecar.
            autoload: Whether to retrieve the sidecar if it is not cached locally.

        Returns:
            Range metadata indexed by ``chromosome:start-end:strand`` identities,
            with one-based inclusive integer coordinates.

        Raises:
            recount3.errors.CompatibilityError: If the bundle does not contain
                exactly one junction matrix.
            recount3.errors.MissingRangesError: If exactly one sidecar matching
                the matrix's project and junction format cannot be identified.
            recount3.errors.RangesCoverageError: If required columns are missing,
                row counts differ, or junction coordinates are duplicated.
            ValueError: If coordinates are invalid.
            recount3.errors.DownloadError: If sidecar retrieval fails.
            OSError: If the cached sidecar cannot be read.
        """
        mm = self.filter(
            resource_type="count_files_junctions", junction_extension="MM"
        ).resources
        rr = self.filter(
            resource_type="count_files_junctions", junction_extension="RR"
        ).resources
        if len(mm) != 1:
            raise errors.CompatibilityError(
                "Each junction matrix needs its own project and RR sidecar."
            )
        rr = [
            res
            for res in rr
            if _project_key(res) == _project_key(mm[0])
            and _description_value(res, "junction_type")
            == _description_value(mm[0], "junction_type")
        ]
        if len(rr) != 1:
            raise errors.MissingRangesError(
                "Expected exactly one matching RR junction coordinate resource."
            )
        res = rr[0]
        if res.is_loaded() and isinstance(res.get_loaded(), pd.DataFrame):
            frame = res.get_loaded()
        else:
            res.ensure_cached(download=autoload)
            frame = pd.read_table(res._cached_path(), compression="infer")
            res._cached_data = frame
        frame = frame.rename(
            columns={
                "chromosome": "seqnames",
                "chrom": "seqnames",
                "chr": "seqnames",
                "start": "starts",
                "end": "ends",
            }
        ).copy()
        if {"seqnames", "starts", "ends"} - set(frame):
            raise errors.RangesCoverageError(
                "RR file missing required columns."
            )
        if "strand" not in frame:
            frame["strand"] = "*"
        frame["strand"] = frame["strand"].replace({"?": "*"}).fillna("*")
        if len(frame) != len(counts):
            raise errors.RangesCoverageError(
                f"RR row count {len(frame)} != MM feature count {len(counts)}."
            )
        _validate_coordinates(frame)
        frame["starts"] = pd.to_numeric(frame["starts"]).astype("int64")
        frame["ends"] = pd.to_numeric(frame["ends"]).astype("int64")
        names = [
            f"{seq}:{start}-{end}:{strand}"
            for seq, start, end, strand in zip(
                frame["seqnames"],
                frame["starts"],
                frame["ends"],
                frame["strand"],
            )
        ]
        if pd.Index(names).has_duplicates:
            raise errors.RangesCoverageError(
                "RR contains duplicate junction coordinates."
            )
        frame.index = names
        return frame

    def _add_bigwig_urls(
        self,
        col_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Add a ``BigWigURL`` column to sample metadata.

        For each sample (identified by ``external_id`` in ``col_df``), this
        method constructs the recount3 mirror URL for its BigWig coverage file
        by matching each sample's ``study`` to a count resource. A bundle with
        only one project supplies the project identity when ``study`` is absent.
        Sample ``file_source`` metadata overrides the source when available;
        each resource's configuration supplies the mirror base URL.

        Args:
          col_df: Sample metadata DataFrame with an ``external_id`` column.

        Returns:
          A copy of ``col_df`` with a ``BigWigURL`` column appended. If
          the necessary resource attributes cannot be inferred, the column
          contains ``None`` for the affected samples.
        """
        out = col_df.copy()
        count_types = {"count_files_gene_or_exon", "count_files_junctions"}
        candidates = [
            r
            for r in self.resources
            if getattr(r.description, "resource_type", None) in count_types
        ]
        urls = []
        for _, row in out.iterrows():
            study = row.get("study")
            matching = [
                r
                for r in candidates
                if pd.notna(study)
                and _description_value(r, "project") == str(study)
            ]
            if not matching and len({_project_key(r) for r in candidates}) == 1:
                matching = candidates[:1]
            res = matching[0] if matching else None
            external = row.get("external_id")
            if res is None or pd.isna(external) or not str(external).strip():
                urls.append(None)
                continue
            org, source, project = _project_key(res)
            source_fields = [
                c for c in out.columns if c.lower().endswith("file_source")
            ]
            if source_fields and isinstance(row[source_fields[0]], str):
                source = row[source_fields[0]].rstrip("/").rsplit("/", 1)[-1]
            if not all((org, source, project)):
                urls.append(None)
                continue
            urls.append(
                resource.build_url(
                    "bigwig_files",
                    organism=org,
                    data_source=source,
                    project=project,
                    sample=str(external),
                    config=res.config,
                )
            )
        out["BigWigURL"] = urls
        return out

    def _normalize_sample_metadata(
        self,
        *,
        sample_ids: Sequence[str],
        autoload: bool = True,
        metadata_join: str = "inner",
    ) -> pd.DataFrame:
        """Merge metadata within each project and align it to count samples.

        Empty tables are ignored when another nonempty table is available.
        External and rail identifiers must define consistent sample mappings.

        Args:
            sample_ids: Count sample identifiers in the desired output order.
            autoload: Whether to load metadata resources that are not already
                in memory.
            metadata_join: ``"inner"`` intersects nonempty tables per project;
                ``"outer"`` retains every count sample with missing metadata
                where no match exists.

        Returns:
            Metadata indexed by selected count sample IDs, with namespaced columns
            and column provenance in ``attrs``. A metadata-free bundle returns
            only an ``external_id`` column for all requested samples.

        Raises:
            ValueError: If the join policy is invalid, all supplied tables are
                empty, identifiers conflict, table columns overlap, required
                samples are absent, or unloaded resources cannot be autoloaded.
            TypeError: If a loaded metadata resource is not a DataFrame.
            recount3.errors.LoadError: If a metadata resource cannot be loaded.
            recount3.errors.DownloadError: If metadata retrieval fails.
        """
        if metadata_join not in ("inner", "outer"):
            raise ValueError("metadata_join must be 'inner' or 'outer'.")
        sample_ids = [str(x) for x in sample_ids]
        groups = {}
        provenance = {}
        for res in self.only_metadata().resources:
            if not res.is_loaded():
                if not autoload:
                    raise ValueError(
                        f"Metadata resource is not loaded: {res.url}"
                    )
                res.load()
            frame = res.get_loaded()
            if not isinstance(frame, pd.DataFrame):
                raise TypeError(
                    f"Metadata resource is not a DataFrame: {res.url}"
                )
            if len(frame) == 0:
                logging.warning("Dropping empty metadata table %s", res.url)
                continue
            frame = _standardize_metadata_frame(frame)
            frame, origin = _namespace_metadata_columns(
                frame, origin=_metadata_origin(res)
            )
            provenance.update(origin)
            if frame.duplicated(list(_METADATA_MERGE_KEYS)).any():
                raise ValueError(f"Duplicate sample metadata keys in {res.url}")
            groups.setdefault(_project_key(res), []).append(frame)
        if not groups:
            result = pd.DataFrame({"external_id": sample_ids}, index=sample_ids)
            if self.only_metadata().resources:
                raise ValueError("All supplied metadata tables are empty.")
            return result
        merged_groups = []
        for frames in groups.values():
            identities = pd.concat(
                [f[list(_METADATA_MERGE_KEYS)] for f in frames]
            ).drop_duplicates()
            for key, other in (
                ("rail_id", "external_id"),
                ("external_id", "rail_id"),
            ):
                if (
                    identities.dropna(subset=[key, other])
                    .groupby(key)[other]
                    .nunique()
                    > 1
                ).any():
                    raise ValueError(
                        f"Conflicting {key}/{other} mappings in metadata."
                    )
            merged = frames[0]
            for frame in frames[1:]:
                overlap = set(merged.columns).intersection(frame.columns) - set(
                    _METADATA_MERGE_KEYS
                )
                if overlap:
                    raise ValueError(
                        f"Repeated metadata table columns: {sorted(overlap)}"
                    )
                merged = pd.merge(
                    merged,
                    frame,
                    on=list(_METADATA_MERGE_KEYS),
                    how=metadata_join,
                    sort=False,
                    validate="one_to_one",
                )
            merged_groups.append(merged)
        merged = pd.concat(merged_groups, ignore_index=True)
        lookup = {}
        for position, row in merged.iterrows():
            for key in ("external_id", "rail_id"):
                value = row[key]
                if pd.notna(value):
                    name = str(value)
                    if name in lookup and lookup[name] != position:
                        raise ValueError(
                            f"Ambiguous sample identifier {name!r} in metadata."
                        )
                    lookup[name] = position
        if metadata_join == "inner":
            represented = {
                lookup[name] for name in sample_ids if name in lookup
            }
            if len(represented) != len(merged):
                raise ValueError(
                    "Metadata contains samples missing from the counts matrix."
                )
            if sample_ids and not represented:
                # len(represented) == len(merged) also holds when both are
                # zero, so the check above cannot see a join that matched
                # nothing. Without this, every sample would be dropped and
                # the caller would receive an experiment with no columns.
                raise ValueError(
                    "No count sample matched the merged sample metadata. "
                    "The metadata tables for this project share no "
                    f"{'/'.join(_METADATA_MERGE_KEYS)} rows, so the inner "
                    "join produced no samples. Check that the tables "
                    "describe the same samples, or pass "
                    "metadata_join='outer' to keep every count sample and "
                    "leave the unmatched metadata missing."
                )
        selected = [
            name
            for name in sample_ids
            if name in lookup or metadata_join == "outer"
        ]
        positions = [lookup.get(name, -1) for name in selected]
        aligned = merged.reindex(positions).copy()
        aligned.index = selected
        missing = aligned["external_id"].isna()
        if missing.any():
            if metadata_join == "inner":
                raise ValueError(
                    "Matched metadata rows have missing external_id values."
                )
            aligned.loc[missing, "external_id"] = pd.Series(
                selected, index=selected
            )[missing]
        aligned.attrs["recount3_metadata_provenance"] = provenance
        return aligned

    def _prepare_experiment(
        self,
        *,
        genomic_unit: str,
        annotation_extension: str | None,
        join_policy: str,
        metadata_join: str,
        autoload: bool,
    ) -> tuple[
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame | None,
        dict[str, Any],
    ]:
        """Prepare aligned counts, metadata, and provenance for either builder.

        Resources are selected by genomic unit and annotation, then metadata is
        merged within each project before combining count tables. Sample URLs use
        each project's source. Duplicate feature occurrences retain their original
        identifiers in row metadata while receiving unique output names.

        Args:
            genomic_unit: ``"gene"``, ``"exon"``, or ``"junction"``.
            annotation_extension: Optional gene/exon annotation code to select.
            join_policy: ``"inner"`` intersects features; ``"outer"`` unions
                features and inserts zeros for structural absence.
            metadata_join: ``"inner"`` intersects nonempty metadata tables within
                each project; ``"outer"`` retains all count samples.
            autoload: Whether to load missing count/metadata resources and retrieve
                coordinate sidecars needed for multi-project junction alignment.

        Returns:
            A tuple containing counts, feature metadata, sample metadata, optional
            aligned multi-project junction ranges, and experiment provenance.
            Gene/exon ranges are resolved separately by the ranged builder.

        Raises:
            ValueError: If no counts match, count values or sample mappings are
                invalid, sample names repeat, or a join policy is unsupported.
            recount3.errors.CompatibilityError: If selected organisms,
                annotations, junction formats, or feature identities conflict.
            recount3.errors.RangesError: If multi-project junctions cannot be
                aligned with their coordinate sidecars.
            recount3.errors.LoadError: If a required metadata resource cannot load.
            recount3.errors.DownloadError: If required resource retrieval fails.
        """
        unit = _utils._normalize_genomic_unit(genomic_unit)
        if metadata_join not in ("inner", "outer"):
            raise ValueError("metadata_join must be 'inner' or 'outer'.")
        if unit == "junction":
            selected = self.filter(
                resource_type="count_files_junctions", junction_extension="MM"
            )
        else:
            selected = self.filter(
                resource_type="count_files_gene_or_exon", genomic_unit=unit
            )
        if annotation_extension and unit != "junction":
            selected = selected.filter(
                annotation_extension=annotation_extension
            )
        if not selected.resources:
            raise ValueError(
                "No count-file resources available for requested annotation/unit."
            )
        feature_key = (
            "junction_type" if unit == "junction" else "annotation_extension"
        )
        signatures = {
            (
                _description_value(res, "organism"),
                _description_value(res, feature_key),
            )
            for res in selected.resources
        }
        if len(signatures) > 1:
            raise errors.CompatibilityError(
                "Ambiguous annotation, organism or junction format; "
                "select compatible counts explicitly."
            )
        annotation = next(iter(signatures))[1] if unit != "junction" else None
        keys = list(dict.fromkeys(_project_key(r) for r in selected.resources))
        frames, columns, range_frames = [], [], []
        used_resources = []
        for key in keys:
            selected_counts = [
                r for r in selected.resources if _project_key(r) == key
            ]
            related = [
                res
                for res in self.resources
                if _project_key(res) == key
                and getattr(res.description, "resource_type", None)
                in {"metadata_files", "count_files_junctions"}
                and res not in selected_counts
                and getattr(res.description, "junction_extension", None) != "MM"
            ]
            group = R3ResourceBundle(selected_counts + related)
            used_resources.extend(group.resources)
            for res in group.only_metadata().resources:
                if not res.is_loaded():
                    if not autoload:
                        raise ValueError(
                            f"Metadata resource is not loaded: {res.url}"
                        )
                    res.load()
            counts = group._stack_counts_for(
                genomic_unit=unit, join_policy=join_policy, autoload=autoload
            )
            col = group._normalize_sample_metadata(
                sample_ids=list(counts.columns),
                autoload=autoload,
                metadata_join=metadata_join,
            )
            counts = counts.loc[:, col.index]
            counts, col = _maybe_relabel_counts_columns_to_external_id(
                counts, col
            )
            col = group._add_bigwig_urls(col)
            if unit == "junction" and len(keys) > 1:
                ranges = group._junction_ranges(counts, autoload=autoload)
                counts.index = ranges.index
                range_frames.append(ranges)
            frames.append(counts)
            columns.append(col)
        counts = _merge_count_frames(frames, join_policy)
        if counts.columns.has_duplicates:
            raise ValueError(
                "Duplicate sample identifiers across selected count resources."
            )
        if counts.shape[1] == 0:
            # A recount3 experiment always has at least one sample: the
            # count resources were selected above and each carries its own
            # samples. Reaching this point means sample alignment discarded
            # all of them, which must never be returned as a valid object.
            raise ValueError(
                "Sample alignment produced no samples; the selected count "
                "resources and their metadata do not describe any sample "
                "in common."
            )
        col = pd.concat(columns, axis=0)
        col.attrs["recount3_metadata_provenance"] = {
            key: value
            for frame in columns
            for key, value in frame.attrs.get(
                "recount3_metadata_provenance", {}
            ).items()
        }
        original_ids = list(counts.index)
        if counts.index.has_duplicates:
            logging.warning(
                "Counts contain duplicate feature IDs; making row names "
                "unique and preserving feature_id."
            )
            counts.index = _make_unique_names(original_ids)
        rows = pd.DataFrame({"feature_id": original_ids}, index=counts.index)
        from recount3.version import __version__

        def collapse(values: Iterable[str | None]) -> str | list[str]:
            """Deduplicate provenance values and simplify a singleton to a string.

            Args:
                values: Ordered provenance values; ``None`` entries are omitted.

            Returns:
                The single distinct value, or a list of distinct values in their
                original order. An empty input produces an empty list.
            """
            values = list(dict.fromkeys(x for x in values if x is not None))
            return values[0] if len(values) == 1 else values

        metadata = {
            "time_created": datetime.now(timezone.utc).isoformat(),
            "recount3_version": __version__,
            "project": collapse(key[2] for key in keys),
            "organism": collapse(key[0] for key in keys),
            "project_home": collapse(
                f"data_sources/{key[1]}" for key in keys if key[1]
            ),
            "type": unit,
            "annotation": annotation,
            "resource_urls": list(
                dict.fromkeys(
                    res.url
                    for res in used_resources
                    if isinstance(res.url, str)
                )
            ),
            "metadata_columns": col.attrs["recount3_metadata_provenance"],
            "metadata_join": metadata_join,
            "join_policy": join_policy,
            "recount3_url": collapse(
                (res.config or resource.default_config()).base_url
                for res in selected.resources
                if isinstance(res, resource.R3Resource)
            ),
        }
        if unit == "junction":
            metadata["jxn_format"] = next(iter(signatures))[1]
        ranges = pd.concat(range_frames) if range_frames else None
        if ranges is not None:
            ranges = ranges.loc[~ranges.index.duplicated()].reindex(
                counts.index
            )
        return counts, rows, col, ranges, metadata

    def to_summarized_experiment(
        self,
        *,
        genomic_unit: str,
        annotation_extension: Optional[str] = None,
        assay_name: str = "raw_counts",
        join_policy: str = "inner",
        metadata_join: str = "inner",
        autoload: bool = True,
    ) -> summarizedexperiment.SummarizedExperiment:
        """Build a BiocPy :class:`SummarizedExperiment` from this bundle.

        This method stacks compatible count matrices, merges available
        sample metadata, and constructs a BiocPy
        :class:`SummarizedExperiment` using the public, validated
        :mod:`summarizedexperiment` constructor (version 0.7.1 or newer).
        Sparse count inputs remain sparse in the assay.

        Args:
          genomic_unit: Genomic unit to summarize, such as ``"gene"``,
            ``"exon"``, or ``"junction"``.
          annotation_extension: Optional annotation code for gene or
            exon summarizations (for example, ``"G026"``). When provided
            and ``genomic_unit`` is gene or exon, only count resources
            with matching annotation are used.
          assay_name: Count assay name. The default ``"raw_counts"`` becomes
            ``"counts"`` for junction experiments.
          join_policy: ``"inner"`` intersects feature rows; ``"outer"``
            unions them and fills only newly introduced rows with zero.
          metadata_join: ``"inner"`` selects the intersection of nonempty
            metadata tables within each project. ``"outer"`` retains all count samples.
          autoload: If :data:`True`, load resources when needed. If False,
            count and metadata resources must already be loaded.

        Returns:
          A BiocPy :class:`SummarizedExperiment` instance.

        Raises:
          ImportError: If BiocPy packages are not installed.
          ValueError: If counts, sample mappings, dimensions, or join policies
            are invalid, or no selected counts are available.
          recount3.errors.CompatibilityError: If selected resources or feature
            identities cannot be combined unambiguously.
          recount3.errors.RangesError: If multi-project junction sidecars cannot
            supply valid coordinate identities.
          recount3.errors.LoadError: If required metadata cannot be loaded.
          recount3.errors.DownloadError: If required resource retrieval fails.
          TypeError: If the underlying
            :class:`SummarizedExperiment` constructor rejects input types.
        """
        counts, rows, columns, _, metadata = self._prepare_experiment(
            genomic_unit=genomic_unit,
            annotation_extension=annotation_extension,
            join_policy=join_policy,
            metadata_join=metadata_join,
            autoload=autoload,
        )
        return _construct_summarized_experiment(
            counts_df=counts,
            row_df=rows,
            col_df=columns,
            assay_name=_default_assay_name(
                _utils._normalize_genomic_unit(genomic_unit), assay_name
            ),
            metadata=metadata,
        )

    def to_ranged_summarized_experiment(
        self,
        *,
        genomic_unit: str,
        annotation_extension: Optional[str] = None,
        prefer_rr_junction_coordinates: bool = True,
        assay_name: str = "raw_counts",
        join_policy: str = "inner",
        metadata_join: str = "inner",
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

        Ranges can fail to resolve for three distinct reasons, which are
        reported separately: the annotation could not be retrieved, it was
        retrieved but could not be parsed, or it parsed cleanly but does
        not describe every feature in the counts matrix. Only the last of
        those is fixed by choosing a different ``annotation_extension``.
        When ranges cannot be resolved and ``allow_fallback_to_se`` is
        :data:`True`, a plain :class:`SummarizedExperiment` is returned
        instead.

        Args:
          genomic_unit: One of ``"gene"``, ``"exon"``, or ``"junction"``.
          annotation_extension: Annotation code for gene/exon
            assays, if desired.
          prefer_rr_junction_coordinates: Whether to use RR sidecars for junction
            ranges. Junction ranges require these sidecars; disabling this
            option raises a range error or uses an explicitly enabled SE
            fallback. Ignored for gene/exon experiments.
          assay_name: Count assay name. The default ``"raw_counts"`` becomes
            ``"counts"`` for junction experiments.
          join_policy: ``"inner"`` intersects features across projects; ``"outer"``
            unions them and inserts zeros only for structurally absent features.
          metadata_join: ``"inner"`` intersects nonempty metadata tables
            within each project. ``"outer"`` retains all count
            samples. This is independent of the feature ``join_policy``.
          autoload: Whether to download and load resources as needed. With
            :data:`False`, counts and metadata must already be loaded and
            range files must be cached locally; no resources are retrieved.
          allow_fallback_to_se: If :data:`True`, return a plain
            :class:`SummarizedExperiment` instead of raising when genomic
            ranges cannot be derived. That object carries no genomic
            ranges, so anything requiring a
            :class:`RangedSummarizedExperiment` -- range queries, overlap
            operations, coordinate-based subsetting -- will not work on
            it. It only changes what is returned on failure: it does not
            retry a failed retrieval and does not repair a mismatched
            annotation. Failures in count or sample preparation, including
            multi-project junction alignment, always propagate.

        Returns:
          A :class:`RangedSummarizedExperiment` instance, or a plain
          :class:`SummarizedExperiment` when ``allow_fallback_to_se`` is
          :data:`True` and ranges are unavailable.

        Raises:
          ImportError: If BiocPy packages are not installed.
          ValueError: If counts, sample mappings, dimensions, or join policies
            are invalid, or no selected counts are available.
          recount3.errors.CompatibilityError: If selected resources or feature
            identities cannot be combined unambiguously.
          recount3.errors.RangesError: If ranges cannot be resolved and fallback
            is disabled, or multi-project junction alignment fails even with
            fallback enabled.
          recount3.errors.LoadError: If required metadata cannot be loaded.
          recount3.errors.DownloadError: If retrieval fails during count or
            metadata preparation.
          TypeError: If the
            :class:`RangedSummarizedExperiment` constructor rejects input types.
        """
        unit = _utils._normalize_genomic_unit(genomic_unit)
        counts, rows, columns, ranges, metadata = self._prepare_experiment(
            genomic_unit=unit,
            annotation_extension=annotation_extension,
            join_policy=join_policy,
            metadata_join=metadata_join,
            autoload=autoload,
        )
        try:
            if unit in {"gene", "exon"}:
                annotation = metadata["annotation"]
                annotations = self
                if isinstance(metadata["organism"], str):
                    annotations = self.filter(organism=metadata["organism"])
                selected = _select_gtf_resource_for_unit(
                    annotations,
                    genomic_unit=unit,
                    annotation_extension=annotation,
                    autoload=autoload,
                )
                if selected is None:
                    raise errors.MissingRangesError(
                        "No matching annotation resource is available."
                    )
                selected.ensure_cached(download=autoload)
                if selected.url not in metadata["resource_urls"]:
                    metadata["resource_urls"].append(selected.url)
                path = selected._cached_path()
                stat = path.stat()
                cache_key = (
                    str(path),
                    stat.st_size,
                    stat.st_mtime_ns,
                    unit,
                    tuple(rows["feature_id"]),
                )
                if (
                    self._range_cache is not None
                    and self._range_cache[0] == cache_key
                ):
                    aligned = self._range_cache[1]
                else:
                    gtf = _read_gtf_dataframe(selected)
                    features = _ranges_from_gtf(gtf, feature_kind=unit)
                    del gtf
                    aligned = _align_ranges_to_features(
                        features, feature_ids=list(rows["feature_id"])
                    )
                    aligned.index = counts.index
                    self._range_cache = (cache_key, aligned)
                if (
                    aligned[["seqnames", "starts", "ends", "strand"]]
                    .isna()
                    .any()
                    .any()
                ):
                    raise errors.RangesCoverageError(
                        "Annotation does not contain ranges for all count feature IDs."
                    )
                ranges = aligned
            elif prefer_rr_junction_coordinates:
                if ranges is None:
                    ranges = self._junction_ranges(counts, autoload=autoload)
                    original = list(rows["feature_id"])
                    counts.index = ranges.index
                    rows = pd.DataFrame(
                        {"feature_id": list(ranges.index), "mm_row": original},
                        index=ranges.index,
                    )
            else:
                raise errors.MissingRangesError(
                    "prefer_rr_junction_coordinates is False; junction "
                    "ranges require RR coordinates."
                )
            extra = ranges.drop(
                columns=["seqnames", "starts", "ends", "strand", "feature_id"],
                errors="ignore",
            )
            rows = rows.join(extra)
            return _construct_ranged_summarized_experiment(
                counts_df=counts,
                row_df=rows,
                col_df=columns,
                ranges_df=ranges,
                assay_name=_default_assay_name(unit, assay_name),
                metadata=metadata,
            )
        except (
            errors.RangesError,
            errors.DownloadError,
            errors.LoadError,
            ValueError,
            OSError,
        ) as exc:
            if not allow_fallback_to_se:
                reason = _classify_ranges_failure(
                    exc,
                    source=_RR_SOURCE if unit == "junction" else "annotation",
                )
                raise errors.RangesError(
                    f"Could not derive genomic ranges: {reason}. "
                    "Pass allow_fallback_to_se=True to receive a plain "
                    "SummarizedExperiment with no genomic ranges; this neither "
                    "retries retrieval nor repairs a mismatched annotation."
                ) from exc
            logging.warning(
                "Falling back to a plain SummarizedExperiment with no genomic ranges: %s (%s)",
                _classify_ranges_failure(
                    exc,
                    source=_RR_SOURCE if unit == "junction" else "annotation",
                ),
                exc,
            )
            metadata["ranges_error"] = str(exc)
            return _construct_summarized_experiment(
                counts_df=counts,
                row_df=rows,
                col_df=columns,
                assay_name=_default_assay_name(unit, assay_name),
                metadata=metadata,
            )

    def download(
        self,
        *,
        dest: str = ".",
        overwrite: bool = False,
        cache: r3_types.CacheMode = "enable",
        max_workers: int = 8,
    ) -> None:
        """Download all resources in the bundle to a local destination.

        This method is a convenience wrapper around
        :meth:`recount3.resource.R3Resource.download` for each contained
        resource. Downloading is I/O-bound, so by default the resources are
        fetched concurrently using a pool of worker threads (the same
        mechanism the ``recount3 download`` CLI command uses). For more
        advanced workflows (per-resource event logs, JSONL progress) prefer
        the command-line interface, ``recount3 download``.

        Concurrency is safe because per-resource downloads are coordinated by
        a shared :class:`threading.Lock` over the on-disk cache and by
        per-path locks for ``.zip`` archives, and files are materialized
        atomically. See :mod:`recount3.resource` for details.

        Args:
          dest: Destination directory or ``.zip`` path. When a directory
            is provided, each resource is materialized as a separate file
            under that directory. When a path ending in ``.zip`` is
            provided, resources are written into that archive.
          overwrite: If :data:`True`, allow overwriting existing files in
            directory mode.
          cache: Cache behavior: ``"enable"``, ``"disable"``, or
            ``"update"`` as defined by :class:`recount3.types.CacheMode`.
          max_workers: Maximum number of parallel download threads. Values
            ``<= 1`` download sequentially. The effective worker count is
            also capped at the number of resources in the bundle.

        Raises:
          ValueError: Propagated from underlying resource download
            failures, for example when an unsupported cache mode is
            selected.
        """
        resources = self.resources
        if max_workers <= 1 or len(resources) <= 1:
            for res in resources:
                res.download(path=dest, cache_mode=cache, overwrite=overwrite)
            return

        workers = min(max_workers, len(resources))
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(
                    res.download,
                    path=dest,
                    cache_mode=cache,
                    overwrite=overwrite,
                )
                for res in resources
            ]
            for fut in concurrent.futures.as_completed(futures):
                fut.result()  # re-raise the first download error, if any
