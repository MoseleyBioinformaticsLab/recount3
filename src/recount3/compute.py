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
"""Read counts, TPM, and scaling for BiocPy experiments.

The functions in this module operate on BiocPy experiments or sample metadata:

* :func:`compute_read_counts`: convert coverage-sum counts to approximate
  read counts using average mapped read length.
* :func:`compute_tpm`: compute Transcripts Per Million from an RSE with
  genomic ranges (uses ``bp_length`` when present, otherwise range widths).
* :func:`compute_scale_factors`: compute per-sample AUC- or
  mapped-reads-based scale factors.
* :func:`transform_counts`: apply scale factors to a count matrix.
* :func:`is_paired_end`: infer paired-end status from sample metadata.

Typical usage example::

  import recount3 as r3

  rse = r3.create_rse(
      project="SRP009615",
      organism="human",
      annotation_label="gencode_v26",
  )
  sf = r3.compute.compute_scale_factors(rse)  # inspect per-sample factors
  scaled = r3.compute.transform_counts(rse, by="auc")

Note:
    Functions working on BiocPy experiments need the ``biocpy`` extra.
    Install it with ``pip install "recount3[biocpy]"``.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
import functools
import logging
import numbers
import sys
from typing import Any, TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import sparse

from recount3 import _utils

if TYPE_CHECKING:  # pragma: no cover
    from numpy.typing import NDArray
    import summarizedexperiment  # type: ignore[import-not-found]


def _assay_frame(
    matrix: NDArray[Any] | sparse.spmatrix | sparse.sparray,
    row_names: Sequence[str] | None,
    column_names: Sequence[str] | None,
    *,
    copy: bool | None = None,
) -> pd.DataFrame:
    """Wrap a dense or sparse assay in a labeled DataFrame.

    Args:
        matrix: Two-dimensional assay values. SciPy sparse inputs are converted
            to pandas sparse columns without densifying the matrix; duplicate
            entries are summed, and ``matrix`` itself is never modified.
        row_names: Feature labels in assay order, or ``None`` for a RangeIndex.
        column_names: Sample labels in assay order, or ``None`` for
            a RangeIndex.
        copy: Whether a dense ``matrix`` is copied. ``None`` uses the pandas
            default (pandas 3 copies; pandas 2 does not). ``False`` wraps the
            array without copying it. Ignored for sparse input.

    Returns:
        A DataFrame preserving assay values and sparse storage where
        applicable. Dense output may share storage with the input array.

    Raises:
        ValueError: If labels do not match the matrix dimensions or the dense
            input is not two-dimensional.

    Note:
        A dense ``F`` by ``S`` assay costs ``O(F * S)`` if copied.
        With ``copy=False``, only supplied labels take time to copy,
        ``O(F + S)`` if both are provided. A sparse assay with ``N`` stored
        values takes ``O(N + F + S)`` when canonical (sorted rows, no
        duplicates), which recount3's loaders guarantee. Otherwise SciPy
        first sorts each column's rows, ``O(N log N)`` at worst, and sums
        duplicates.
    """
    kwargs = {
        "index": list(row_names) if row_names is not None else None,
        "columns": (list(column_names) if column_names is not None else None),
    }
    if sparse.issparse(matrix):
        matrix = sparse.csc_matrix(matrix)
        if not matrix.has_canonical_format:
            matrix = matrix.copy()
            matrix.sum_duplicates()
        return pd.DataFrame.sparse.from_spmatrix(matrix, **kwargs)
    return pd.DataFrame(matrix, copy=copy, **kwargs)


def _sample_metadata_values(rse: Any, column_name: str) -> NDArray[Any]:
    """Return one sample-metadata column's values.

    Resolving the column on ``col_data.to_pandas()`` converts every metadata
    column (recount3 has about two hundred) to read one. When the column data
    exposes BiocFrame's ``column_names`` and ``get_column``, the column is
    resolved with the same rules and read on its own, which gives the same
    values. Other objects, and names that only exist once nested frames are
    flattened, go through ``to_pandas()``.

    Args:
        rse: Experiment whose ``get_column_data()`` holds the column.
        column_name: Column name, resolved as by
            :func:`recount3._utils._resolve_metadata_column`.

    Returns:
        The column as a one-dimensional array.

    Raises:
        ValueError: If the column cannot be found.

    Note:
        ``O(M + S)`` for ``M`` metadata columns and ``S`` samples when the
        column is read on its own, and ``O(S * M)`` through ``to_pandas()``.
    """
    col_data = rse.get_column_data()
    column_names = getattr(col_data, "column_names", None)
    if column_names is not None and hasattr(col_data, "get_column"):
        try:
            resolved = _utils._resolve_metadata_column_name(
                column_names, column_name
            )
        except ValueError:
            resolved = None
        if resolved is not None:
            values = col_data.get_column(resolved)
            if not hasattr(values, "column_names"):
                return pd.Series(values).to_numpy()
    return _utils._resolve_metadata_column(
        col_data.to_pandas(), column_name
    ).to_numpy()


def _is_plain_biocframe(col_data: Any) -> bool:
    """Return whether ``col_data`` converts to pandas column by column.

    ``BiocFrame.to_pandas()`` builds a :class:`pandas.DataFrame` from a dict
    of the frame's columns. The fast path uses the same constructor for
    individual columns of an exact :class:`biocframe.BiocFrame` whose
    columns are plain lists or one-dimensional unstructured arrays. The
    full conversion flattens nested frames into ``parent.child`` columns,
    converts ``Factor`` columns to categoricals, and validates all arrays,
    including unused columns.

    Note:
        ``O(M)`` for ``M`` columns; no column is converted.

    Args:
        col_data: Sample metadata of an experiment.

    Returns:
        :data:`True` if every column can be read on its own.
    """
    frame_cls = getattr(sys.modules.get("biocframe"), "BiocFrame", None)
    if not isinstance(frame_cls, type) or type(col_data) is not frame_cls:
        return False
    for name in col_data.get_column_names():
        column = col_data.get_column(name)
        if type(column) is list:
            continue
        if (
            type(column) is not np.ndarray
            or column.ndim != 1
            or column.dtype.fields is not None
        ):
            return False
    return True


def _sample_metadata_columns(
    sample_metadata_source: Any,
) -> Callable[[str], pd.Series]:
    """Return a function that reads one sample-metadata column by name.

    Each column resolves as :func:`recount3._utils._resolve_metadata_column`
    resolves it on ``col_data.to_pandas()``, and reads the same values and
    errors. Converting all of ``col_data`` costs every column (recount3 has
    about two hundred), while scale factors and paired-end status need four
    or five, so plain BiocFrames are read one column at a time.

    :func:`_sample_metadata_values` serves the coverage helpers, which need
    only an array of values. The scale-factor functions return Series whose
    index and name come from these columns, so here a column is read on its
    own only when the Series provably equals the one ``to_pandas()`` holds.

    Steps

    1. A DataFrame is read directly. Any other source whose ``col_data`` is
       not plain (:func:`_is_plain_biocframe`) is converted with
       ``to_pandas()`` once, and its columns are read from that frame.
    2. For a plain BiocFrame, the labels ``to_pandas()`` would produce are
       the column names, plus the ``rownames`` column that it appends when
       the samples are named. A requested name is resolved among those
       labels with the same rules: case-insensitively, ``__`` standing in
       for the first ``.``, the last match winning.
    3. A resolved column is read through a one-column DataFrame with the
       sample names as its index. This uses the same dtype inference as
       ``to_pandas()``, including for empty columns.
    4. A name that resolves to ``rownames``, or does not resolve, is looked
       up on ``to_pandas()``, converted at most once. It supplies the sample
       names in place of any column of that name, or raises the same error
       for a missing column.

    For example, with columns ``external_id`` and ``recount_qc__star.x`` and
    named samples::

        labels:    external_id, recount_qc__star.x, rownames
        "recount_qc.star.x" -> recount_qc__star.x   (read on its own)
        "ROWNAMES"          -> rownames             (via to_pandas())
        "missing"           -> no match             (to_pandas() raises)

    Note:
        ``O(M)`` to set up for ``M`` metadata columns and ``O(M + S)`` per
        column read for ``S`` samples, where ``to_pandas()`` costs
        ``O(S * M)``. A source that is not plain costs ``O(S * M)`` once.

    Args:
        sample_metadata_source: A :class:`pandas.DataFrame`, or an object
            with ``col_data.to_pandas()``.

    Returns:
        A function from a column name to that column as a Series.

    Raises:
        TypeError: If the source is neither form.

    Examples:
        >>> column = _sample_metadata_columns(
        ...     pd.DataFrame({"External_ID": ["S1"], "recount_qc__x": [2.0]})
        ... )
        >>> column("external_id").tolist(), column("recount_qc.x").tolist()
        (['S1'], [2.0])
    """
    col_data = (
        None
        if isinstance(sample_metadata_source, pd.DataFrame)
        else getattr(sample_metadata_source, "col_data", None)
    )
    if not _is_plain_biocframe(col_data):
        metadata = _utils._coerce_col_data_to_pandas(sample_metadata_source)
        return functools.partial(_utils._resolve_metadata_column, metadata)

    row_names = col_data.get_row_names()
    labels = list(col_data.get_column_names())
    if row_names is not None:
        labels.append("rownames")
    index = None if row_names is None else pd.Index(list(row_names))

    @functools.cache
    def converted() -> pd.DataFrame:
        return col_data.to_pandas()

    def column(name: str) -> pd.Series:
        try:
            resolved = _utils._resolve_metadata_column_name(labels, name)
        except ValueError:
            resolved = None
        if resolved is None or (
            row_names is not None and resolved == "rownames"
        ):
            return _utils._resolve_metadata_column(converted(), name)
        values = col_data.get_column(resolved)
        return pd.DataFrame({resolved: values}, index=index)[resolved]

    return column


_REAL_DTYPE_KINDS = "biuf"


def _coverage_inputs(
    rse: Any,
    avg_mapped_read_length_column: str,
) -> tuple[str, Any, NDArray[np.float64]]:
    """Resolve and validate a coverage assay and per-sample read lengths.

    Shared by :func:`compute_read_counts` and :func:`compute_tpm`.

    Args:
        rse: Experiment holding a coverage-sum assay and sample metadata.
        avg_mapped_read_length_column: Metadata column with each sample's
            average mapped read length.

    Returns:
        The resolved assay name; the assay values, as the stored SciPy matrix
        for a sparse assay, or as an array for a dense one (boolean, integer,
        and float arrays keep their dtype; anything else is converted to
        float); and the read lengths as a float array (missing values are
        NaN).

    Raises:
        ValueError: If the assay or the read-length column is missing, a read
            length is not numeric, the assay is not two-dimensional, or its
            sample count differs from the number of read lengths.

    Note:
        ``O(M + S)`` for ``M`` metadata columns and ``S`` samples (see
        :func:`_sample_metadata_values`). The assay is returned as stored,
        except that a dense assay of a non-numeric dtype is converted to
        float in ``O(F * S)``.
    """
    assay_name = _utils._resolve_counts_assay_name(rse)

    try:
        avg_len_values = _sample_metadata_values(
            rse, avg_mapped_read_length_column
        )
    except ValueError as exc:
        raise ValueError(
            "Required metadata column not found in col_data: "
            f"{avg_mapped_read_length_column!r}."
        ) from exc

    assay = rse.get_assay(assay_name)
    if sparse.issparse(assay):
        counts = assay
    else:
        counts = np.asarray(assay)
        if counts.dtype.kind not in _REAL_DTYPE_KINDS:
            counts = np.asarray(assay, dtype=float)
    if counts.ndim != 2:
        raise ValueError(
            f"{assay_name!r} assay must be a 2D matrix "
            f"(got shape {counts.shape})."
        )

    for value in avg_len_values:
        if pd.isna(value):
            continue
        if isinstance(value, (str, bytes)):
            raise ValueError(
                f"Metadata column {avg_mapped_read_length_column!r} "
                f"must be numeric; found {type(value).__name__} "
                f"value {value!r}."
            )
        if not isinstance(value, numbers.Real):
            raise ValueError(
                f"Metadata column {avg_mapped_read_length_column!r} "
                f"must be numeric; found {type(value).__name__}."
            )

    avg_len = np.asarray(avg_len_values, dtype=float)

    if counts.shape[1] != avg_len.shape[0]:
        raise ValueError(
            f"Mismatch between number of samples in {assay_name!r} "
            f"({counts.shape[1]}) and length of "
            f"{avg_mapped_read_length_column!r} ({avg_len.shape[0]})."
        )
    return assay_name, counts, avg_len


def compute_read_counts(
    rse: summarizedexperiment.RangedSummarizedExperiment,
    round_to_integers: bool = True,
    avg_mapped_read_length_column: str = (
        "recount_qc.star.average_mapped_length"
    ),
) -> pd.DataFrame:
    """Convert coverage-sum counts into approximate read/fragment counts.

    The "raw_counts" assay used by recount3-style resources represents summed
    per-base coverage over each feature (for example, total coverage across all
    bases in a gene), not the number of reads/fragments overlapping the feature.
    A common approximation for converting coverage-sum values into read/fragment
    counts is to divide by the sample's average mapped read length.

    For each feature i and sample j:

      read_counts[i, j] = raw_counts[i, j] / avg_mapped_read_length_column[j]

    This operation is performed column-wise: each sample column in the assay is
    divided by a single scalar derived from that sample's metadata.

    Rounding is optional. When enabled, values are rounded to 0 decimals to
    produce integer-like counts, which is convenient for downstream methods that
    assume counts are integer-valued.

    Args:
      rse: A :class:`~summarizedexperiment.RangedSummarizedExperiment`-like
        object containing a "raw_counts" assay and sample metadata in
        `col_data`.
      round_to_integers: If True, round the resulting values to 0 decimals.
      avg_mapped_read_length_column: Name of the metadata column containing
        average mapped read length per sample. The column must be present in
        `col_data` and must contain numeric values.

    Returns:
      A :class:`~pandas.DataFrame` of approximate read counts with the same
      shape as the "raw_counts" assay (features x samples). Row and column
      names are preserved when available. Sparse inputs produce a sparse
      DataFrame, including missing-value propagation.

    Raises:
      TypeError: If `rse` is not a
        :class:`~summarizedexperiment.RangedSummarizedExperiment`, or if
        `round_to_integers` is not a bool.
      ValueError: If the "raw_counts" assay is missing, if
        `avg_mapped_read_length_column` is missing from `col_data`, or if the
        assay and metadata dimensions do not align.

    Note:
        ``O(F * S + M)`` for a dense assay of ``F`` features and ``S``
        samples with ``M`` metadata columns. A sparse assay with ``N`` stored
        values takes ``O(N + F + S + M)``: division and rounding touch only
        the stored values, and the result stays sparse.
    """
    ranged_summarized_experiment_cls = (
        _utils.get_ranged_summarizedexperiment_class()
    )

    if not isinstance(rse, ranged_summarized_experiment_cls):
        raise TypeError(
            "rse must be a RangedSummarizedExperiment to compute read counts."
        )
    if not isinstance(round_to_integers, bool):
        raise TypeError("round_to_integers must be a bool.")

    _, counts, avg_len = _coverage_inputs(rse, avg_mapped_read_length_column)
    raw_counts = (
        counts.astype(float)
        if sparse.issparse(counts)
        else np.asarray(counts, dtype=float)
    )

    read_counts = _assay_frame(
        raw_counts, rse.get_row_names(), rse.get_column_names()
    ).div(avg_len, axis=1)
    return read_counts.apply(np.rint) if round_to_integers else read_counts


def _feature_lengths(rse: Any) -> NDArray[np.float64]:
    """Return ``bp_length`` from row data when present, else range widths.

    Note:
        ``O(F)`` for ``F`` features.

    Raises:
        ValueError: If neither is available.
    """
    try:
        row_data = rse.get_row_data() if hasattr(rse, "get_row_data") else None
        if row_data is not None and "bp_length" in row_data.column_names:
            return np.asarray(row_data["bp_length"], dtype=float)
        return np.asarray(rse.width, dtype=float)
    except AttributeError as exc:
        raise ValueError(
            "Could not determine feature widths from rowRanges."
        ) from exc


def _dense_tpm(
    counts: NDArray[Any],
    lengths: NDArray[np.float64],
    usable_samples: NDArray[np.bool_],
) -> NDArray[np.float64]:
    """Return TPM values for a dense assay as a new float array.

    The result is column-major, the layout pandas uses for DataFrame columns,
    so each sample is contiguous. The division converts integer counts as it
    writes, so the result is the only matrix-sized allocation. Callers
    silence floating-point warnings; see :func:`compute_tpm`.

    Steps

    1. Divide each count by its feature's length into the output array.
    2. Sum each sample. Only if a sum is NaN, sum again skipping missing
       values, which needs a boolean mask the size of the assay.
    3. Multiply each sample by ``1e6 / sum``, or by NaN where the sample
       cannot be normalized.

    Note:
        ``O(F * S)`` time for ``F`` features and ``S`` samples, with one
        float64 output array and, only when values are missing, one boolean
        mask of the same shape.
    """
    tpm = np.empty(counts.shape, order="F")
    np.divide(counts, lengths[:, np.newaxis], out=tpm)
    sample_sums = tpm.sum(axis=0)
    if np.isnan(sample_sums).any():
        observed = np.isnan(tpm)
        np.logical_not(observed, out=observed)
        sample_sums = np.sum(tpm, axis=0, where=observed)
    tpm *= np.where(usable_samples, 1_000_000.0 / sample_sums, np.nan)
    return tpm


def _sparse_tpm(
    counts: sparse.spmatrix | sparse.sparray,
    lengths: NDArray[np.float64],
    usable_samples: NDArray[np.bool_],
) -> sparse.csc_matrix:
    """Return TPM values for a sparse assay as a new CSC matrix.

    Only stored values are scaled. An implicit zero becomes
    ``0 / length * scale``, which is 0 unless the length is missing or zero
    or the scale is not finite; those entries are stored as NaN. Callers
    silence floating-point warnings; see :func:`compute_tpm`.

    Steps

    1. Copy the assay to a float CSC matrix with sorted rows and no
       duplicate entries, and record the sample of each stored value.
    2. Divide each stored value by its feature's length.
    3. Sum each sample's stored values, skipping missing ones, with one
       ``bincount`` over the values' samples.
    4. Multiply each stored value by its sample's ``1e6 / sum``, or by NaN
       where the sample cannot be normalized.
    5. If any feature has a missing or zero length, or any sample cannot be
       normalized, store NaN at those lines' implicit zeros with
       :func:`_store_nan_lines`.

    Note:
        ``O(N + F + S + Z)`` for ``N`` stored values, ``F`` features, ``S``
        samples and ``Z`` stored NaN entries. The assay is never densified.
        An assay whose rows are not sorted within each column is sorted
        first, ``O(N log N)`` at worst.
    """
    tpm = sparse.csc_matrix(counts, dtype=float, copy=True)
    tpm.sum_duplicates()
    n_samples = tpm.shape[1]
    sample_of_value = np.repeat(np.arange(n_samples), np.diff(tpm.indptr))

    tpm.data /= lengths[tpm.indices]
    sample_sums = np.bincount(
        sample_of_value,
        weights=np.where(np.isnan(tpm.data), 0.0, tpm.data),
        minlength=n_samples,
    )
    scale = np.where(usable_samples, 1_000_000.0 / sample_sums, np.nan)
    tpm.data *= scale[sample_of_value]

    nan_features = np.isnan(lengths) | (lengths == 0)
    nan_samples = ~np.isfinite(scale)
    if nan_features.any() or nan_samples.any():
        tpm = _store_nan_lines(tpm, sample_of_value, nan_features, nan_samples)
    return tpm


def _store_nan_lines(
    matrix: sparse.csc_matrix,
    sample_of_value: NDArray[np.intp],
    nan_features: NDArray[np.bool_],
    nan_samples: NDArray[np.bool_],
) -> sparse.csc_matrix:
    """Store NaN at every implicit zero of the flagged features and samples.

    The NaN entries are laid out directly in sorted CSC order, skipping
    positions that already hold a stored value, and merged with ``matrix``,
    in time linear in the size of the result.

    1. Lay out the NaN entries sample by sample, in CSC order: a flagged
       sample gets every feature, and every other sample gets the flagged
       features in ascending order.
    2. Find the stored values that fall on a flagged feature or sample.
       Their slots in the layout are already filled, so they are marked:
       a value's slot is its feature in a flagged sample's block, and its
       feature's rank among the flagged features in any other block.
    3. Drop the marked slots and shift each sample's start by the number
       of slots dropped before it.
    4. Add the NaN matrix to ``matrix``. The two store disjoint positions,
       so the sum keeps every stored value and adds the NaN entries.

    Note:
        ``O(N + Z + F + S)`` for ``N`` stored values, ``Z`` NaN entries,
        ``F`` features and ``S`` samples: each step is one vectorized pass,
        and SciPy merges two canonical matrices in linear time.

    Examples:
        Feature 1 has no usable length and sample 1 cannot be normalized,
        so every implicit zero on either becomes NaN, while stored values,
        such as the 2 of feature 2 in sample 1, are kept:

        >>> matrix = sparse.csc_matrix(
        ...     np.array([[5.0, 0.0], [0.0, 0.0], [0.0, 2.0]])
        ... )
        >>> _store_nan_lines(
        ...     matrix,
        ...     sample_of_value=np.array([0, 1]),
        ...     nan_features=np.array([False, True, False]),
        ...     nan_samples=np.array([False, True]),
        ... ).toarray()
        array([[ 5., nan],
               [nan, nan],
               [ 0.,  2.]])

    Args:
        matrix: Canonical CSC matrix. It is not modified.
        sample_of_value: Sample (column) index of each stored value.
        nan_features: Features whose implicit zeros become NaN.
        nan_samples: Samples whose implicit zeros become NaN.

    Returns:
        A canonical CSC matrix with the stored values of ``matrix`` and NaN at
        the implicit zeros of the flagged lines.
    """
    n_features, n_samples = matrix.shape
    flagged_features = np.flatnonzero(nan_features)
    entries_per_sample = np.where(
        nan_samples, n_features, flagged_features.size
    )
    indptr = np.concatenate(([0], np.cumsum(entries_per_sample)))
    indices = np.empty(indptr[-1], dtype=np.int64)
    full = indptr[:-1][nan_samples]
    indices[full[:, np.newaxis] + np.arange(n_features)] = np.arange(n_features)
    partial = indptr[:-1][~nan_samples]
    indices[partial[:, np.newaxis] + np.arange(flagged_features.size)] = (
        flagged_features
    )

    feature_of_value = matrix.indices
    in_lines = nan_features[feature_of_value] | nan_samples[sample_of_value]
    features = feature_of_value[in_lines]
    samples = sample_of_value[in_lines]
    rank = np.zeros(n_features, dtype=np.int64)
    rank[flagged_features] = np.arange(flagged_features.size)
    offsets = np.where(nan_samples[samples], features, rank[features])
    occupied = np.zeros(indptr[-1], dtype=bool)
    occupied[indptr[samples] + offsets] = True
    indptr -= np.concatenate(
        ([0], np.cumsum(np.bincount(samples, minlength=n_samples)))
    )
    nan_lines = sparse.csc_matrix(
        (np.full(indptr[-1], np.nan), indices[~occupied], indptr),
        shape=matrix.shape,
    )
    return nan_lines + matrix


def compute_tpm(
    rse: summarizedexperiment.RangedSummarizedExperiment,
) -> pd.DataFrame:
    """Compute Transcripts Per Million (TPM) from raw coverage sums.

    TPM is defined as:
      1. Approximate Read Counts = Coverage AUC / Avg Read Length
      2. RPK = Read Counts / (Feature Length / 1000)
      3. Scale Factor = Sum(RPK) / 1,000,000
      4. TPM = RPK / Scale Factor

    The read length and the factor of 1000 are the same for every feature in
    a sample, so they cancel in step 4, and each value is computed directly
    as ``1e6 * (coverage / length) / sum(coverage / length)`` over the
    sample. The read length (``recount_qc.star.average_mapped_length``) still
    decides whether a sample can be normalized: a missing, zero, or infinite
    value makes the whole sample NaN. Sums skip missing values, so a feature
    with a missing length is NaN while the other features in each sample
    still total one million. A zero length is NaN too; where that feature has
    coverage, its RPK is infinite and the sample's other values become 0. A
    sample with no coverage is NaN.

    Note:
        A dense assay with ``F`` features and ``S`` samples takes ``O(F * S)``
        time, with one output-sized array and a boolean mask of the same
        shape when needed to skip NaNs. A canonical sparse assay takes
        ``O(N + F + S + Z)`` for ``N`` stored input values and ``Z`` implicit
        zeros written as NaNs. ``Z`` can reach ``F * S`` if every sample
        has zero coverage. Sorting unsorted indices can add ``O(N log N)``
        time. Reading the read-length column directly takes ``O(M + S)``
        for ``M`` metadata columns.

    Args:
      rse: A :class:`~summarizedexperiment.RangedSummarizedExperiment`
        object containing raw coverage sums. Uses annotated ``bp_length``
        when available, otherwise the genomic widths in ``row_ranges``.

    Returns:
      A :class:`~pandas.DataFrame` of float64 TPM values in feature-by-sample
      order, preserving row and column names. A sparse count assay gives
      sparse columns with a fill value of 0 that store at most the assay's
      stored values plus the NaN entries. The assay itself is not modified.

    Raises:
      TypeError: If rse is not a
        :class:`~summarizedexperiment.RangedSummarizedExperiment`
        (needs rowRanges).
      ValueError: If feature lengths cannot be obtained or their number
        differs from the number of features, the count assay or read-length
        metadata column is absent, a read length is not numeric, or the
        dimensions do not align.

    Examples:
        Compute TPM from an RSE built with :func:`create_rse`::

            rse = r3.create_rse(project="SRP009615", organism="human",
                             annotation_label="gencode_v26")
            tpm_df = r3.compute.compute_tpm(rse)
    """
    ranged_summarized_experiment_cls = (
        _utils.get_ranged_summarizedexperiment_class()
    )

    if not isinstance(rse, ranged_summarized_experiment_cls):
        raise TypeError(
            "compute_tpm requires a RangedSummarizedExperiment with "
            "defined rowRanges to determine feature lengths."
        )

    assay_name, counts, read_lengths = _coverage_inputs(
        rse, "recount_qc.star.average_mapped_length"
    )
    lengths = _feature_lengths(rse)
    if lengths.shape != (counts.shape[0],):
        raise ValueError(
            f"Expected one feature length per row of the {assay_name!r} "
            f"assay ({counts.shape[0]}); got {lengths.size}."
        )
    usable_samples = np.isfinite(read_lengths) & (read_lengths != 0)

    with np.errstate(divide="ignore", invalid="ignore"):
        if sparse.issparse(counts):
            tpm = _sparse_tpm(counts, lengths, usable_samples)
        else:
            tpm = _dense_tpm(counts, lengths, usable_samples)
    return _assay_frame(
        tpm, rse.get_row_names(), rse.get_column_names(), copy=False
    )


_AVG_READ_LENGTH_COLUMN = "recount_seq_qc.avg_len"


def is_paired_end(
    sample_metadata_source: Any,
    avg_mapped_read_length_column: str = (
        "recount_qc.star.average_mapped_length"
    ),
    avg_read_length_column: str = _AVG_READ_LENGTH_COLUMN,
) -> pd.Series:
    """Infer paired-end status, matching recount3::is_paired_end().

    In recount3 (R), paired-end status is inferred via:
      ratio <- round(avg_mapped_read_length / avg_read_length, 0)
      ratio must be 1 (single-end) or 2 (paired-end), otherwise NA with warning.
      result <- ratio == 2, with names(result) = external_id.

    Args:
        sample_metadata_source: Sample metadata (DataFrame) or a
          ``(Ranged)SummarizedExperiment``-like object with
          `col_data.to_pandas()`.
        avg_mapped_read_length_column: Metadata column with average mapped
          length.
        avg_read_length_column: Metadata column containing average read length.

    Returns:
        A pandas Series of dtype "boolean" (True/False/pd.NA), indexed by
        external_id.

    Raises:
        ValueError: If required metadata columns are missing or non-numeric.

    Note:
        ``O(S + M)`` for ``S`` samples and ``M`` metadata columns: the three
        columns used are found among the ``M`` names and read on their own,
        also from a BiocPy experiment whose column data holds plain lists
        and one-dimensional unstructured arrays. Other experiment-like
        objects are converted with ``col_data.to_pandas()``, ``O(S * M)``.

    Examples:
        >>> metadata = pd.DataFrame(
        ...     {
        ...         "external_id": ["SRR1", "SRR2", "SRR3"],
        ...         "recount_qc.star.average_mapped_length": [99.0, 202.0, 50],
        ...         "recount_seq_qc.avg_len": [100.0, 101.0, 100.0],
        ...     }
        ... )
        >>> is_paired_end(metadata).tolist()
        [False, True, <NA>]
    """
    return _paired_end_status(
        _sample_metadata_columns(sample_metadata_source),
        avg_mapped_read_length_column,
        avg_read_length_column,
    )


def _paired_end_status(
    column: Callable[[str], pd.Series],
    avg_mapped_read_length_column: str,
    avg_read_length_column: str,
) -> pd.Series:
    """Infer paired-end status from sample-metadata columns.

    Steps

    1. Read ``external_id`` as text, and both lengths as numbers.
    2. Round each sample's mapped-to-read length ratio half to even, as R's
       ``round(x, 0)`` does.
    3. A ratio of 1 is single-end, 2 is paired-end; any other ratio is
       unknown and is logged once for all such samples.

    Note:
        ``O(S)`` for ``S`` samples, plus the cost of reading three columns.

    Args:
        column: Reads a metadata column by name; see
            :func:`_sample_metadata_columns`.
        avg_mapped_read_length_column: Column with average mapped length.
        avg_read_length_column: Column with average read length.

    Returns:
        The result documented by :func:`is_paired_end`.

    Raises:
        ValueError: If required metadata columns are missing or non-numeric.
    """
    external_id = column("external_id").map(str)
    avg_mapped = _utils._coerce_numeric_column(
        column(avg_mapped_read_length_column), avg_mapped_read_length_column
    )
    avg_len = _utils._coerce_numeric_column(
        column(avg_read_length_column), avg_read_length_column
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


def compute_scale_factors(  # noqa: PLR0917 - Existing positional API.
    sample_metadata_source: Any,
    by: str = "auc",
    target_read_count: float = 4e7,
    target_read_length_bp: float = 100,
    auc_column: str = "recount_qc.bc_auc.all_reads_all_bases",
    avg_mapped_read_length_column: str = (
        "recount_qc.star.average_mapped_length"
    ),
    mapped_reads_column: str = "recount_qc.star.all_mapped_reads",
    paired_end_status: Sequence[bool] | pd.Series | None = None,
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
       common target_read_count:

         s[j] = target_read_count / auc[j]

       This method preserves relative feature coverage within each sample while
       adjusting overall sample magnitude to be comparable across samples.

    2) by="mapped_reads"
       Uses mapped read counts and read length to normalize samples to a common
       target_read_count and a common target read length target_read_length_bp:

         s[j] = target_read_count * target_read_length_bp *
                paired_multiplier[j] /
                (mapped_reads[j] * (avg_mapped_read_length[j] ** 2))

       paired_multiplier is:
         - 2 for paired-end samples
         - 1 for single-end samples
         - missing for samples whose paired-end status cannot be inferred

       If `paired_end_status` is not provided, paired-end status is inferred
       from metadata by comparing average mapped length to average read
       length:

         ratio = round(avg_mapped_read_length / avg_read_length)

       ratio==2 indicates paired-end, ratio==1 indicates single-end. Other
       ratios are treated as unknown and produce missing paired multipliers.

    Paired-end status is resolved only for by="mapped_reads". It cancels out
    of the AUC formula, so by="auc" never reads `avg_read_length_column` and
    keeps working on projects whose ``recount_seq_qc`` table is empty or was
    not published. This mirrors recount3's R implementation, where the
    ``paired_end`` default argument is an ``is_paired_end()`` call that the
    AUC branch never forces.

    Missing values in required metadata propagate to missing scale factors.
    Non-numeric metadata values raise an error.

    Args:
      sample_metadata_source: Sample metadata as a
        :class:`~pandas.DataFrame`, or an object with
        `col_data.to_pandas()` that yields sample metadata.
      by: Scaling method: "auc" or "mapped_reads".
      target_read_count: Target library size used to compute scale factors.
        Interpreted as the number of single-end reads to scale each sample to.
      target_read_length_bp: Target read length used only when
        by="mapped_reads".
      auc_column: Metadata column name for the per-sample AUC metric.
      avg_mapped_read_length_column: Metadata column name for average mapped
        read length per sample.
      mapped_reads_column: Metadata column name for mapped read counts per
        sample.
      paired_end_status: Optional paired-end indicator per sample. If
        provided, it must align with the samples in `external_id`. If
        omitted, paired-end status is inferred from metadata. Used only
        when by="mapped_reads"; ignored for by="auc", which does not
        depend on it.

    Returns:
      A pandas Series of scale factors indexed by `external_id`. The Series name
      is "scale_factor".

    Raises:
      ValueError: If `by` is invalid, required metadata columns are missing, or
        non-numeric metadata values are present. `auc_column`,
        `avg_mapped_read_length_column`, `mapped_reads_column` and
        `external_id` are required for both methods;
        `avg_read_length_column` is required only by by="mapped_reads",
        and only when `paired_end_status` is not supplied.
      TypeError: If `target_read_count` or `target_read_length_bp` are not
        numeric scalars.

    Note:
        ``O(S + M)`` for ``S`` samples and ``M`` metadata columns: the four
        or five columns used are found among the ``M`` names and read on
        their own, and every factor is one vectorized operation per sample.
        A BiocPy experiment whose column data holds plain lists and
        one-dimensional unstructured arrays is read the same way. Other
        experiment-like objects are converted with ``col_data.to_pandas()``,
        ``O(S * M)``.

    Examples:
        AUC-based scaling (default)::

            sf = r3.compute.compute_scale_factors(rse)  # inspect factors
            scaled = r3.compute.transform_counts(rse, by="auc")  # apply scaling

        Mapped-reads-based scaling::

            sf = r3.compute.compute_scale_factors(rse, by="mapped_reads")
    """
    if by not in ("auc", "mapped_reads"):
        raise ValueError(
            f"'by' must be either 'auc' or 'mapped_reads'; got {by!r}."
        )

    if not isinstance(target_read_count, numbers.Real) or isinstance(
        target_read_count, bool
    ):
        raise TypeError("target_read_count must be a numeric scalar.")
    if not isinstance(target_read_length_bp, numbers.Real) or isinstance(
        target_read_length_bp, bool
    ):
        raise TypeError("target_read_length_bp must be a numeric scalar.")

    column = _sample_metadata_columns(sample_metadata_source)

    # Match recount3's stopifnot(): all of these must exist even if by="auc".
    external_id = column("external_id").map(str)

    auc_values = _utils._coerce_numeric_column(column(auc_column), auc_column)
    avg_mapped_values = _utils._coerce_numeric_column(
        column(avg_mapped_read_length_column), avg_mapped_read_length_column
    )
    mapped_reads_values = _utils._coerce_numeric_column(
        column(mapped_reads_column), mapped_reads_column
    )

    auc_values.index = external_id
    avg_mapped_values.index = external_id
    mapped_reads_values.index = external_id

    if by == "auc":
        scale_factor = float(target_read_count) / auc_values
    else:
        if paired_end_status is None:
            paired_end_series = _paired_end_status(
                column, avg_mapped_read_length_column, _AVG_READ_LENGTH_COLUMN
            )
        else:
            paired_end_series = pd.Series(
                paired_end_status,
                index=external_id,
                dtype="boolean",
            )

        pe_multiplier = pd.Series(np.nan, index=external_id, dtype=float)
        # Compare elementwise; the nullable Series also contains unknowns.
        pe_multiplier.loc[paired_end_series.eq(other=True)] = 2.0
        pe_multiplier.loc[paired_end_series.eq(other=False)] = 1.0

        denom = mapped_reads_values * (avg_mapped_values**2)
        scale_factor = (
            float(target_read_count)
            * float(target_read_length_bp)
            * pe_multiplier
            / denom
        )

    scale_factor.name = "scale_factor"
    return scale_factor


def transform_counts(
    rse: summarizedexperiment.RangedSummarizedExperiment,
    by: str = "auc",
    target_read_count: float = 4e7,
    target_read_length_bp: float = 100,
    round_to_integers: bool = True,
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
        scale_factor[j] = target_read_count / auc[j]

      where auc is a per-sample total coverage metric. This method scales each
      sample to have total coverage approximately equal to target_read_count.

    - by="mapped_reads":
        scale_factor[j] = (
            target_read_count * target_read_length_bp * paired_multiplier[j]
            / (mapped_reads[j] * (avg_mapped_read_length[j] ** 2))

        )

      where paired_multiplier is 2 for paired-end samples, 1 for single-end
      samples, and missing when paired-end status cannot be inferred. This
      method incorporates mapped reads and read length so that samples with
      different read lengths are normalized onto the same target read length
      target_read_length_bp.

    The returned values remain in the same feature-by-sample shape as the
    input. If round_to_integers is True, values are rounded to integer-like
    counts.

    Args:
      rse: A :class:`~summarizedexperiment.RangedSummarizedExperiment`-like
        object containing a "raw_counts" assay and sample metadata in
        `col_data`.
      by: Scaling method: "auc" or "mapped_reads".
      target_read_count: Target library size used to compute scale factors.
        Interpreted as the number of single-end reads to scale each sample to.
      target_read_length_bp: Target read length used only when
        by="mapped_reads".
      round_to_integers: If True, round scaled values to 0 decimals.
      **kwargs: Additional parameters forwarded to
        :func:`~recount3.compute.compute_scale_factors`. Use this to override
        metadata column names (for example, `auc_column=...`,
        `mapped_reads_column=...`, `avg_mapped_read_length_column=...`) or to
        provide `paired_end_status=...` when paired-end status should not be
        inferred.

    Returns:
      A :class:`~pandas.DataFrame` of scaled counts with the same dimensions
      as `assay("raw_counts")`. Row and column names are preserved when
      available. Sparse assays produce sparse DataFrames; missing scale
      factors propagate to every value in the corresponding sample column.

    Raises:
      ValueError: If the required assay or metadata columns are missing,
        ``by`` is invalid, or assay and metadata dimensions do not align.
      TypeError: If ``rse`` is not a RangedSummarizedExperiment,
        ``round_to_integers`` is not a bool, or numeric parameters are not
        valid scalars.

    Note:
        ``O(F * S + M)`` for a dense assay of ``F`` features and ``S``
        samples with ``M`` metadata columns, and ``O(N + F + S + M)`` for a
        sparse assay with ``N`` stored values, which stays sparse: scaling
        and rounding touch only stored values. The scale factors cost
        ``O(S + M)`` (see :func:`compute_scale_factors`).
    """
    ranged_summarized_experiment_cls = (
        _utils.get_ranged_summarizedexperiment_class()
    )

    if not isinstance(rse, ranged_summarized_experiment_cls):
        raise TypeError(
            "rse must be a RangedSummarizedExperiment (BiocPy "
            "summarizedexperiment)."
        )

    assay_name = _utils._resolve_counts_assay_name(rse)

    if not isinstance(round_to_integers, bool):
        raise TypeError("round_to_integers must be a bool.")

    counts = rse.get_assay(assay_name)
    counts_array = (
        counts.astype(float)
        if sparse.issparse(counts)
        else np.asarray(counts, dtype=float)
    )

    scale_factor = compute_scale_factors(
        rse,
        by=by,
        target_read_count=target_read_count,
        target_read_length_bp=target_read_length_bp,
        **kwargs,
    )

    if counts_array.ndim != 2:
        raise ValueError(
            f"{assay_name!r} assay must be 2D (got shape {counts_array.shape})."
        )
    if counts_array.shape[1] != len(scale_factor):
        raise ValueError(
            f"Mismatch between number of samples in {assay_name!r} "
            f"({counts_array.shape[1]}) and number of scale factors "
            f"({len(scale_factor)})."
        )

    scaled = _assay_frame(
        counts_array, rse.get_row_names(), rse.get_column_names()
    ).mul(scale_factor.to_numpy(dtype=float), axis=1)
    return scaled.apply(np.rint) if round_to_integers else scaled
