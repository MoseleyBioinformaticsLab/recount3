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
"""High-level builders and utilities for SummarizedExperiment objects.

This module provides helpers for constructing BiocPy
:class:`~summarizedexperiment.SummarizedExperiment` and
:class:`~summarizedexperiment.RangedSummarizedExperiment` objects from
recount3 projects, plus utilities for working with SRA-style sample
attributes. See :mod:`recount3.compute` for normalization and scaling.

The heavy lifting is implemented as methods on
:class:`~recount3.bundle.R3ResourceBundle`:

* :meth:`~recount3.bundle.R3ResourceBundle.to_summarized_experiment`
* :meth:`~recount3.bundle.R3ResourceBundle.to_ranged_summarized_experiment`

The functions in this module are thin wrappers around those methods (for
example, :func:`create_rse`) and convenience utilities that operate on
BiocPy objects directly:

* :func:`expand_sra_attributes`: parse SRA ``key;;value|...`` attribute
  strings into separate columns on a metadata DataFrame or SE/RSE object.
* :func:`to_anndata`: convert an experiment to an ``anndata.AnnData``
  with samples in rows, for Python machine-learning tooling.

Note:
    Most functions in this module require BiocPy packages
    (``biocframe``, ``summarizedexperiment``, ``genomicranges``).
    Install them with::

        pip install "recount3[biocpy]"
"""

from __future__ import annotations

from collections.abc import Sequence
import logging
from typing import Any, TYPE_CHECKING
import warnings

import pandas as pd

from recount3 import _utils
from recount3 import bundle as r3_bundle
from recount3 import compute
from recount3 import search

if TYPE_CHECKING:  # pragma: no cover
    import summarizedexperiment  # type: ignore[import-not-found]


def _resolve_sra_attributes_column(
    columns: pd.Index | Sequence[str],
    requested: str,
) -> str | None:
    """Pick which column actually holds the SRA attribute strings.

    R-style metadata frames spell this column ``'sra.sample_attributes'``,
    while the Python bundle layer namespaces metadata columns as
    ``'{table}__{column}'`` and produces ``'sra__sample_attributes'`` for
    the same field. Callers may pass either form; this helper tries the
    requested name first, then the ``.`` <-> ``__`` variant.

    Args:
      columns: Available column names of the metadata frame.
      requested: The user-supplied (or default) column name.

    Returns:
      The matching name from ``columns`` if either form is present, else
      :data:`None`.
    """
    cols = set(columns)
    if requested in cols:
        return requested
    if "." in requested:
        candidate = requested.replace(".", "__", 1)
        if candidate in cols:
            return candidate
    if "__" in requested:
        candidate = requested.replace("__", ".", 1)
        if candidate in cols:
            return candidate
    return None


def _expand_sra_attributes_df(
    col_df: pd.DataFrame,
    *,
    sra_attributes_column: str = "sra.sample_attributes",
    attribute_column_prefix: str = "sra_attribute.",
) -> pd.DataFrame:
    """Return a copy of column metadata with SRA attributes expanded.

    recount3 stores SRA sample attributes in a single string column such
    as ``'sra.sample_attributes'`` using the encoding::

        "age;;67.78|biomaterial_provider;;LIBD|disease;;Control|..."

    This helper mirrors the behavior of the recount3 R function
    ``expand_sra_attributes()``: each entry is split on ``'|'`` into
    individual attributes, then on ``';;'`` into ``key`` and ``value``
    pairs. New columns named ``{attribute_column_prefix}{key}`` (with
    spaces in ``key`` replaced by ``'_'``) are appended to the metadata.

    Args:
      col_df: Column metadata DataFrame.
      sra_attributes_column: Name of the column carrying the raw SRA
        attribute strings. If this exact name is absent, the helper tries the
        ``.`` <-> ``__`` variant so the same call works against both
        R-style (``sra.sample_attributes``) and Python-namespaced
        (``sra__sample_attributes``) metadata. If neither form is
        present, ``col_df`` is returned unchanged.
      attribute_column_prefix: Prefix to prepend to attribute names when
        forming new columns.

    Returns:
      A new DataFrame with the same index as ``col_df`` and additional
      columns containing one parsed SRA attribute per column.

    Note:
        ``O(L + S * (M + K))`` for ``S`` samples, ``L`` characters of
        attribute text, ``M`` existing columns and ``K`` distinct attribute
        keys: each string is split once, and the result is an ``S`` by
        ``M + K`` table. That is the size of the input plus the output.

    Examples:
        >>> metadata = pd.DataFrame(
        ...     {"sra.sample_attributes": ["age;;67|tissue;;brain", None]}
        ... )
        >>> _expand_sra_attributes_df(metadata).iloc[:, 1:]
          sra_attribute.age sra_attribute.tissue
        0                67                brain
        1               NaN                  NaN
    """
    resolved = _resolve_sra_attributes_column(
        col_df.columns, sra_attributes_column
    )
    if resolved is None:
        return col_df.copy()

    ser = col_df[resolved]
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
                key_slug = key_clean.replace(" ", "_")
                col_name = f"{attribute_column_prefix}{key_slug}"
                row_map[col_name] = val.strip()
        rows.append(row_map)

    if not rows:
        return col_df.copy()

    attrs_df = pd.DataFrame(rows, index=col_df.index)
    return pd.concat([col_df, attrs_df], axis=1)


def _resolve_annotation_extension(
    *,
    organism: str,
    genomic_unit: str,
    annotation_label: str | None,
    annotation_extension: str | None,
) -> str | None:
    """Resolve an annotation extension for :func:`~recount3.create_rse`.

    For gene and exon assays, this helper maps a human-readable
    annotation name (for example, ``"gencode_v26"``) or an explicit
    extension (for example, ``"G026"``) to the underlying recount3
    extension. When both ``annotation_label`` and
    ``annotation_extension`` are provided, the explicit extension
    wins.

    Args:
      organism: Organism identifier (``"human"`` or ``"mouse"``).
      genomic_unit: Normalized genomic unit (``"gene"``, ``"exon"``, or
        ``"junction"``).
      annotation_label: Optional human-readable annotation label.
      annotation_extension: Optional explicit annotation extension.

    Returns:
      A single annotation extension string for gene/exon assays, or
      :data:`None` when no annotation is required (for example,
      junction-level assays).
    """
    if genomic_unit not in {"gene", "exon"}:
        return None

    if annotation_extension:
        return str(annotation_extension).strip()

    if annotation_label:
        return search.annotation_ext(
            organism=organism, annotation=annotation_label
        )

    org = organism.strip().lower()
    if org == "human":
        return "G026"
    if org == "mouse":
        return "M023"

    raise ValueError(
        f"Unsupported organism {organism!r}; expected 'human' or 'mouse'."
    )


def build_summarized_experiment(
    bundle: r3_bundle.R3ResourceBundle,
    *,
    genomic_unit: str,
    annotation_extension: str | None = None,
    assay_name: str = "raw_counts",
    join_policy: str = "inner",
    metadata_join: str = "inner",
    autoload: bool = True,
) -> summarizedexperiment.SummarizedExperiment:
    """Create a :class:`~summarizedexperiment.SummarizedExperiment`.

    Builds the experiment from a resource bundle.

    This is a convenience wrapper around
    :meth:`recount3.bundle.R3ResourceBundle.to_summarized_experiment`.

    Args:
      bundle: Resource bundle containing counts and metadata.
      genomic_unit: One of ``"gene"``, ``"exon"``, or ``"junction"``.
      annotation_extension: Optional annotation code for gene/exon
        assays (for example, ``"G026"``).
      assay_name: Name for the count assay in the SE.
      join_policy: ``"inner"`` intersects features; ``"outer"`` unions them
        and inserts zeros for features absent from an input project.
      metadata_join: ``"inner"`` intersects nonempty
        metadata tables within each project; ``"outer"`` retains every
        count sample.
      autoload: If :data:`True`, load resources transparently. If False,
        counts and metadata must already be loaded.

    Returns:
      A :class:`summarizedexperiment.SummarizedExperiment` instance.

    Raises:
      ImportError: If the required BiocPy packages are unavailable.
      ValueError: If counts, identifiers, dimensions, or join policies
        are invalid, or no counts match the selection.
      recount3.errors.CompatibilityError: If resources cannot be combined
        unambiguously.
      recount3.errors.RangesError: If required junction coordinates cannot
        be aligned, or requested ranges fail without an enabled fallback.
      recount3.errors.LoadError: If required metadata cannot be loaded.
      recount3.errors.DownloadError: If required resource retrieval fails.
    """
    unit = _utils._normalize_genomic_unit(genomic_unit)
    return bundle.to_summarized_experiment(
        genomic_unit=unit,
        annotation_extension=annotation_extension,
        assay_name=assay_name,
        join_policy=join_policy,
        metadata_join=metadata_join,
        autoload=autoload,
    )


def build_ranged_summarized_experiment(
    bundle: r3_bundle.R3ResourceBundle,
    *,
    genomic_unit: str,
    annotation_extension: str | None = None,
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
    """Create a :class:`~summarizedexperiment.RangedSummarizedExperiment`.

    Returns a ranged experiment when genomic ranges can be resolved.
    This is a convenience wrapper around
    :meth:`recount3.bundle.R3ResourceBundle.to_ranged_summarized_experiment`.

    Args:
      bundle: Resources for counts, metadata, and (optionally)
        annotations.
      genomic_unit: One of ``"gene"``, ``"exon"``, or ``"junction"``.
      annotation_extension: Annotation code for gene/exon, if
        desired.
      prefer_rr_junction_coordinates: Whether to use RR junction sidecars.
        Disabling this option prevents junction range construction and raises
        a range error unless SE fallback is enabled. Ignored for genes/exons.
      assay_name: Name for the count assay in the output.
      join_policy: ``"inner"`` intersects features; ``"outer"`` unions them
        and inserts zeros for features absent from an input project.
      metadata_join: ``"inner"`` intersects nonempty
        metadata tables within each project; ``"outer"`` retains every
        count sample.
      autoload: If :data:`True`, download and load resources as needed.
        If False, counts and metadata must already be loaded, and range
        files must be cached.
      allow_fallback_to_se: If :data:`True`, return a plain SE when
        ranges are unavailable. That object has no genomic ranges, so
        operations requiring an RSE will not work on it; the flag changes
        the return type on failure and neither retries a failed retrieval
        nor repairs a mismatched annotation.

    Returns:
      A :class:`summarizedexperiment.RangedSummarizedExperiment` object,
      or a plain :class:`summarizedexperiment.SummarizedExperiment` when
      ``allow_fallback_to_se`` is :data:`True` and ranges cannot be
      resolved.

    Raises:
      ImportError: If the required BiocPy packages are unavailable.
      ValueError: If counts, identifiers, dimensions, or join policies
        are invalid, or no counts match the selection.
      recount3.errors.CompatibilityError: If resources cannot be combined
        unambiguously.
      recount3.errors.RangesError: If required junction coordinates cannot
        be aligned, or requested ranges fail without an enabled fallback.
      recount3.errors.LoadError: If required metadata cannot be loaded.
      recount3.errors.DownloadError: If required resource retrieval fails.
    """
    unit = _utils._normalize_genomic_unit(genomic_unit)
    return bundle.to_ranged_summarized_experiment(
        genomic_unit=unit,
        annotation_extension=annotation_extension,
        prefer_rr_junction_coordinates=prefer_rr_junction_coordinates,
        assay_name=assay_name,
        join_policy=join_policy,
        metadata_join=metadata_join,
        autoload=autoload,
        allow_fallback_to_se=allow_fallback_to_se,
    )


def create_ranged_summarized_experiment(
    *,
    project: str,
    genomic_unit: str = "gene",
    organism: str = "human",
    data_source: str = "sra",
    annotation_label: str | None = None,
    annotation_extension: str | None = None,
    junction_type: str = "ALL",
    junction_extensions: Sequence[str] | None = None,
    include_metadata: bool = True,
    include_bigwig: bool = False,
    assay_name: str = "raw_counts",
    join_policy: str = "inner",
    metadata_join: str = "inner",
    autoload: bool = True,
    allow_fallback_to_se: bool = False,
    strict: bool = True,
) -> (
    summarizedexperiment.RangedSummarizedExperiment
    | summarizedexperiment.SummarizedExperiment
):
    """Discover project resources and construct a ranged expression experiment.

    Count matrices and sample metadata are aligned by identity. Gene/exon
    ranges come from the selected GTF annotation; junction ranges use the
    matching RR sidecar. Sparse counts remain sparse in the output.

    Note:
        Discovery reads the project's resource listings. See
        :meth:`~recount3.bundle.R3ResourceBundle.to_ranged_summarized_experiment`
        for the cost of building the experiment. Downloads usually dominate
        when files are uncached.

    Args:
      project: Study or project identifier (for example, ``"SRP009615"``).
      genomic_unit: One of ``{"gene", "exon", "junction"}``.
      organism: Organism identifier (``"human"`` or ``"mouse"``).
      data_source: Data source (``"sra"``, ``"gtex"``, or ``"tcga"``).
      annotation_label: Optional human-readable annotation label for gene/exon
        assays, such as ``"gencode_v26"`` or ``"gencode_v29"``. Ignored
        for junction-level assays.
      annotation_extension: Optional explicit annotation extension
        (for example, ``"G026"``). When provided, this takes precedence
        over ``annotation_label``. Ignored for junction-level assays.
      junction_type: Junction type; typically ``"ALL"``.
      junction_extensions: Iterable of junction artifact extensions
        to include (for example, ``("MM",)`` or ``("MM", "RR")``). If
        :data:`None`, the default depends on ``genomic_unit``:
        ``("MM", "RR")`` for ``"junction"`` (so genomic ranges can be
        attached to each row), and ``("MM",)`` for ``"gene"`` / ``"exon"``
        (the RR sidecar does not apply to those units).
      include_metadata: Whether to include project metadata tables in
        the bundle (five for SRA, four for GTEx/TCGA; recommended).
      include_bigwig: Whether to include per-sample BigWig coverage
        resources in the bundle. These can be large.
      assay_name: Count assay name. The default ``"raw_counts"`` becomes
        ``"counts"`` for junction experiments.
      join_policy: Join policy across projects when stacking count matrices
        (passed to :func:`build_ranged_summarized_experiment`).
      metadata_join: ``"inner"`` intersects nonempty
        metadata tables within each project; ``"outer"`` retains every
        count sample.
      autoload: If :data:`True`, download and load resources as needed.
        If False, counts and metadata must already be loaded, and range
        files must be cached.
      allow_fallback_to_se: If :data:`True`, construct a plain SE when
        genomic ranges cannot be derived for the requested combination.
        That object has no genomic ranges, so operations requiring an RSE
        will not work on it; the flag changes the return type on failure
        and neither retries a failed retrieval nor repairs a mismatched
        annotation.
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
      recount3.errors.CompatibilityError: If selected resources cannot be
        combined unambiguously.
      recount3.errors.RangesError: If ranges cannot be resolved and fallback
        is disabled, or required junction alignment fails.
      recount3.errors.LoadError: If required metadata cannot be loaded.
      recount3.errors.DownloadError: If retrieval fails during count or
        metadata preparation.
      TypeError: If the underlying
        :class:`RangedSummarizedExperiment` constructor rejects input types.

    Examples:
        Build an RSE for a human SRA gene-count project (GENCODE 26)::

            rse = create_ranged_summarized_experiment(
                project="SRP009615",
                organism="human",
                annotation_label="gencode_v26",
            )

        Use the raw annotation extension instead of a label::

            rse = create_ranged_summarized_experiment(
                project="SRP009615",
                organism="human",
                annotation_extension="G026",
            )

        Build a junction-level RSE::

            rse = create_ranged_summarized_experiment(
                project="SRP009615",
                organism="human",
                genomic_unit="junction",
            )
    """
    unit = _utils._normalize_genomic_unit(genomic_unit)
    ann_ext = _resolve_annotation_extension(
        organism=organism,
        genomic_unit=unit,
        annotation_label=annotation_label,
        annotation_extension=annotation_extension,
    )

    if junction_extensions is None:
        junction_exts: tuple[str, ...] = (
            ("MM", "RR") if unit == "junction" else ("MM",)
        )
    else:
        junction_exts = tuple(str(ext).strip() for ext in junction_extensions)

    project_bundle = r3_bundle.R3ResourceBundle.discover(
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
        annotation_extension=ann_ext,
        prefer_rr_junction_coordinates=True,
        assay_name=assay_name,
        join_policy=join_policy,
        metadata_join=metadata_join,
        autoload=autoload,
        allow_fallback_to_se=allow_fallback_to_se,
    )


def create_rse(
    *,
    project: str,
    genomic_unit: str = "gene",
    organism: str = "human",
    data_source: str = "sra",
    annotation_label: str | None = None,
    annotation_extension: str | None = None,
    junction_type: str = "ALL",
    junction_extensions: Sequence[str] | None = None,
    include_metadata: bool = True,
    include_bigwig: bool = False,
    assay_name: str = "raw_counts",
    join_policy: str = "inner",
    metadata_join: str = "inner",
    autoload: bool = True,
    allow_fallback_to_se: bool = False,
    strict: bool = True,
) -> (
    summarizedexperiment.RangedSummarizedExperiment
    | summarizedexperiment.SummarizedExperiment
):
    """Discover project resources and construct a ranged expression experiment.

    Convenience alias for :func:`create_ranged_summarized_experiment`.

    Count matrices and sample metadata are aligned by identity. Gene/exon
    ranges come from the selected GTF annotation; junction ranges use the
    matching RR sidecar. Sparse counts remain sparse in the output.

    Note:
        Discovery reads the project's resource listings. See
        :meth:`~recount3.bundle.R3ResourceBundle.to_ranged_summarized_experiment`
        for the cost of building the experiment. Downloads usually dominate
        when files are uncached.

    Args:
      project: Study or project identifier (for example, ``"SRP009615"``).
      genomic_unit: One of ``{"gene", "exon", "junction"}``.
      organism: Organism identifier (``"human"`` or ``"mouse"``).
      data_source: Data source (``"sra"``, ``"gtex"``, or ``"tcga"``).
      annotation_label: Optional human-readable annotation label for gene/exon
        assays, such as ``"gencode_v26"`` or ``"gencode_v29"``. Ignored
        for junction-level assays.
      annotation_extension: Optional explicit annotation extension
        (for example, ``"G026"``). When provided, this takes precedence
        over ``annotation_label``. Ignored for junction-level assays.
      junction_type: Junction type; typically ``"ALL"``.
      junction_extensions: Iterable of junction artifact extensions
        to include (for example, ``("MM",)`` or ``("MM", "RR")``). If
        :data:`None`, the default depends on ``genomic_unit``:
        ``("MM", "RR")`` for ``"junction"`` (so genomic ranges can be
        attached to each row), and ``("MM",)`` for ``"gene"`` / ``"exon"``
        (the RR sidecar does not apply to those units).
      include_metadata: Whether to include project metadata tables in
        the bundle (five for SRA, four for GTEx/TCGA; recommended).
      include_bigwig: Whether to include per-sample BigWig coverage
        resources in the bundle. These can be large.
      assay_name: Count assay name. The default ``"raw_counts"`` becomes
        ``"counts"`` for junction experiments.
      join_policy: Join policy across projects when stacking count matrices
        (passed to :func:`build_ranged_summarized_experiment`).
      metadata_join: ``"inner"`` intersects nonempty
        metadata tables within each project; ``"outer"`` retains every
        count sample.
      autoload: If :data:`True`, download and load resources as needed.
        If False, counts and metadata must already be loaded, and range
        files must be cached.
      allow_fallback_to_se: If :data:`True`, construct a plain SE when
        genomic ranges cannot be derived for the requested combination.
        That object has no genomic ranges, so operations requiring an RSE
        will not work on it; the flag changes the return type on failure
        and neither retries a failed retrieval nor repairs a mismatched
        annotation.
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
      recount3.errors.CompatibilityError: If selected resources cannot be
        combined unambiguously.
      recount3.errors.RangesError: If ranges cannot be resolved and fallback
        is disabled, or required junction alignment fails.
      recount3.errors.LoadError: If required metadata cannot be loaded.
      recount3.errors.DownloadError: If retrieval fails during count or
        metadata preparation.
      TypeError: If the underlying
        :class:`RangedSummarizedExperiment` constructor rejects input types.

    Examples:
        Build an RSE for a human SRA gene-count project (GENCODE 26)::

            rse = create_ranged_summarized_experiment(
                project="SRP009615",
                organism="human",
                annotation_label="gencode_v26",
            )

        Use the raw annotation extension instead of a label::

            rse = create_ranged_summarized_experiment(
                project="SRP009615",
                organism="human",
                annotation_extension="G026",
            )

        Build a junction-level RSE::

            rse = create_ranged_summarized_experiment(
                project="SRP009615",
                organism="human",
                genomic_unit="junction",
            )
    """
    return create_ranged_summarized_experiment(
        project=project,
        genomic_unit=genomic_unit,
        organism=organism,
        data_source=data_source,
        annotation_label=annotation_label,
        annotation_extension=annotation_extension,
        junction_type=junction_type,
        junction_extensions=junction_extensions,
        include_metadata=include_metadata,
        include_bigwig=include_bigwig,
        assay_name=assay_name,
        join_policy=join_policy,
        metadata_join=metadata_join,
        autoload=autoload,
        allow_fallback_to_se=allow_fallback_to_se,
        strict=strict,
    )


def expand_sra_attributes(
    experiment_or_coldata: (
        summarizedexperiment.RangedSummarizedExperiment
        | summarizedexperiment.SummarizedExperiment
        | pd.DataFrame
    ),
    *,
    sra_attributes_column: str = "sra.sample_attributes",
    attribute_column_prefix: str = "sra_attribute.",
) -> (
    summarizedexperiment.RangedSummarizedExperiment
    | summarizedexperiment.SummarizedExperiment
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
    ``{attribute_column_prefix}{key}`` (with spaces in ``key`` replaced by
    ``'_'``), and the parsed values are stored per sample. The original column
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
      experiment_or_coldata: Column metadata DataFrame or a BiocPy SE/RSE
        object.
      sra_attributes_column: Name of the column that holds the raw SRA attribute
        strings.
      attribute_column_prefix: Prefix to prepend to the generated attribute
        column names.

    Returns:
      A new object of the same type as ``experiment_or_coldata`` with
      additional columns corresponding to parsed SRA attributes. If the
      requested column is missing, ``experiment_or_coldata`` is returned
      unchanged.

    Raises:
      ImportError: If a :class:`~summarizedexperiment.SummarizedExperiment` /
        :class:`~summarizedexperiment.RangedSummarizedExperiment` is supplied
        but BiocPy packages are not installed.
      AttributeError: If the BiocPy object does not expose
        ``column_data`` or ``set_column_data`` in the expected API.
      TypeError: If ``experiment_or_coldata`` is neither a
        :class:`~pandas.DataFrame` nor a supported BiocPy experiment object.

    Note:
        ``O(L + S * (M + K))`` for ``S`` samples, ``M`` metadata columns,
        ``L`` characters of attribute text and ``K`` distinct attribute
        keys: the size of the input plus the output. For an experiment the
        sample metadata is converted to pandas and back, ``O(S * M)`` each
        way.

    Examples:
        Expand attributes on a metadata DataFrame::

            expanded_df = r3.se.expand_sra_attributes(col_data_df)

        Expand attributes directly on an RSE (returns a new RSE)::

            rse2 = r3.se.expand_sra_attributes(rse)

        Inspect the new attribute columns::

            attr_cols = [c for c in rse2.column_data.column_names
                         if c.startswith("sra_attribute.")]
    """
    # DataFrame mode: no BiocPy dependency.
    if isinstance(experiment_or_coldata, pd.DataFrame):
        return _expand_sra_attributes_df(
            experiment_or_coldata,
            sra_attributes_column=sra_attributes_column,
            attribute_column_prefix=attribute_column_prefix,
        )

    bioc_frame_cls = _utils.get_biocframe_class()
    summarized_experiment_cls = _utils.get_summarizedexperiment_class()

    if not isinstance(
        experiment_or_coldata,
        summarized_experiment_cls,
    ):
        raise TypeError(
            "expand_sra_attributes() expects a pandas DataFrame or a BiocPy "
            "SummarizedExperiment/RangedSummarizedExperiment; got "
            f"{type(experiment_or_coldata)!r}."
        )

    col_bf = experiment_or_coldata.column_data
    try:
        col_df = col_bf.to_pandas()  # type: ignore[attr-defined]
    except AttributeError as exc:
        raise AttributeError(
            "Expected column_data to implement to_pandas(); got "
            f"{type(col_bf)!r} instead."
        ) from exc

    resolved = _resolve_sra_attributes_column(
        col_df.columns, sra_attributes_column
    )
    if resolved is None:
        logging.warning(
            "expand_sra_attributes: column %r (and its ./__ variant) "
            "not present in column_data; returning object unchanged.",
            sra_attributes_column,
        )
        return experiment_or_coldata

    expanded_df = _expand_sra_attributes_df(
        col_df,
        sra_attributes_column=resolved,
        attribute_column_prefix=attribute_column_prefix,
    )
    expanded_bf = bioc_frame_cls.from_pandas(expanded_df)
    return experiment_or_coldata.set_column_data(expanded_bf)


def compute_read_counts(
    rse: summarizedexperiment.RangedSummarizedExperiment,
    round_to_integers: bool = True,
    avg_mapped_read_length_column: str = (
        "recount_qc.star.average_mapped_length"
    ),
) -> pd.DataFrame:
    """Forward to :func:`recount3.compute.compute_read_counts`.

    .. deprecated:: 1.4.0
       Use :func:`recount3.compute.compute_read_counts`.
    """
    warnings.warn(
        "recount3.se.compute_read_counts is deprecated and will be removed "
        "in v2.0.0. Use recount3.compute.compute_read_counts instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return compute.compute_read_counts(
        rse,
        round_to_integers,
        avg_mapped_read_length_column,
    )


def compute_tpm(
    rse: summarizedexperiment.RangedSummarizedExperiment,
) -> pd.DataFrame:
    """Forward to :func:`recount3.compute.compute_tpm`.

    .. deprecated:: 1.4.0
       Use :func:`recount3.compute.compute_tpm`.
    """
    warnings.warn(
        "recount3.se.compute_tpm is deprecated and will be removed "
        "in v2.0.0. Use recount3.compute.compute_tpm instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return compute.compute_tpm(rse)


def is_paired_end(
    sample_metadata_source: Any,
    avg_mapped_read_length_column: str = (
        "recount_qc.star.average_mapped_length"
    ),
    avg_read_length_column: str = compute._AVG_READ_LENGTH_COLUMN,
) -> pd.Series:
    """Forward to :func:`recount3.compute.is_paired_end`.

    .. deprecated:: 1.4.0
       Use :func:`recount3.compute.is_paired_end`.
    """
    warnings.warn(
        "recount3.se.is_paired_end is deprecated and will be removed "
        "in v2.0.0. Use recount3.compute.is_paired_end instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return compute.is_paired_end(
        sample_metadata_source,
        avg_mapped_read_length_column,
        avg_read_length_column,
    )


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
    """Forward to :func:`recount3.compute.compute_scale_factors`.

    .. deprecated:: 1.4.0
       Use :func:`recount3.compute.compute_scale_factors`.
    """
    warnings.warn(
        "recount3.se.compute_scale_factors is deprecated and will be removed "
        "in v2.0.0. Use recount3.compute.compute_scale_factors instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return compute.compute_scale_factors(
        sample_metadata_source,
        by,
        target_read_count,
        target_read_length_bp,
        auc_column,
        avg_mapped_read_length_column,
        mapped_reads_column,
        paired_end_status,
    )


def transform_counts(
    rse: summarizedexperiment.RangedSummarizedExperiment,
    by: str = "auc",
    target_read_count: float = 4e7,
    target_read_length_bp: float = 100,
    round_to_integers: bool = True,
    **kwargs: Any,
) -> pd.DataFrame:
    """Forward to :func:`recount3.compute.transform_counts`.

    .. deprecated:: 1.4.0
       Use :func:`recount3.compute.transform_counts`.
    """
    warnings.warn(
        "recount3.se.transform_counts is deprecated and will be removed "
        "in v2.0.0. Use recount3.compute.transform_counts instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return compute.transform_counts(
        rse,
        by,
        target_read_count,
        target_read_length_bp,
        round_to_integers,
        **kwargs,
    )


def to_anndata(experiment: Any, *, sanitize_for_hdf5: bool = False) -> Any:
    """Convert a SummarizedExperiment to AnnData, keeping its provenance.

    AnnData is the container most Python single-cell and machine-learning
    tooling reads, and it is transposed relative to BiocPy: samples are rows
    (``obs``) and features are columns (``var``), with each recount3 assay
    kept as a layer.

    BiocPy's own ``experiment.to_anndata()`` cannot be used on a recount3
    experiment. It forwards experiment metadata straight to
    ``AnnData(uns=...)``, and BiocPy stores that metadata as a
    :class:`~biocutils.NamedList.NamedList`, which AnnData rejects with
    ``Only mutable mapping types (e.g. dict) are allowed for `.uns`.``. This
    function converts the provenance to plain dicts and lists first, so the
    resulting object carries the same ``project``, ``annotation``,
    ``resource_urls`` and ``metadata_columns`` record.

    Writing the result to HDF5 has a separate requirement, which
    ``sanitize_for_hdf5`` handles: recount3 STAR QC fields are named after
    splice motifs and contain ``/``, which HDF5 reads as a path separator,
    both in the sample column names and in the ``uns`` provenance map keyed
    by them. Renaming them is a real change to what an analysis indexes by,
    so it is opt-in rather than automatic.

    Args:
        experiment: A :class:`~summarizedexperiment.SummarizedExperiment` or
          :class:`~summarizedexperiment.RangedSummarizedExperiment`, for
          example from :func:`create_rse`. The object is not modified.
        sanitize_for_hdf5: If :data:`True`, make the result writable by
          :meth:`anndata.AnnData.write_h5ad`: cast all-missing object columns
          to NaN, and replace ``/`` with ``_`` in ``obs``/``var`` column names
          and in nested ``uns`` keys. Each renamed provenance entry keeps its
          original name in its value, so nothing is lost. This is what
          ``recount3 bundle rse --sanitize-columns`` applies.

    Returns:
        An ``anndata.AnnData`` with samples in rows, features in columns,
        each assay as a layer, and the experiment's provenance in ``uns``.

    Raises:
        ImportError: If the ``anndata`` extra is not installed.
        ValueError: If ``sanitize_for_hdf5`` is requested and a sanitized
          name would collide with an existing one, which would silently
          merge two distinct fields.

    Note:
        The size of the experiment: its assays (``O(F * S)`` dense, or
        ``O(N)`` for ``N`` stored values of a sparse assay), its row and
        column data, ``O(F * K + S * M)``, and its provenance. Sanitizing
        adds one pass over every column and ``uns`` key.

    Examples:
        Convert an RSE and inspect the orientation::

            rse = r3.create_rse(project="SRP009615", organism="human",
                                annotation_label="gencode_v26")
            adata = r3.se.to_anndata(rse)
            adata.shape            # (12, 63856): samples x genes
            adata.uns["project"]   # 'SRP009615'

        Convert and write an HDF5 file::

            adata = r3.se.to_anndata(rse, sanitize_for_hdf5=True)
            adata.write_h5ad("rse.h5ad")
    """
    _utils.ensure_anndata_support()
    adata = _utils.experiment_to_anndata(experiment)
    if sanitize_for_hdf5:
        _utils.normalize_anndata_for_hdf5(adata)
        _utils.sanitize_anndata_column_names(adata)
        _utils.sanitize_anndata_uns_keys(adata)
    return adata
