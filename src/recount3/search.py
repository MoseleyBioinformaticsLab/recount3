"""Discovery and search helpers.

These functions generate :class:`R3Resource` objects for common queries and
preserve the original behavior.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable
from typing import Final

import pandas as pd

from ._descriptions import (
    R3ResourceDescription,
    _VALID_DATA_SOURCES,
)
from .resource import R3Resource
from .types import FieldSpec, StringOrIterable


def _as_tuple(value: StringOrIterable) -> tuple[str, ...]:
    """Normalize a string or iterable of strings into a tuple of strings."""
    if isinstance(value, str):
        return (value,)
    return tuple(value)


def _match_spec(value: object | None, spec: FieldSpec) -> bool:
    """Return True if ``value`` satisfies the selection ``spec``."""
    if spec is None:
        return True
    if callable(spec):
        return bool(spec(value))
    return _as_tuple(spec) and value in _as_tuple(spec)  # type: ignore[reportReturnType]


def _build_param_grid(resource_type: str, **required_values: StringOrIterable) -> list[dict[str, str]]:
    """Cartesian product over provided values."""
    keys = list(required_values.keys())
    vals = [_as_tuple(required_values[k]) for k in keys]
    grid: list[dict[str, str]] = []
    for combo in itertools.product(*vals):
        d: dict[str, str] = {"resource_type": resource_type}
        d.update({k: v for k, v in zip(keys, combo)})
        grid.append(d)
    return grid


def _make_resources(
    param_dicts: Iterable[dict[str, str]],
    *,
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Instantiate R3Resource for each param dict, with optional de-duplication."""
    out: list[R3Resource] = []
    seen: set[str] = set()
    for params in param_dicts:
        try:
            desc = R3ResourceDescription(**params)
            res = R3Resource(desc)
        except Exception:
            if strict:
                raise
            continue
        if deduplicate:
            if res.url in seen:
                continue
            if res.url is not None:
                seen.add(res.url)
        out.append(res)
    return out


# ---- Type-specific search wrappers -----------------------------------

def search_annotations(
    *,
    organism: StringOrIterable,
    genomic_unit: StringOrIterable,
    annotation_file_extension: StringOrIterable,
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Return resources for annotation files."""
    grid = _build_param_grid(
        "annotations",
        organism=organism,
        genomic_unit=genomic_unit,
        annotation_file_extension=annotation_file_extension,
    )
    return _make_resources(grid, strict=strict, deduplicate=deduplicate)


def search_count_files_gene_or_exon(
    *,
    organism: StringOrIterable,
    data_source: StringOrIterable,
    genomic_unit: StringOrIterable,
    project: StringOrIterable,
    annotation_file_extension: StringOrIterable = ("G026",),
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Return resources for per-project gene/exon counts."""
    grid = _build_param_grid(
        "count_files_gene_or_exon",
        organism=organism,
        data_source=data_source,
        genomic_unit=genomic_unit,
        project=project,
        annotation_file_extension=annotation_file_extension,
    )
    return _make_resources(grid, strict=strict, deduplicate=deduplicate)


def search_count_files_junctions(
    *,
    organism: StringOrIterable,
    data_source: StringOrIterable,
    project: StringOrIterable,
    junction_type: StringOrIterable = "ALL",
    junction_file_extension: StringOrIterable = "MM",
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Return resources for per-project junction count files."""
    grid = _build_param_grid(
        "count_files_junctions",
        organism=organism,
        data_source=data_source,
        project=project,
        junction_type=junction_type,
        junction_file_extension=junction_file_extension,
    )
    return _make_resources(grid, strict=strict, deduplicate=deduplicate)


def search_metadata_files(
    *,
    organism: StringOrIterable,
    data_source: StringOrIterable,
    table_name: StringOrIterable,
    project: StringOrIterable,
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Return resources for per-project metadata tables."""
    grid = _build_param_grid(
        "metadata_files",
        organism=organism,
        data_source=data_source,
        table_name=table_name,
        project=project,
    )
    return _make_resources(grid, strict=strict, deduplicate=deduplicate)


def search_bigwig_files(
    *,
    organism: StringOrIterable,
    data_source: StringOrIterable,
    project: StringOrIterable,
    sample: StringOrIterable,
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Return resources for BigWig files (download-only; load() returns a reader)."""
    grid = _build_param_grid(
        "bigwig_files",
        organism=organism,
        data_source=data_source,
        project=project,
        sample=sample,
    )
    return _make_resources(grid, strict=strict, deduplicate=deduplicate)


def search_data_sources(
    *,
    organism: StringOrIterable,
    strict: bool = True,
) -> list[R3Resource]:
    """Return an R3Resource representing the organism-level data-source index."""
    grid = _build_param_grid("data_sources", organism=organism)
    return _make_resources(grid, strict=strict, deduplicate=True)


def search_data_source_metadata(
    *,
    organism: StringOrIterable,
    data_source: StringOrIterable,
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Return resources for data-source-level metadata listings."""
    grid = _build_param_grid(
        "data_source_metadata",
        organism=organism,
        data_source=data_source,
    )
    return _make_resources(grid, strict=strict, deduplicate=deduplicate)


# ---------------------------------------------------------------------------

def create_sample_project_lists(organism: str = "") -> tuple[list[str], list[str]]:
    """Return (samples, projects) discovered from per-data-source metadata.

    This follows the original approach by iterating known data sources and
    extracting sample/project identifiers from their metadata tables.

    Args:
      organism: Optional organism filter ("human" or "mouse", case-insensitive).

    Returns:
      (samples, projects): Sorted unique IDs.
    """
    import pandas as pd  # local, to keep module import cost small

    data_sources = sorted(_VALID_DATA_SOURCES)

    samples: set[str] = set()
    projects: set[str] = set()

    def _rows(obj) -> list[dict]:
        if pd is not None and isinstance(obj, pd.DataFrame):
            return obj.to_dict("records")
        return list(obj)

    proj_keys = ("project", "project_id", "study", "study_accession")
    samp_keys = ("sample", "sample_id", "run", "run_accession", "external_id")
    org_keys = ("organism", "species")

    for ds in data_sources:
        try:
            res = R3Resource(
                R3ResourceDescription(
                    resource_type="data_source_metadata",
                    organism="human",
                    data_source=ds,
                )
            )
            meta = res.load()
        except Exception:
            continue

        for row in _rows(meta):
            if organism:
                ov = None
                for k in org_keys:
                    if k in row and row[k]:
                        ov = str(row[k]).lower()
                        break
                if ov and ov != organism.lower():
                    continue

            pv = None
            for k in proj_keys:
                if k in row and row[k]:
                    pv = str(row[k])
                    break
            if pv:
                projects.add(pv)

            sv = None
            for k in samp_keys:
                if k in row and row[k]:
                    sv = str(row[k])
                    break
            if sv:
                samples.add(sv)

    return sorted(samples), sorted(projects)

#   https://rna.recount.bio/docs/raw-files.html  (Section 6.2)
# Human:  G026 (Gencode v26), G029 (Gencode v29), F006 (FANTOM6_cat),
#         R109 (RefSeq), ERCC (ERCC), SIRV (SIRV)
# Mouse:  M023 (Gencode v23)
_ANN_EXT_HUMAN: Final[tuple[str, ...]] = (
    "G026",
    "G029",
    "F006",
    "R109",
    "ERCC",
    "SIRV",
)
_ANN_EXT_MOUSE: Final[tuple[str, ...]] = ("M023",)


def _resolve_annotation_exts(
    organism: str,
    annotations: str | Iterable[str] | None,
    explicit_exts: Iterable[str] | None,
) -> tuple[str, ...]:
    """Resolve which annotation file extensions to include.

    Args:
      organism: "human" or "mouse".
      annotations: Either "default", "all", or an explicit iterable of
        extensions (for example, ("G026", "G029")).
      explicit_exts: Optional explicit extensions (typically from CLI's
        ``annotation_file_extension=...``). If provided, these win.

    Returns:
      A tuple of annotation extensions.

    Raises:
      ValueError: If the organism is invalid.
    """
    org = organism.strip().lower()
    if org not in ("human", "mouse"):
        raise ValueError(f"Invalid organism: {organism!r}")

    if explicit_exts:
        return tuple(str(x).strip() for x in explicit_exts)

    if annotations is None or annotations == "default":
        return ("G026",) if org == "human" else ("M023",)

    if isinstance(annotations, str):
        if annotations == "all":
            return _ANN_EXT_HUMAN if org == "human" else _ANN_EXT_MOUSE
        # Allow comma-separated strings.
        parts = [p.strip() for p in annotations.split(",") if p.strip()]
        return tuple(parts)

    return tuple(str(x).strip() for x in annotations)


def samples_for_project(
    *,
    organism: str,
    data_source: str,
    project: str,
) -> list[str]:
    """Return sample identifiers for a given project.

    This helper reads the *data-source-level* metadata table and extracts
    sample IDs for the requested project by consulting common column names:
    "sample", "sample_id", "run", "run_accession", and "external_id".

    The logic mirrors the R docs recommendation to consult project-level
    metadata first, while mining samples from data-source metadata when
    assembling coverage files per sample.

    Args:
      organism: "human" or "mouse".
      data_source: "sra", "gtex", or "tcga".
      project: Study identifier (for example, "SRP009615").

    Returns:
      Sorted list of unique sample identifiers.

    Raises:
      ValueError: If the project cannot be validated against the metadata.
    """
    # Search and load the *source-level* metadata listing (single resource).
    ds_meta = search_data_source_metadata(
        organism=organism,
        data_source=data_source,
    )
    if not ds_meta:
        raise ValueError(
            "No data-source metadata found; cannot resolve samples."
        )
    df = ds_meta[0].load()
    if not isinstance(df, pd.DataFrame) or df.empty:
        raise ValueError("Empty data-source metadata; cannot resolve samples.")

    # Accept a few common project keys and sample keys. These cover SRA/GTEx.
    proj_keys = ("project", "study", "project_id", "project_accession")
    samp_keys = ("sample", "sample_id", "run", "run_accession", "external_id")

    # Filter to rows matching the project using the first available key.
    proj_col = next((k for k in proj_keys if k in df.columns), None)
    if proj_col is None:
        raise ValueError("Project column not found in metadata table.")
    mask = df[proj_col].astype(str) == str(project)
    sub = df.loc[mask]

    if sub.empty:
        raise ValueError(f"Project {project!r} not found in metadata.")

    # Extract sample identifiers from the first matching sample key.
    samp_col = next((k for k in samp_keys if k in sub.columns), None)
    if samp_col is None:
        # Fall back to any column named like "*run*" or "*sample*".
        candidates = [c for c in sub.columns if "run" in c or "sample" in c]
        if not candidates:
            return []
        samp_col = candidates[0]

    values = sorted({str(x) for x in sub[samp_col].dropna().astype(str)})
    return values


def search_project_all(
    *,
    organism: str,
    data_source: str,
    project: str,
    genomic_units: Iterable[str] = ("gene", "exon"),
    annotations: str | Iterable[str] = "default",
    junction_type: str = "ALL",
    junction_file_extension: Iterable[str] = ("MM",),
    include_metadata: bool = True,
    include_bigwig: bool = False,
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Enumerate all files for a project (counts, junctions, metadata, bw).

    This function composes the existing search helpers to implement a
    one-shot, project-scoped discovery routine. It adheres to recount3's
    raw file layout documented at:
    https://rna.recount.bio/docs/raw-files.html (Sections 6.2–6.4). :contentReference[oaicite:1]{index=1}

    Args:
      organism: "human" or "mouse".
      data_source: "sra", "gtex", or "tcga".
      project: Study identifier (for example, "SRP009615").
      genomic_units: Which expression levels to include. Defaults to both.
      annotations: "default", "all", comma-separated string, or an iterable
        of annotation file extensions (for example, ("G026", "G029")).
      junction_type: Junction type; typically "ALL".
      junction_file_extension: Iterable of junction artifacts to include:
        "MM" (counts), "RR" (coordinates), "ID" (sample IDs).
      include_metadata: Whether to include the five metadata tables.
      include_bigwig: Whether to include per-sample BigWig coverage files.
      strict: If True, raise on invalid parameters; else skip broken items.
      deduplicate: If True, drop duplicates across resource families.

    Returns:
      A list of :class:`R3Resource` objects covering the requested bundle.

    Raises:
      ValueError: If validation fails (for example, missing project).
    """
    # Validate project existence by consulting the data-source metadata.
    _ = samples_for_project(
        organism=organism,
        data_source=data_source,
        project=project,
    )

    # Resolve which annotation extensions to include.
    ann_exts = _resolve_annotation_exts(
        organism=organism,
        annotations=annotations,
        explicit_exts=None,
    )

    found: list[R3Resource] = []

    # Counts: gene/exon matrices, one per requested annotation.
    include_gene = "gene" in set(u.lower() for u in genomic_units)
    include_exon = "exon" in set(u.lower() for u in genomic_units)
    if include_gene or include_exon:
        units: list[str] = []
        if include_gene:
            units.append("gene")
        if include_exon:
            units.append("exon")

        found += search_count_files_gene_or_exon(
            organism=organism,
            data_source=data_source,
            genomic_unit=tuple(units),
            project=project,
            annotation_file_extension=ann_exts,
            strict=strict,
            deduplicate=deduplicate,
        )

    # Junctions: counts (MM) and optional RR/ID artifacts.
    if junction_file_extension:
        found += search_count_files_junctions(
            organism=organism,
            data_source=data_source,
            project=project,
            junction_type=junction_type,
            junction_file_extension=tuple(junction_file_extension),
            strict=strict,
            deduplicate=deduplicate,
        )

    # Project-level metadata (five tables).
    if include_metadata:
        # The CLI allows a single table_name, but project mode wants all.
        table_names = (
            "recount_project",
            "recount_qc",
            "recount_seq_qc",
            "recount_pred",
            data_source,  # the source-specific "project_meta" table
        )
        found += search_metadata_files(
            organism=organism,
            data_source=data_source,
            table_name=table_names,
            project=project,
            strict=strict,
            deduplicate=deduplicate,
        )

    # Per-sample BigWigs (optional; can be very large).
    if include_bigwig:
        samp_ids = samples_for_project(
            organism=organism,
            data_source=data_source,
            project=project,
        )
        if samp_ids:
            found += search_bigwig_files(
                organism=organism,
                data_source=data_source,
                project=project,
                sample=tuple(samp_ids),
                strict=strict,
                deduplicate=deduplicate,
            )

    # Return in a deterministic, human-friendly order.
    def _sort_key(r: R3Resource) -> tuple:
        d = r.description
        return (
            getattr(d, "resource_type", ""),
            getattr(d, "genomic_unit", ""),
            getattr(d, "annotation_file_extension", ""),
            getattr(d, "junction_type", ""),
            getattr(d, "junction_file_extension", ""),
            getattr(d, "table_name", ""),
            getattr(d, "sample", ""),
        )

    return sorted(found, key=_sort_key)