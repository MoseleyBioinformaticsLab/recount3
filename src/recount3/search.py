"""Discovery and search helpers.

These functions generate :class:`R3Resource` objects for common queries and
preserve the original behavior.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable
from typing import Callable

from .descriptions import (
    R3ResourceDescription,
    _VALID_DATA_SOURCES,  # internal, but used as before
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
