"""Discovery and search helpers.

These functions generate :class:`R3Resource` objects for common queries and
preserve the original behavior.
"""

from __future__ import annotations

import itertools
import posixpath
from collections.abc import Iterable
from typing import Final

import pandas as pd

from recount3._descriptions import (
    VALID_DATA_SOURCES,
    VALID_ORGANISMS,
    R3ResourceDescription,
)
from recount3.resource import R3Resource
from recount3.types import FieldSpec, StringOrIterable


def as_tuple(value: StringOrIterable) -> tuple[str, ...]:
    """Normalize a string or iterable of strings into a tuple of strings."""
    if isinstance(value, str):
        return (value,)
    return tuple(value)


def _build_param_grid(
    resource_type: str, **required_values: StringOrIterable
) -> list[dict[str, str]]:
    """Cartesian product over provided values."""
    keys = list(required_values.keys())
    vals = [as_tuple(required_values[k]) for k in keys]
    grid: list[dict[str, str]] = []
    for combo in itertools.product(*vals):
        d: dict[str, str] = {"resource_type": resource_type}
        d.update({k: v for k, v in zip(keys, combo, strict=True)})
        grid.append(d)
    return grid


def _make_resources(
    param_dicts: Iterable[dict[str, str]],
    *,
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Build R3Resource objects from param dicts; deduplicate by URL."""
    out: list[R3Resource] = []
    seen: set[str] = set()
    for params in param_dicts:
        try:
            desc = R3ResourceDescription(**params)
            res = R3Resource(desc)
        except Exception:  # pylint: disable=broad-exception-caught
            if strict:
                raise
            continue
        if deduplicate:
            if res.url in seen:
                continue
            # Resources with url=None are never added to `seen`, so multiple
            # such resources all pass through even with deduplicate=True.
            if res.url is not None:
                seen.add(res.url)
        out.append(res)
    return out


def match_spec(value: object | None, spec: FieldSpec) -> bool:
    """Return True if ``value`` satisfies the selection ``spec``."""
    if spec is None:
        return True
    if callable(spec):
        return bool(spec(value))
    return value in as_tuple(spec)


# ---- Type-specific search wrappers -----------------------------------

def search_annotations(
    *,
    organism: StringOrIterable,
    genomic_unit: StringOrIterable,
    annotation_extension: StringOrIterable,
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Return resources for annotation files."""
    grid = _build_param_grid(
        "annotations",
        organism=organism,
        genomic_unit=genomic_unit,
        annotation_extension=annotation_extension,
    )
    return _make_resources(grid, strict=strict, deduplicate=deduplicate)


def search_count_files_gene_or_exon(
    *,
    organism: StringOrIterable,
    data_source: StringOrIterable,
    genomic_unit: StringOrIterable,
    project: StringOrIterable,
    annotation_extension: StringOrIterable = ("G026",),
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
        annotation_extension=annotation_extension,
    )
    return _make_resources(grid, strict=strict, deduplicate=deduplicate)


def search_count_files_junctions(
    *,
    organism: StringOrIterable,
    data_source: StringOrIterable,
    project: StringOrIterable,
    junction_type: StringOrIterable = "ALL",
    junction_extension: StringOrIterable = "MM",
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
        junction_extension=junction_extension,
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
    """Return resources for per-sample BigWig coverage files."""
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
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Return resources for the organism-level data-source index."""
    grid = _build_param_grid("data_sources", organism=organism)
    return _make_resources(grid, strict=strict, deduplicate=deduplicate)


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

def create_sample_project_lists(
    organism: str = "",
) -> tuple[list[str], list[str]]:
    """Return (samples, projects) discovered from metadata tables.

    This is a compatibility wrapper around :func:`available_samples` and
    :func:`available_projects`. It preserves the original ``ids`` CLI
    behavior (simple ID lists) but now benefits from the richer and more
    robust metadata parsing.

    Args:
      organism: Optional organism filter. Accepts "human" or "mouse"
        (case-insensitive). An empty string means "all supported
        organisms".

    Returns:
      A tuple (samples, projects), where each element is a sorted list
      of unique identifier strings.
    """
    # Collect per-organism sample and project tables.
    if organism:
        org = _normalize_organism_name(organism)
        samples_df = available_samples(organism=org, strict=False)
        projects_df = available_projects(organism=org, strict=False)
    else:
        sample_frames: list[pd.DataFrame] = []
        project_frames: list[pd.DataFrame] = []
        for org in sorted(VALID_ORGANISMS):
            sample_frames.append(
                available_samples(organism=org, strict=False)
            )
            project_frames.append(
                available_projects(organism=org, strict=False)
            )
        samples_df = (
            pd.concat(sample_frames, ignore_index=True)
            if sample_frames
            else pd.DataFrame()
        )
        projects_df = (
            pd.concat(project_frames, ignore_index=True)
            if project_frames
            else pd.DataFrame()
        )

    # Extract sample IDs, preferring "external_id" where available.
    sample_ids: list[str] = []
    if not samples_df.empty:
        for col in (
            "external_id",
            "sample",
            "sample_id",
            "run",
            "run_accession",
        ):
            if col in samples_df.columns:
                sample_ids = sorted(
                    {
                        str(x)
                        for x in samples_df[col]
                        .dropna()
                        .astype(str)
                    }
                )
                break

    # Extract project IDs from the project-level table.
    project_ids: list[str] = []
    if not projects_df.empty and "project" in projects_df.columns:
        project_ids = sorted(
            {
                str(x)
                for x in projects_df["project"]
                .dropna()
                .astype(str)
            }
        )

    return sample_ids, project_ids


def _normalize_organism_name(organism: str) -> str:
    """Return a canonical lowercase organism name.

    Args:
      organism: Organism name, typically "human" or "mouse".

    Returns:
      Canonical lowercase organism name ("human" or "mouse").

    Raises:
      ValueError: If the organism is not supported.
    """
    org = organism.strip().lower()
    if org not in VALID_ORGANISMS:
        expected = sorted(VALID_ORGANISMS)
        msg = (
            f"Unsupported organism {organism!r}; expected one of "
            f"{expected!r}."
        )
        raise ValueError(msg)
    return org


def _strip_md_prefix(df: pd.DataFrame) -> pd.DataFrame:
    """Strip recount3-style prefixes from metadata column names.

    Per-data-source metadata tables usually prefix columns with a
    source-specific string such as "sra." or "gtex.". This helper removes
    everything up to the last dot.

    Args:
      df: Input DataFrame.

    Returns:
      A copy of ``df`` with cleaned column names.
    """
    out = df.copy()
    out.columns = [str(c).rsplit(".", 1)[-1] for c in out.columns]
    return out


def _path_basename(path: str) -> str:
    """Return the basename of a POSIX-style path without trailing slash.

    Args:
      path: POSIX-like path string.

    Returns:
      Basename component of the path. If ``path`` is empty, returns an
      empty string.
    """
    if not path:
        return ""
    return posixpath.basename(path.rstrip("/"))


def _path_dirname(path: str) -> str:
    """Return the dirname of a POSIX-style path without trailing slash.

    Args:
      path: POSIX-like path string.

    Returns:
      Directory component of the path, or an empty string for top-level
      paths or empty input.
    """
    if not path:
        return ""
    return posixpath.dirname(path.rstrip("/"))


def available_samples(
    *,
    organism: str = "human",
    data_sources: StringOrIterable | None = None,
    strict: bool = True,
) -> pd.DataFrame:
    """Return a sample overview similar to recount3::available_samples().

    This reads per-data-source ``*.recount_project.MD.gz`` tables and
    returns a normalized DataFrame describing all samples for the requested
    organism.

    Args:
      organism: Organism to query ("human" or "mouse").
      data_sources: Optional subset of data sources to include (for example,
        "sra", "gtex", "tcga"). By default all known data sources are used.
      strict: If True, raise a ValueError when no metadata can be found.
        If False, return an empty DataFrame instead.

    Returns:
      A DataFrame with at least the following columns when available:

      * ``external_id``: Sample identifier in the original source.
      * ``project``: Project or study identifier.
      * ``organism``: Canonical organism label ("human" or "mouse").
      * ``file_source``: Origin of the raw data (basename only).
      * ``date_processed``: Processing date in YYYY-MM-DD format.
      * ``project_home``: recount3 project home path.
      * ``project_type``: High-level project type (for example,
        "data_sources" or "collections").

      Additional columns present in the raw metadata are preserved.

    Raises:
      ValueError: If inputs are invalid or no metadata resources are found
        and ``strict`` is True.
      RuntimeError: If metadata resources are found but all fail to load.
    """
    org = _normalize_organism_name(organism)

    if data_sources is None:
        selected_sources = sorted(VALID_DATA_SOURCES)
    else:
        selected_sources = list(as_tuple(data_sources))
        invalid = [
            source
            for source in selected_sources
            if source not in VALID_DATA_SOURCES
        ]
        if invalid:
            msg = (
                "Unsupported data_sources "
                f"{invalid!r}; expected a subset of "
                f"{sorted(VALID_DATA_SOURCES)!r}."
            )
            raise ValueError(msg)

    if not selected_sources:
        if strict:
            raise ValueError("No data_sources specified.")
        return pd.DataFrame()

    resources = search_data_source_metadata(
        organism=org,
        data_source=tuple(selected_sources),
        strict=strict,
        deduplicate=True,
    )

    if not resources:
        if strict:
            msg = (
                "No data-source metadata resources found for organism "
                f"{org!r} and data_sources {tuple(selected_sources)!r}."
            )
            raise ValueError(msg)
        return pd.DataFrame()

    frames: list[pd.DataFrame] = []
    load_errors: list[tuple[str | None, Exception]] = []

    for res in resources:
        try:
            obj = res.load()
        except Exception as exc:  # pylint: disable=broad-except
            load_errors.append((res.url, exc))
            continue

        if isinstance(obj, pd.DataFrame):
            frames.append(_strip_md_prefix(obj))

    if not frames:
        if load_errors:
            first_url, first_exc = load_errors[0]
            msg = (
                "Failed to load any data-source metadata table. "
                f"First error: {first_exc!r} (while loading {first_url!r})."
            )
            raise RuntimeError(msg) from first_exc
        if strict:
            raise ValueError(
                "Metadata resources were located but produced no usable "
                "tables."
            )
        return pd.DataFrame()

    samples = pd.concat(frames, axis=0, ignore_index=True)

    samples.columns = [str(c) for c in samples.columns]

    if "organism" in samples.columns:
        samples["organism"] = samples["organism"].astype(str).replace(
            {"Homo sapiens": "human", "Mus musculus": "mouse"}
        )
    else:
        samples["organism"] = org

    if "project" not in samples.columns and "study" in samples.columns:
        samples = samples.rename(columns={"study": "project"})
    elif "project" in samples.columns and "study" in samples.columns:
        if not samples["project"].equals(samples["study"]):
            msg = (
                "Metadata columns 'project' and 'study' are not identical; "
                "this violates recount3 expectations."
            )
            raise ValueError(msg)
        samples = samples.drop(columns=["study"])

    if "file_source" in samples.columns:
        samples["file_source"] = (
            samples["file_source"].astype(str).map(_path_basename)
        )

    if (
        "project_home" not in samples.columns
        and "metadata_source" in samples.columns
    ):
        samples["project_home"] = samples["metadata_source"]
        samples = samples.drop(columns=["metadata_source"])

    if "project_home" in samples.columns:
        samples["project_type"] = (
            samples["project_home"].astype(str).map(_path_dirname)
        )

    for col in ("rail_id",):
        if col in samples.columns:
            samples = samples.drop(columns=[col])

    preferred = [
        "external_id",
        "project",
        "organism",
        "file_source",
        "date_processed",
        "project_home",
        "project_type",
    ]
    present = [c for c in preferred if c in samples.columns]
    remaining = [c for c in samples.columns if c not in present]

    return samples.loc[:, present + remaining]


def available_projects(
    *,
    organism: str = "human",
    data_sources: StringOrIterable | None = None,
    strict: bool = True,
) -> pd.DataFrame:
    """Return a project overview like recount3::available_projects().

    This aggregates the sample-level metadata from :func:`available_samples`
    and summarizes it at the project level.

    Args:
      organism: Organism to query ("human" or "mouse").
      data_sources: Optional subset of data sources to include.
      strict: Passed through to :func:`available_samples`.

    Returns:
      A DataFrame with one row per project and at least:

      * ``project``: Project or study identifier.
      * ``organism``: Canonical organism label.
      * ``file_source``: Origin of the raw data (basename only).
      * ``project_home``: recount3 project home path.
      * ``project_type``: High-level project type (for example,
        "data_sources" or "collections").
      * ``n_samples``: Number of samples in the project.

      Additional project-level columns are preserved.
    """
    samples = available_samples(
        organism=organism,
        data_sources=data_sources,
        strict=strict,
    )

    if samples.empty:
        columns = [
            "project",
            "organism",
            "file_source",
            "project_home",
            "project_type",
            "n_samples",
        ]
        return pd.DataFrame(columns=columns)

    df = samples.copy()

    for col in ("external_id", "date_processed"):
        if col in df.columns:
            df = df.drop(columns=[col])

    if "project_type" in df.columns:
        df = df.drop(columns=["project_type"])

    key_cols = [
        col
        for col in ("project", "organism", "project_home")
        if col in df.columns
    ]

    if key_cols:
        projects = (
            df.groupby(key_cols, as_index=False, dropna=False)
            .first()
            .reset_index(drop=True)
        )
    else:
        projects = df.drop_duplicates().reset_index(drop=True)

    if "project_home" in projects.columns:
        projects["project_type"] = (
            projects["project_home"].astype(str).map(_path_dirname)
        )

    if key_cols:
        sample_keys = df[key_cols].astype(str).agg("_".join, axis=1)
        counts = sample_keys.value_counts()

        project_keys = (
            projects[key_cols]
            .astype(str)
            .agg("_".join, axis=1)
        )
        projects["n_samples"] = project_keys.map(counts).astype("Int64")
    else:
        projects["n_samples"] = pd.Series(dtype="Int64")

    preferred = [
        "project",
        "organism",
        "file_source",
        "project_home",
        "project_type",
        "n_samples",
    ]
    present = [c for c in preferred if c in projects.columns]
    remaining = [c for c in projects.columns if c not in present]

    return projects.loc[:, present + remaining]


def project_homes(
    *,
    organism: str = "human",
    data_sources: StringOrIterable | None = None,
    strict: bool = True,
) -> pd.DataFrame:
    """Return a project home summary similar to recount3::project_homes().

    This is a thin layer on top of :func:`available_projects` that
    collapses projects down to unique ``project_home`` paths.

    Args:
      organism: Organism to query ("human" or "mouse").
      data_sources: Optional subset of data sources to include.
      strict: Passed through to :func:`available_projects`.

    Returns:
      A DataFrame with one row per project home and at least the
      following columns:

      * ``project_home``: recount3 project home path.
      * ``project_type``: High-level project type (for example,
        "data_sources" or "collections").
      * ``organism``: Canonical organism label.
      * ``file_source``: Data source label when available.
      * ``n_projects``: Number of projects using this home.

      Additional columns from :func:`available_projects` may appear.
    """
    projects = available_projects(
        organism=organism,
        data_sources=data_sources,
        strict=strict,
    )

    if projects.empty or "project_home" not in projects.columns:
        columns = [
            "project_home",
            "project_type",
            "organism",
            "file_source",
            "n_projects",
        ]
        return pd.DataFrame(columns=columns)

    if "project_type" not in projects.columns:
        projects["project_type"] = (
            projects["project_home"].astype(str).map(_path_dirname)
        )

    if "file_source" not in projects.columns:
        projects["file_source"] = pd.NA

    group_cols = ["project_home", "project_type", "organism", "file_source"]
    grouped = (
        projects.groupby(group_cols, dropna=False)["project"]
        .nunique()
        .reset_index(name="n_projects")
    )

    return grouped


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

_ANNOTATION_NAME_TO_EXT_HUMAN: Final[dict[str, str]] = {
    "gencode_v26": "G026",
    "gencode_v29": "G029",
    "fantom6_cat": "F006",
    "refseq": "R109",
    "ercc": "ERCC",
    "sirv": "SIRV",
}

_ANNOTATION_NAME_TO_EXT_MOUSE: Final[dict[str, str]] = {
    "gencode_v23": "M023",
}


def annotation_options(organism: str) -> dict[str, str]:
    """Return annotation options for a given organism.

    Args:
      organism: Organism name ("human" or "mouse"), case-insensitive.

    Returns:
      A new dict mapping canonical annotation names (for example,
      "gencode_v26") to recount3 annotation file extensions (for
      example, "G026").

    Raises:
      ValueError: If the organism is not recognized.
    """
    org = _normalize_organism_name(organism)
    if org == "human":
        mapping = _ANNOTATION_NAME_TO_EXT_HUMAN
    else:  # "mouse"
        mapping = _ANNOTATION_NAME_TO_EXT_MOUSE
    return dict(mapping)


def annotation_ext(organism: str, annotation: str) -> str:
    """Return the recount3 annotation extension for a given annotation.

    This helper is analogous to the R recount3::annotation_ext() function.
    It accepts either a canonical annotation name (for example,
    "gencode_v26") or a raw extension code (for example, "G026").

    Args:
      organism: Organism name ("human" or "mouse"), case-insensitive.
      annotation: Annotation name or extension code.

    Returns:
      The recount3 annotation file extension (for example, "G026").

    Raises:
      ValueError: If the organism or annotation is not recognized.
    """
    if not annotation or not annotation.strip():
        raise ValueError("annotation must be a non-empty string.")

    org = _normalize_organism_name(organism)
    name_key = annotation.strip().lower()
    ext_candidate = annotation.strip().upper()

    if org == "human":
        valid_exts = _ANN_EXT_HUMAN
        mapping = _ANNOTATION_NAME_TO_EXT_HUMAN
    else:  # "mouse"
        valid_exts = _ANN_EXT_MOUSE
        mapping = _ANNOTATION_NAME_TO_EXT_MOUSE

    if ext_candidate in valid_exts:
        return ext_candidate

    try:
        return mapping[name_key]
    except KeyError as exc:
        valid_names = sorted(mapping.keys())
        valid_ext_list = ", ".join(valid_exts)
        valid_name_list = ", ".join(valid_names)
        message = (
            f"Unknown annotation {annotation!r} for organism {organism!r}. "
            f"Valid names are: {valid_name_list}. "
            f"Valid extensions are: {valid_ext_list}."
        )
        raise ValueError(message) from exc


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
        ``annotation_extension=...``). If provided, these win.

    Returns:
      A tuple of annotation extensions.

    Raises:
      ValueError: If the organism is invalid.
    """
    org = _normalize_organism_name(organism)

    if explicit_exts:
        return tuple(str(x).strip() for x in explicit_exts)

    if annotations is None or annotations == "default":
        return ("G026",) if org == "human" else ("M023",)

    if isinstance(annotations, str):
        if annotations == "all":
            return _ANN_EXT_HUMAN if org == "human" else _ANN_EXT_MOUSE
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
    df = _strip_md_prefix(df)

    proj_keys = ("project", "study", "project_id", "project_accession")
    samp_keys = ("sample", "sample_id", "run", "run_accession", "external_id")

    proj_col = next((k for k in proj_keys if k in df.columns), None)
    if proj_col is None:
        raise ValueError("Project column not found in metadata table.")
    mask = df[proj_col].astype(str) == str(project)
    sub = df.loc[mask]

    if sub.empty:
        raise ValueError(f"Project {project!r} not found in metadata.")

    samp_col = next((k for k in samp_keys if k in sub.columns), None)
    if samp_col is None:
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
    junction_extension: Iterable[str] = ("MM",),
    include_metadata: bool = True,
    include_bigwig: bool = False,
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Enumerate all files for a project (counts, junctions, metadata, bw).

    This function composes the existing search helpers to implement a
    one-shot, project-scoped discovery routine. It adheres to recount3's
    raw file layout documented at:
    https://rna.recount.bio/docs/raw-files.html (Sections 6.2-6.4).

    Args:
      organism: "human" or "mouse".
      data_source: "sra", "gtex", or "tcga".
      project: Study identifier (for example, "SRP009615").
      genomic_units: Which expression levels to include. Defaults to both.
      annotations: "default", "all", comma-separated string, or an iterable
        of annotation file extensions (for example, ("G026", "G029")).
      junction_type: Junction type; typically "ALL".
      junction_extension: Iterable of junction artifacts to include:
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
    samp_ids = samples_for_project(
        organism=organism,
        data_source=data_source,
        project=project,
    )

    ann_exts = _resolve_annotation_exts(
        organism=organism,
        annotations=annotations,
        explicit_exts=None,
    )

    found: list[R3Resource] = []

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
            annotation_extension=ann_exts,
            strict=strict,
            deduplicate=deduplicate,
        )

        found += search_annotations(
            organism=organism,
            genomic_unit=tuple(units),
            annotation_extension=ann_exts,
            strict=strict,
            deduplicate=deduplicate,
        )

    if junction_extension:
        found += search_count_files_junctions(
            organism=organism,
            data_source=data_source,
            project=project,
            junction_type=junction_type,
            junction_extension=tuple(junction_extension),
            strict=strict,
            deduplicate=deduplicate,
        )

    if include_metadata:
        table_names = (
            "recount_project",
            "recount_qc",
            "recount_seq_qc",
            "recount_pred",
            data_source,
        )
        found += search_metadata_files(
            organism=organism,
            data_source=data_source,
            table_name=table_names,
            project=project,
            strict=strict,
            deduplicate=deduplicate,
        )

    if include_bigwig and samp_ids:
        found += search_bigwig_files(
            organism=organism,
            data_source=data_source,
            project=project,
            sample=tuple(samp_ids),
            strict=strict,
            deduplicate=deduplicate,
        )

    def _sort_key(r: R3Resource) -> tuple:
        d = r.description
        return (
            getattr(d, "resource_type", ""),
            getattr(d, "genomic_unit", ""),
            getattr(d, "annotation_extension", ""),
            getattr(d, "junction_type", ""),
            getattr(d, "junction_extension", ""),
            getattr(d, "table_name", ""),
            getattr(d, "sample", ""),
        )

    return sorted(found, key=_sort_key)
