"""Project-scoped workflows and convenience APIs.

This module defines :class:`R3Project`, a thin subclass of
:class:`recount3.bundle.R3ResourceBundle` that binds the trio
(organism, data_source, project) and offers project-aware discovery
helpers plus convenience methods for common slices (counts, metadata,
bigWigs).

Typical usage example:
  from recount3.project import R3Project

  proj = R3Project.discover(
      organism="human",
      data_source="sra",
      project="SRP009615",
      include_bigwig=False,
  )
  proj.download_all(dest="downloads/")
  counts = proj.counts().stack_count_matrices(out="counts.parquet")

The public API intentionally mirrors the "discover -> manifest -> download"
workflow used elsewhere in this package; here it is simply bundled around
a project identity to reduce boilerplate for end users.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Optional

from recount3.bundle import R3ResourceBundle
from recount3.resource import R3Resource
from recount3.types import CacheMode
from recount3.search import (
    samples_for_project,
    search_project_all,
)


class R3Project(R3ResourceBundle):
    """A project-scoped bundle of recount3 resources.

    Instances bind ``organism``, ``data_source``, and ``project`` and
    provide convenience selectors for counts, metadata, and bigWigs.

    Attributes:
      organism: Organism identifier ("human" or "mouse").
      data_source: Where the data was sourced from ("sra", "gtex", "tcga").
      project: Study/project identifier (for example, "SRP009615").
    """

    def __init__(
        self,
        *,
        organism: str,
        data_source: str,
        project: str,
        resources: Optional[Iterable[R3Resource]] = None,
    ) -> None:
        """Initialize a new :class:`R3Project`.

        Args:
          organism: Organism identifier ("human" or "mouse").
          data_source: Data source ("sra", "gtex", or "tcga").
          project: Project/study accession (for example, "SRP009615").
          resources: Optional initial resources to seed the bundle.

        Raises:
          ValueError: If any identifier is empty.
        """
        super().__init__()  # Initialize base container
        if not organism or not data_source or not project:
            raise ValueError("organism, data_source, and project are required.")
        self.organism = organism
        self.data_source = data_source
        self.project = project
        if resources:
            self.extend(list(resources))

    # ------------------------------------------------------------------ #
    # Factory
    # ------------------------------------------------------------------ #

    @classmethod
    def discover(
        cls,
        *,
        organism: str,
        data_source: str,
        project: str,
        genomic_units: tuple[str, ...] = ("gene", "exon"),
        annotations: str | tuple[str, ...] = "default",
        junction_exts: tuple[str, ...] = ("MM",),
        junction_type: str = "ALL",
        include_metadata: bool = True,
        include_bigwig: bool = False,
        strict: bool = True,
        deduplicate: bool = True,
    ) -> "R3Project":
        """Validate a project, enumerate resources, and return a bundle.

        Args:
          organism: Organism ("human" or "mouse").
          data_source: Data source ("sra", "gtex", or "tcga").
          project: Study/project identifier (for example, "SRP009615").
          genomic_units: Expression feature levels to include. Defaults to
            both ("gene", "exon").
          annotations: Either "default", "all", or a tuple of specific
            annotation file extensions (for example, ("G026", "G029")).
          junction_exts: Junction artifact file extensions to include. Use
            "MM" for counts; add "RR" and/or "ID" to include coordinates and
            sample IDs.
          junction_type: Junction type; typically "ALL".
          include_metadata: Whether to include the 5 project metadata tables.
          include_bigwig: Whether to include per-sample BigWig coverage files.
          strict: If True, raise on invalid parameters; otherwise skip.
          deduplicate: If True, drop duplicate resources.

        Returns:
          An :class:`R3Project` containing discovered resources.

        Raises:
          ValueError: If the project is not found or inputs are invalid.
        """
        resources = search_project_all(
            organism=organism,
            data_source=data_source,
            project=project,
            genomic_units=genomic_units,
            annotations=annotations,
            junction_file_extension=junction_exts,
            junction_type=junction_type,
            include_metadata=include_metadata,
            include_bigwig=include_bigwig,
            strict=strict,
            deduplicate=deduplicate,
        )
        return cls(
            organism=organism,
            data_source=data_source,
            project=project,
            resources=resources,
        )

    # ------------------------------------------------------------------ #
    # Convenience accessors
    # ------------------------------------------------------------------ #

    def samples(self) -> list[str]:
        """Return the list of sample identifiers associated with the project.

        The sample IDs are mined from the data-source-level metadata tables
        and typically correspond to columns such as "external_id", "run",
        or "run_accession".

        Returns:
          Sorted list of unique sample identifiers for this project.

        Raises:
          ValueError: If the project cannot be validated.
        """
        return samples_for_project(
            organism=self.organism,
            data_source=self.data_source,
            project=self.project,
        )

    def counts(self) -> R3ResourceBundle:
        """Return a sub-bundle containing gene/exon and junction counts.

        This helper selects resources whose type is either
        "count_files_gene_or_exon" or "count_files_junctions".

        Returns:
          A new :class:`R3ResourceBundle` with count resources only.
        """
        a = self.filter(resource_type="count_files_gene_or_exon")
        b = self.filter(resource_type="count_files_junctions")
        out = R3ResourceBundle(resources=[])
        out.extend(a.resources)
        out.extend(b.resources)
        return out

    def metadata(self) -> R3ResourceBundle:
        """Return a sub-bundle containing the five project metadata tables.

        Returns:
          A new :class:`R3ResourceBundle` with metadata resources only.
        """
        return self.filter(resource_type="metadata")

    def bigwigs(self) -> R3ResourceBundle:
        """Return a sub-bundle containing BigWig coverage files.

        Returns:
          A new :class:`R3ResourceBundle` with BigWig resources only.
        """
        return self.filter(resource_type="bigwig_files")

    # ------------------------------------------------------------------ #
    # One-shot materialization
    # ------------------------------------------------------------------ #

    def download_all(
        self,
        *,
        dest: str = ".",
        overwrite: bool = False,
        cache: CacheMode = "enable",
    ) -> None:
        """Download all resources in the project bundle.

        This is a convenience wrapper around :meth:`R3Resource.download`
        for each resource in the bundle. For fine-grained control (events,
        JSON logs, or zips), prefer using the CLI ``recount3 download``.

        Args:
          dest: Destination directory or ".zip" path.
          overwrite: If True, overwrite files in directory mode.
          cache: Cache behavior: "enable", "disable", or "update".

        Raises:
          ValueError: If any resource download fails in strict mode. Failures
            when using a cache policy will raise in the underlying call.
        """
        for res in self.resources:
            res.download(path=dest, cache_mode=cache, overwrite=overwrite)
