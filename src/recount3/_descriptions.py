"""Resource descriptions and URL construction.

This module implements the multi-factory :class:`R3ResourceDescription`
and dataclasses that validate parameters and construct stable duffel
URL paths.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

_VALID_ORGANISMS = {"human", "mouse"}
_VALID_DATA_SOURCES = {"sra", "gtex", "tcga"}
_VALID_GENOMIC_UNITS = {"gene", "exon"}


def _p2(s: str | None) -> str:
    """Return the last two characters of ``s`` for duffel shard directories.

    Examples:
      * None / "" -> ""
      * "A" -> "A"
      * "SRP107565" -> "65"
    """
    if not s:
        return ""
    return s[-2:]


@dataclass(slots=True)
class _R3CommonFields:
    """Uniform field schema shared by all description subclasses.

    The multi-factory returns concrete subclasses that still expose all of these
    attributes so existing code paths remain source-compatible.
    """

    resource_type: str
    organism: Optional[str] = None
    data_source: Optional[str] = None
    genomic_unit: Optional[str] = None
    project: Optional[str] = None
    sample: Optional[str] = None
    annotation_file_extension: Optional[str] = None
    junction_type: Optional[str] = None
    junction_file_extension: Optional[str] = None
    table_name: Optional[str] = None


class R3ResourceDescription:
    """Abstract base class and multi-factory for recount3 resource descriptors.

    Usage:
      >>> desc = R3ResourceDescription(resource_type="annotations", organism="human", ...)

    The base class overrides :py:meth:`__new__` to return an instance of the
    registered subclass for the given ``resource_type``. Subclasses should
    inherit from :class:`_R3CommonFields` and :class:`R3ResourceDescription`,
    and implement :py:meth:`url_path()` plus any needed validation.
    """

    _TYPE_REGISTRY: dict[str, type["R3ResourceDescription"]] = {}

    def __new__(cls, *args, **kwargs):
        """Return an instance of the appropriate concrete subclass."""
        if cls is not R3ResourceDescription:
            return super().__new__(cls)

        rtype = kwargs.get("resource_type")
        if rtype is None and args:
            rtype = args[0]

        if not rtype:
            raise KeyError("resource_type is required to construct R3ResourceDescription.")

        subcls = cls._TYPE_REGISTRY.get(str(rtype))
        if subcls is None:
            raise ValueError(f"Unsupported resource_type: {rtype!r}")

        return super().__new__(subcls)  # type: ignore

    @classmethod
    def register_type(cls, resource_type: str):
        """Decorator to register a concrete subclass in the factory."""
        def _decorator(subcls: type["R3ResourceDescription"]):
            cls._TYPE_REGISTRY[resource_type] = subcls
            setattr(subcls, "_RESOURCE_TYPE", resource_type)
            return subcls

        return _decorator

    # ---- shared validation helpers ------------------------------------

    def _require(self, *names: str) -> None:
        for n in names:
            if getattr(self, n, None) in (None, ""):
                raise KeyError(f"Missing required field: {n}")

    def _check_organism(self) -> None:
        if getattr(self, "organism", None) not in _VALID_ORGANISMS:
            raise ValueError(f"Invalid organism: {getattr(self, 'organism', None)!r}")

    def _check_data_source(self) -> None:
        if getattr(self, "data_source", None) not in _VALID_DATA_SOURCES:
            raise ValueError(f"Invalid data_source: {getattr(self, 'data_source', None)!r}")

    def _check_genomic_unit(self) -> None:
        if getattr(self, "genomic_unit", None) not in _VALID_GENOMIC_UNITS:
            raise ValueError(f"Invalid genomic_unit: {getattr(self, 'genomic_unit', None)!r}")

    # ---- abstract API --------------------------------------------------

    def url_path(self) -> str:  # pragma: no cover
        """Return the URL path (including leading slash)."""
        raise NotImplementedError


@R3ResourceDescription.register_type("annotations")
@dataclass(slots=True)
class R3Annotations(_R3CommonFields, R3ResourceDescription):
    """Descriptor for gene/exon annotation tables."""

    def __post_init__(self) -> None:
        self._require("organism", "genomic_unit", "annotation_file_extension")
        self._check_organism()
        self._check_genomic_unit()

    def url_path(self) -> str:
        return (
            f"{self.organism}/annotations/{self.genomic_unit}_sums/"
            f"{self.organism}.{self.genomic_unit}_sums.{self.annotation_file_extension}.gtf.gz"
        )


@R3ResourceDescription.register_type("count_files_gene_or_exon")
@dataclass(slots=True)
class R3GeneOrExonCounts(_R3CommonFields, R3ResourceDescription):
    """Descriptor for per-project gene/exon count matrices."""

    def __post_init__(self) -> None:
        self._require("organism", "data_source", "genomic_unit", "project", "annotation_file_extension")
        self._check_organism()
        self._check_data_source()
        self._check_genomic_unit()

    def url_path(self) -> str:
        return (
            f"{self.organism}/data_sources/{self.data_source}/{self.genomic_unit}_sums/"
            f"{_p2(self.project)}/{self.project}/"
            f"{self.data_source}.{self.genomic_unit}_sums.{self.project}.{self.annotation_file_extension}.gz"
        )


@R3ResourceDescription.register_type("count_files_junctions")
@dataclass(slots=True)
class R3JunctionCounts(_R3CommonFields, R3ResourceDescription):
    """Descriptor for per-project junction count files."""

    def __post_init__(self) -> None:
        self._require("organism", "data_source", "project", "junction_type", "junction_file_extension")
        self._check_organism()
        self._check_data_source()

    def url_path(self) -> str:
        return (
            f"{self.organism}/data_sources/{self.data_source}/junctions/"
            f"{_p2(self.project)}/{self.project}/"
            f"{self.data_source}.junctions.{self.project}.{self.junction_type}.{self.junction_file_extension}.gz"
        )


@R3ResourceDescription.register_type("metadata_files")
@dataclass(slots=True)
class R3ProjectMetadata(_R3CommonFields, R3ResourceDescription):
    """Descriptor for per-project metadata tables."""

    def __post_init__(self) -> None:
        self._require("organism", "data_source", "project", "table_name")
        self._check_organism()
        self._check_data_source()

    def url_path(self) -> str:
        return (
            f"{self.organism}/data_sources/{self.data_source}/metadata/"
            f"{_p2(self.project)}/{self.project}/"
            f"{self.data_source}.{self.table_name}.{self.project}.MD.gz"
        )


@R3ResourceDescription.register_type("bigwig_files")
@dataclass(slots=True)
class R3BigWig(_R3CommonFields, R3ResourceDescription):
    """Descriptor for per-sample BigWig coverage files (download-only)."""

    def __post_init__(self) -> None:
        self._require("organism", "data_source", "project", "sample")
        self._check_organism()
        self._check_data_source()

    def url_path(self) -> str:
        return (
            f"{self.organism}/data_sources/{self.data_source}/base_sums/"
            f"{_p2(self.project)}/{self.project}/{_p2(self.sample)}/"
            f"{self.data_source}.base_sums.{self.project}_{self.sample}.ALL.bw"
        )


@R3ResourceDescription.register_type("data_sources")
@dataclass(slots=True)
class R3DataSources(_R3CommonFields, R3ResourceDescription):
    """Descriptor for the organism-level data-source index (homes_index)."""

    def __post_init__(self) -> None:
        self._require("organism")
        self._check_organism()

    def url_path(self) -> str:
        return f"{self.organism}/homes_index"


@R3ResourceDescription.register_type("data_source_metadata")
@dataclass(slots=True)
class R3DataSourceMetadata(_R3CommonFields, R3ResourceDescription):
    """Descriptor for source-level metadata listing tables."""

    def __post_init__(self) -> None:
        self._require("organism", "data_source")
        self._check_organism()
        self._check_data_source()

    def url_path(self) -> str:
        return (
            f"{self.organism}/data_sources/{self.data_source}/metadata/"
            f"{self.data_source}.recount_project.MD.gz"
        )
