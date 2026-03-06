"""Resource descriptions and duffel URL-path construction.

This module defines a small set of *resource description* dataclasses for the
recount3 duffel layout. A description is a validated, immutable-ish bundle of
parameters (organism, project, etc.) that can deterministically construct the
relative path to a resource in the duffel repository.

The main entry point is :class:`R3ResourceDescription`, which acts as a
multi-factory: instantiating ``R3ResourceDescription(resource_type=...)`` returns
an instance of the registered concrete subclass for that ``resource_type``.

Typical usage:

  desc = R3ResourceDescription(
      resource_type="count_files_gene_or_exon",
      organism="human",
      data_source="sra",
      genomic_unit="gene",
      project="SRP107565",
      annotation_extension="G026",
  )
  path = desc.url_path()

You can also pass the resource type positionally as the first argument:

  desc = R3ResourceDescription("data_sources", organism="human")

Design notes:
- The duffel layout often uses 2-character "shard" directories derived from an
  identifier (e.g., project ID) to reduce directory fanout. The helper `_p2()`
  implements that convention.
- Concrete description classes are dataclasses with ``slots=True`` to minimize
  per-instance overhead and reduce accidental attribute creation.
- ``_TYPE_REGISTRY`` is mutable by design: it is the factory registry used to
  bind resource-type strings to concrete classes. See `register_type()`.
"""

from __future__ import annotations

import dataclasses
from typing import Any, Callable, Optional, cast

VALID_ORGANISMS: frozenset[str] = frozenset({"human", "mouse"})
VALID_DATA_SOURCES: frozenset[str] = frozenset({"sra", "gtex", "tcga"})
_VALID_GENOMIC_UNITS: frozenset[str] = frozenset({"gene", "exon"})


def _p2(value: str | None) -> str:
    """Returns the last two characters of an identifier for duffel shard paths.

    Many duffel directories are sharded by a 2-character suffix to reduce the
    number of entries in a single directory. This helper returns that suffix.

    Args:
        value: Identifier string to shard (e.g., a project ID) or None.

    Returns:
        The final two characters of `value`. If `value` is None or empty,
        returns the empty string.

    Examples:
        >>> _p2(None)
        ''
        >>> _p2('')
        ''
        >>> _p2('A')
        'A'
        >>> _p2('SRP107565')
        '65'
    """
    if not value:
        return ""
    return value[-2:]


@dataclasses.dataclass(slots=True)
class _R3CommonFields:
    """Shared field schema for all resource description subclasses.

    Concrete resource description classes inherit from this dataclass to provide
    a uniform attribute surface area across all resource types. This enables
    generic code (e.g., serialization or filtering) to access common attributes
    without special-casing individual description subclasses.

    Only a subset of fields are required for a given resource type; each
    concrete subclass enforces its own required fields in `__post_init__`.

    Attributes:
        resource_type: The resource-type discriminator used by the factory.
        organism: Organism identifier (e.g., "human", "mouse").
        data_source: Source identifier (e.g., "sra", "gtex", "tcga").
        genomic_unit: Unit for feature-count matrices (e.g., "gene", "exon").
        project: Project identifier (e.g., SRP/TCGA project).
        sample: Sample/run identifier used by some resources (e.g., BigWig).
        annotation_extension: Annotation build tag (e.g., "G026").
        junction_type: Junction category selector (e.g., "ALL").
        junction_extension: Junction file format tag (e.g., "MM", "ID", "RR").
        table_name: Metadata table name (e.g., "recount_qc").
    """

    resource_type: str
    organism: Optional[str] = None
    data_source: Optional[str] = None
    genomic_unit: Optional[str] = None
    project: Optional[str] = None
    sample: Optional[str] = None
    annotation_extension: Optional[str] = None
    junction_type: Optional[str] = None
    junction_extension: Optional[str] = None
    table_name: Optional[str] = None


class R3ResourceDescription:
    """Abstract base class and multi-factory for recount3 resource descriptors.

    Instantiate this class with a `resource_type` to obtain an instance of the
    registered concrete subclass.

    Example:
        >>> desc = R3ResourceDescription(
        ...     resource_type="annotations",
        ...     organism="human",
        ...     genomic_unit="gene",
        ...     annotation_extension="G026",
        ... )
        >>> isinstance(desc, R3Annotations)
        True

    Concrete subclasses should:
    - Inherit from both `_R3CommonFields` and `R3ResourceDescription`.
    - Validate required parameters in `__post_init__`.
    - Implement `url_path()` to return the duffel-relative path.

    Attributes:
        _TYPE_REGISTRY: Mapping from resource-type string to concrete class.
            This is mutable by design to support dynamic registration.
    """

    _TYPE_REGISTRY: dict[str, type["R3ResourceDescription"]] = {}

    def __new__(
        cls,
        *args: Any,
        **kwargs: Any,
    ) -> "R3ResourceDescription":
        """Constructs and returns an instance of the appropriate subclass.

        The factory accepts the resource type either as:
        - keyword argument `resource_type=...`, or
        - the first positional argument.

        Args:
            *args: Optional positional arguments; if present, `args[0]` may
                supply the resource type.
            **kwargs: Keyword arguments used to initialize the selected
                dataclass subclass.

        Returns:
            An instance of the concrete subclass registered for the selected
            resource type.

        Raises:
            KeyError: If `resource_type` is missing or empty.
            ValueError: If `resource_type` is not registered.
        """
        if cls is not R3ResourceDescription:
            return super().__new__(cls)

        resource_type = kwargs.get("resource_type")
        if resource_type is None and args:
            resource_type = args[0]

        if not resource_type:
            raise KeyError(
                "resource_type is required to construct R3ResourceDescription."
            )

        subcls = cls._TYPE_REGISTRY.get(str(resource_type))
        if subcls is None:
            raise ValueError(f"Unsupported resource_type: {resource_type!r}")

        return cast(
            "R3ResourceDescription",
            object.__new__(subcls),
        )

    @classmethod
    def register_type(
        cls,
        resource_type: str,
    ) -> Callable[[type["R3ResourceDescription"]], type["R3ResourceDescription"]]:
        """Registers a concrete subclass for a given `resource_type`.

        This decorator binds `resource_type` to the provided subclass and also
        sets an internal `_RESOURCE_TYPE` attribute on that subclass.

        Args:
            resource_type: Resource-type string used as the factory key.

        Returns:
            A decorator that registers the decorated subclass.

        Example:
            >>> @R3ResourceDescription.register_type("annotations")
            ... @dataclasses.dataclass(slots=True)
            ... class R3Annotations(_R3CommonFields, R3ResourceDescription):
            ...     ...
        """

        def _decorator(
            subcls: type["R3ResourceDescription"],
        ) -> type["R3ResourceDescription"]:
            cls._TYPE_REGISTRY[resource_type] = subcls
            setattr(subcls, "_RESOURCE_TYPE", resource_type)
            return subcls

        return _decorator


    def _require(self, *names: str) -> None:
        """Ensures that the given attributes are present and non-empty.

        Args:
            *names: Attribute names to validate.

        Raises:
            KeyError: If any named attribute is None or the empty string.
        """
        for name in names:
            if getattr(self, name, None) in (None, ""):
                raise KeyError(f"Missing required field: {name}")

    def _check_organism(self) -> None:
        """Validates that `organism` is one of the supported values.

        Raises:
            ValueError: If `organism` is not in `VALID_ORGANISMS`.
        """
        organism = getattr(self, "organism", None)
        if organism not in VALID_ORGANISMS:
            raise ValueError(f"Invalid organism: {organism!r}")

    def _check_data_source(self) -> None:
        """Validates that `data_source` is one of the supported values.

        Raises:
            ValueError: If `data_source` is not in `VALID_DATA_SOURCES`.
        """
        data_source = getattr(self, "data_source", None)
        if data_source not in VALID_DATA_SOURCES:
            raise ValueError(f"Invalid data_source: {data_source!r}")

    def _check_genomic_unit(self) -> None:
        """Validates that `genomic_unit` is one of the supported values.

        Raises:
            ValueError: If `genomic_unit` is not in `_VALID_GENOMIC_UNITS`.
        """
        genomic_unit = getattr(self, "genomic_unit", None)
        if genomic_unit not in _VALID_GENOMIC_UNITS:
            raise ValueError(f"Invalid genomic_unit: {genomic_unit!r}")


    def url_path(self) -> str:
        """Return the duffel-relative URL path for this resource.

        `R3ResourceDescription` is a multi-factory and abstract interface: calling
        `R3ResourceDescription(resource_type=...)` returns an instance of a concrete
        registered subclass (e.g., `R3Annotations`, `R3GeneOrExonCounts`). Those
        subclasses implement `url_path()` to construct a deterministic path within
        the duffel repository layout.

        This base-class method is therefore not expected to be called directly; it
        exists to document the interface and to provide a consistent method name
        across all description types.

        Implementations must return a path without a leading slash. Callers should
        join this value onto the configured base URL (and add any separators as
        needed).

        Returns:
            A duffel-relative path string, no leading slash.

        Raises:
            NotImplementedError: Always, in the base class. Concrete subclasses
                override this method.
        """
        raise NotImplementedError


@R3ResourceDescription.register_type("annotations")
@dataclasses.dataclass(slots=True)
class R3Annotations(_R3CommonFields, R3ResourceDescription):
    """Descriptor for annotation GTF files.

    Required fields:
      - organism
      - genomic_unit
      - annotation_extension

    Duffel layout:
      {organism}/annotations/{genomic_unit}_sums/
        {organism}.{genomic_unit}_sums.{annotation_extension}.gtf.gz
    """

    def __post_init__(self) -> None:
        self._require("organism", "genomic_unit", "annotation_extension")
        self._check_organism()
        self._check_genomic_unit()

    def url_path(self) -> str:
        return (
            f"{self.organism}/annotations/{self.genomic_unit}_sums/"
            f"{self.organism}.{self.genomic_unit}_sums."
            f"{self.annotation_extension}.gtf.gz"
        )


@R3ResourceDescription.register_type("count_files_gene_or_exon")
@dataclasses.dataclass(slots=True)
class R3GeneOrExonCounts(_R3CommonFields, R3ResourceDescription):
    """Descriptor for per-project gene/exon count matrices.

    Required fields:
      - organism
      - data_source
      - genomic_unit
      - project
      - annotation_extension

    Duffel layout:
      {organism}/data_sources/{data_source}/{genomic_unit}_sums/
        {p2(project)}/{project}/
        {data_source}.{genomic_unit}_sums.{project}.{annotation_extension}.gz
    """

    def __post_init__(self) -> None:
        self._require(
            "organism",
            "data_source",
            "genomic_unit",
            "project",
            "annotation_extension",
        )
        self._check_organism()
        self._check_data_source()
        self._check_genomic_unit()

    def url_path(self) -> str:
        return (
            f"{self.organism}/data_sources/{self.data_source}/"
            f"{self.genomic_unit}_sums/{_p2(self.project)}/{self.project}/"
            f"{self.data_source}.{self.genomic_unit}_sums.{self.project}."
            f"{self.annotation_extension}.gz"
        )


@R3ResourceDescription.register_type("count_files_junctions")
@dataclasses.dataclass(slots=True)
class R3JunctionCounts(_R3CommonFields, R3ResourceDescription):
    """Descriptor for per-project junction count files.

    Required fields:
      - organism
      - data_source
      - project
      - junction_type
      - junction_extension

    Duffel layout:
      {organism}/data_sources/{data_source}/junctions/
        {p2(project)}/{project}/
        {data_source}.junctions.{project}.{junction_type}.{junction_extension}.gz
    """

    def __post_init__(self) -> None:
        self._require(
            "organism",
            "data_source",
            "project",
            "junction_type",
            "junction_extension",
        )
        self._check_organism()
        self._check_data_source()

    def url_path(self) -> str:
        return (
            f"{self.organism}/data_sources/{self.data_source}/junctions/"
            f"{_p2(self.project)}/{self.project}/"
            f"{self.data_source}.junctions.{self.project}."
            f"{self.junction_type}.{self.junction_extension}.gz"
        )


@R3ResourceDescription.register_type("metadata_files")
@dataclasses.dataclass(slots=True)
class R3ProjectMetadata(_R3CommonFields, R3ResourceDescription):
    """Descriptor for per-project metadata tables.

    Required fields:
      - organism
      - data_source
      - project
      - table_name

    Duffel layout:
      {organism}/data_sources/{data_source}/metadata/
        {p2(project)}/{project}/
        {data_source}.{table_name}.{project}.MD.gz
    """

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
@dataclasses.dataclass(slots=True)
class R3BigWig(_R3CommonFields, R3ResourceDescription):
    """Descriptor for per-sample BigWig coverage files.

    Required fields:
      - organism
      - data_source
      - project
      - sample

    Duffel layout:
      {organism}/data_sources/{data_source}/base_sums/
        {p2(project)}/{project}/{p2(sample)}/
        {data_source}.base_sums.{project}_{sample}.ALL.bw
    """

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
@dataclasses.dataclass(slots=True)
class R3DataSources(_R3CommonFields, R3ResourceDescription):
    """Descriptor for the organism-level data-source index (homes_index).

    Required fields:
      - organism

    Duffel layout:
      {organism}/homes_index
    """

    def __post_init__(self) -> None:
        self._require("organism")
        self._check_organism()

    def url_path(self) -> str:
        return f"{self.organism}/homes_index"


@R3ResourceDescription.register_type("data_source_metadata")
@dataclasses.dataclass(slots=True)
class R3DataSourceMetadata(_R3CommonFields, R3ResourceDescription):
    """Descriptor for source-level metadata listings.

    Required fields:
      - organism
      - data_source

    Duffel layout:
      {organism}/data_sources/{data_source}/metadata/
        {data_source}.recount_project.MD.gz
    """

    def __post_init__(self) -> None:
        self._require("organism", "data_source")
        self._check_organism()
        self._check_data_source()

    def url_path(self) -> str:
        return (
            f"{self.organism}/data_sources/{self.data_source}/metadata/"
            f"{self.data_source}.recount_project.MD.gz"
        )
