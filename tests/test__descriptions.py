from __future__ import annotations

import dataclasses

import pytest

import recount3._descriptions


COMMON_FIELDS = recount3._descriptions._R3CommonFields  # pylint: disable=protected-access


@pytest.mark.parametrize(
    "value,expected",
    [
        (None, ""),
        ("", ""),
        ("A", "A"),
        ("SR", "SR"),
        ("SRP107565", "65"),
    ],
)
def test_p2(value: str | None, expected: str) -> None:
    descriptions = recount3._descriptions
    result = descriptions._p2(value)  # pylint: disable=protected-access
    assert result == expected


def test_factory_requires_resource_type() -> None:
    with pytest.raises(KeyError, match=r"resource_type is required"):
        recount3._descriptions.R3ResourceDescription()

    with pytest.raises(KeyError, match=r"resource_type is required"):
        recount3._descriptions.R3ResourceDescription(resource_type="")

    with pytest.raises(KeyError, match=r"resource_type is required"):
        recount3._descriptions.R3ResourceDescription("")


def test_factory_rejects_unsupported_resource_type() -> None:
    with pytest.raises(ValueError, match=r"Unsupported resource_type"):
        recount3._descriptions.R3ResourceDescription(resource_type="nope")

    with pytest.raises(ValueError, match=r"Unsupported resource_type"):
        recount3._descriptions.R3ResourceDescription(resource_type=123)


def test_registry_contains_expected_types() -> None:
    registry = recount3._descriptions.R3ResourceDescription._TYPE_REGISTRY
    expected = {
        "annotations",
        "count_files_gene_or_exon",
        "count_files_junctions",
        "metadata_files",
        "bigwig_files",
        "data_sources",
        "data_source_metadata",
    }
    assert set(registry.keys()) == expected

    for resource_type, subcls in registry.items():
        assert subcls._RESOURCE_TYPE == resource_type  # pylint: disable=protected-access


@pytest.mark.parametrize(
    "resource_type,kwargs,expected_class",
    [
        (
            "annotations",
            {
                "organism": "human",
                "genomic_unit": "gene",
                "annotation_extension": "G026",
            },
            recount3._descriptions.R3Annotations,
        ),
        (
            "count_files_gene_or_exon",
            {
                "organism": "human",
                "data_source": "sra",
                "genomic_unit": "exon",
                "project": "SRP107565",
                "annotation_extension": "G026",
            },
            recount3._descriptions.R3GeneOrExonCounts,
        ),
        (
            "count_files_junctions",
            {
                "organism": "human",
                "data_source": "sra",
                "project": "SRP000001",
                "junction_type": "ALL",
                "junction_extension": "MM",
            },
            recount3._descriptions.R3JunctionCounts,
        ),
        (
            "metadata_files",
            {
                "organism": "mouse",
                "data_source": "gtex",
                "project": "PRJ12345",
                "table_name": "recount_qc",
            },
            recount3._descriptions.R3ProjectMetadata,
        ),
        (
            "bigwig_files",
            {
                "organism": "human",
                "data_source": "sra",
                "project": "SRP000001",
                "sample": "SAMP0001",
            },
            recount3._descriptions.R3BigWig,
        ),
        (
            "data_sources",
            {
                "organism": "human",
            },
            recount3._descriptions.R3DataSources,
        ),
        (
            "data_source_metadata",
            {
                "organism": "human",
                "data_source": "sra",
            },
            recount3._descriptions.R3DataSourceMetadata,
        ),
    ],
)
def test_factory_dispatches_to_correct_type(
    resource_type: str,
    kwargs: dict[str, str],
    expected_class: type[recount3._descriptions.R3ResourceDescription],
) -> None:
    desc = recount3._descriptions.R3ResourceDescription(
        resource_type=resource_type,
        **kwargs,
    )
    assert isinstance(desc, expected_class)

    assert desc.resource_type == resource_type


def test_factory_accepts_positional_resource_type() -> None:
    desc = recount3._descriptions.R3ResourceDescription(
        "data_sources",
        organism="human",
    )
    assert isinstance(desc, recount3._descriptions.R3DataSources)
    assert desc.url_path() == "human/homes_index"


def test_direct_subclass_construction_bypasses_factory() -> None:
    desc = recount3._descriptions.R3Annotations(
        resource_type="annotations",
        organism="human",
        genomic_unit="gene",
        annotation_extension="G026",
    )
    assert isinstance(desc, recount3._descriptions.R3Annotations)


def test_base_url_path_raises_not_implemented() -> None:
    desc = recount3._descriptions.R3Annotations(
        resource_type="annotations",
        organism="human",
        genomic_unit="gene",
        annotation_extension="G026",
    )
    with pytest.raises(NotImplementedError):
        recount3._descriptions.R3ResourceDescription.url_path(desc)  # pylint: disable=protected-access


@pytest.mark.parametrize(
    "resource_type,kwargs,missing_field",
    [
        (
            "annotations",
            {
                "organism": "human",
                "genomic_unit": "gene",
                "annotation_extension": "",
            },
            "annotation_extension",
        ),
        (
            "count_files_gene_or_exon",
            {
                "organism": "human",
                "data_source": "sra",
                "genomic_unit": "gene",
                "project": "SRP000001",
            },
            "annotation_extension",
        ),
        (
            "count_files_junctions",
            {
                "organism": "human",
                "data_source": "sra",
                "project": "SRP000001",
                "junction_type": "ALL",
            },
            "junction_extension",
        ),
        (
            "metadata_files",
            {
                "organism": "human",
                "data_source": "sra",
                "project": "",
                "table_name": "recount_qc",
            },
            "project",
        ),
        (
            "bigwig_files",
            {
                "organism": "human",
                "data_source": "sra",
                "project": "SRP000001",
                "sample": "",
            },
            "sample",
        ),
        (
            "data_sources",
            {},
            "organism",
        ),
        (
            "data_source_metadata",
            {
                "organism": "human",
            },
            "data_source",
        ),
    ],
)
def test_required_fields_enforced(
    resource_type: str,
    kwargs: dict[str, str],
    missing_field: str,
) -> None:
    pattern = rf"Missing required field: {missing_field}"
    with pytest.raises(KeyError, match=pattern):
        recount3._descriptions.R3ResourceDescription(
            resource_type=resource_type,
            **kwargs,
        )


@pytest.mark.parametrize(
    "resource_type,kwargs,expected_error",
    [
        (
            "annotations",
            {
                "organism": "dog",
                "genomic_unit": "gene",
                "annotation_extension": "G026",
            },
            r"Invalid organism",
        ),
        (
            "count_files_gene_or_exon",
            {
                "organism": "human",
                "data_source": "ena",
                "genomic_unit": "gene",
                "project": "SRP000001",
                "annotation_extension": "G026",
            },
            r"Invalid data_source",
        ),
        (
            "count_files_gene_or_exon",
            {
                "organism": "human",
                "data_source": "sra",
                "genomic_unit": "transcript",
                "project": "SRP000001",
                "annotation_extension": "G026",
            },
            r"Invalid genomic_unit",
        ),
        (
            "data_source_metadata",
            {
                "organism": "mouse",
                "data_source": "nope",
            },
            r"Invalid data_source",
        ),
    ],
)
def test_validation_rejects_invalid_enum_values(
    resource_type: str,
    kwargs: dict[str, str],
    expected_error: str,
) -> None:
    with pytest.raises(ValueError, match=expected_error):
        recount3._descriptions.R3ResourceDescription(
            resource_type=resource_type,
            **kwargs,
        )


@pytest.mark.parametrize(
    "resource_type,kwargs,expected_path",
    [
        (
            "annotations",
            {
                "organism": "human",
                "genomic_unit": "gene",
                "annotation_extension": "G026",
            },
            "human/annotations/gene_sums/human.gene_sums.G026.gtf.gz",
        ),
        (
            "count_files_gene_or_exon",
            {
                "organism": "human",
                "data_source": "sra",
                "genomic_unit": "exon",
                "project": "SRP107565",
                "annotation_extension": "G026",
            },
            "human/data_sources/sra/exon_sums/65/SRP107565/"
            "sra.exon_sums.SRP107565.G026.gz",
        ),
        (
            "count_files_junctions",
            {
                "organism": "human",
                "data_source": "sra",
                "project": "SRP000001",
                "junction_type": "ALL",
                "junction_extension": "MM",
            },
            "human/data_sources/sra/junctions/01/SRP000001/"
            "sra.junctions.SRP000001.ALL.MM.gz",
        ),
        (
            "metadata_files",
            {
                "organism": "mouse",
                "data_source": "gtex",
                "project": "PRJ12345",
                "table_name": "recount_qc",
            },
            "mouse/data_sources/gtex/metadata/45/PRJ12345/"
            "gtex.recount_qc.PRJ12345.MD.gz",
        ),
        (
            "bigwig_files",
            {
                "organism": "human",
                "data_source": "sra",
                "project": "SRP000001",
                "sample": "SAMP0001",
            },
            "human/data_sources/sra/base_sums/01/SRP000001/01/"
            "sra.base_sums.SRP000001_SAMP0001.ALL.bw",
        ),
        (
            "bigwig_files",
            {
                "organism": "human",
                "data_source": "gtex",
                "project": "BLADDER",
                "sample": "GTEX-11DXZ-0526-SM-5EQRP",
            },
            "human/data_sources/gtex/base_sums/ER/BLADDER/EQ/"
            "gtex.base_sums.BLADDER_GTEX-11DXZ-0526-SM-5EQRP.ALL.bw",
        ),
        (
            "data_sources",
            {
                "organism": "human",
            },
            "human/homes_index",
        ),
        (
            "data_source_metadata",
            {
                "organism": "human",
                "data_source": "sra",
            },
            "human/data_sources/sra/metadata/sra.recount_project.MD.gz",
        ),
    ],
)
def test_url_path_construction(
    resource_type: str,
    kwargs: dict[str, str],
    expected_path: str,
) -> None:
    desc = recount3._descriptions.R3ResourceDescription(
        resource_type=resource_type,
        **kwargs,
    )
    assert desc.url_path() == expected_path


def test_common_fields_is_dataclass_with_slots() -> None:
    assert dataclasses.is_dataclass(COMMON_FIELDS)
    assert hasattr(COMMON_FIELDS, "__slots__")

    field_names = [f.name for f in dataclasses.fields(COMMON_FIELDS)]
    assert field_names[0] == "resource_type"


@pytest.mark.parametrize(
    "sample,data_source,expected",
    [
        (None, "sra", ""),
        ("", "sra", ""),
        (None, None, ""),
        ("SRR001", "sra", "01"),
        ("SAMP0001", "sra", "01"),
        ("AB", "sra", "AB"),
        ("A", "sra", "A"),
        ("GTEX-11DXZ-0526-SM-5EQRP", "gtex", "EQ"),
        ("AB", "gtex", "AB"),
        ("A", "gtex", "A"),
        ("SRR001", None, "01"),
        ("GTEX-SAMPLE-XYZ1", "GTEX", "XY"),
    ],
)
def test_bigwig_sample_shard(
    sample: str | None,
    data_source: str | None,
    expected: str,
) -> None:
    result = recount3._descriptions._bigwig_sample_shard(  # pylint: disable=protected-access
        sample, data_source
    )
    assert result == expected