"""Beta command-line entrypoint for basic tasks.

  1) Build human `samplist.txt` and `projlist.txt` via
     :func:`create_sample_project_lists`.
  2) Run a small smoke test that prints example URLs and attempts downloads
     into the current directory.

The CLI is minimal to keep the library import cost low.
"""

from __future__ import annotations

import logging
import traceback
from pathlib import Path

from .descriptions import R3ResourceDescription
from .resource import R3Resource
from .search import create_sample_project_lists


def _example_descriptions() -> list[R3ResourceDescription]:
    """Construct the example descriptions used by the smoke test."""
    return [
        # 1) Annotations
        R3ResourceDescription(
            resource_type="annotations",
            organism="human",
            genomic_unit="gene",
            annotation_file_extension="G026",
        ),
        # 2) Gene counts (needs annotation_file_extension on duffel)
        R3ResourceDescription(
            resource_type="count_files_gene_or_exon",
            organism="human",
            data_source="sra",
            genomic_unit="gene",
            project="SRP107565",
            annotation_file_extension="G026",
        ),
        # 3) Junction counts
        R3ResourceDescription(
            resource_type="count_files_junctions",
            organism="human",
            data_source="sra",
            project="SRP107565",
            junction_type="ALL",
            junction_file_extension="MM",
        ),
        # 4) Per-project metadata
        R3ResourceDescription(
            resource_type="metadata_files",
            organism="human",
            data_source="sra",
            project="SRP096765",
            table_name="recount_pred",
        ),
        # 5) BigWig
        R3ResourceDescription(
            resource_type="bigwig_files",
            organism="mouse",
            data_source="sra",
            project="DRP001299",
            sample="DRR014697",
        ),
        # 6) Data sources index (duffel requires organism)
        R3ResourceDescription(
            resource_type="data_sources",
            organism="human",
        ),
        # 7) Data-source-level metadata (duffel requires organism)
        R3ResourceDescription(
            resource_type="data_source_metadata",
            organism="human",
            data_source="sra",
        ),
    ]


def test_download() -> None:
    """Smoke-test connectivity and download logic for exemplar resources.

    Prints each URL and attempts a download (to `./downloads/<type>`), catching
    and printing exceptions rather than aborting the run (parity with original).
    """
    for desc in _example_descriptions():
        try:
            res = R3Resource(desc)
            print(res.url)
            folder_name = Path("./downloads") / str(getattr(desc, "resource_type", "unknown"))
            res.download(path=str(folder_name), cache_mode="enable")
        except Exception:
            print(traceback.format_exc())


def main() -> None:
    """Run the default CLI flow.

    - Generate sample and project lists for human data.
    - Write lists to samplist.txt / projlist.txt.
    - Run example download smoke tests.
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")

    samplist, projlist = create_sample_project_lists("human")

    with open("samplist.txt", "w", encoding="utf-8") as samples_file:
        samples_file.write("\n".join(samplist))

    with open("projlist.txt", "w", encoding="utf-8") as projects_file:
        projects_file.write("\n".join(projlist))

    test_download()
