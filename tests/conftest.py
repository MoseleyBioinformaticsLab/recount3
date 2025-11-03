"""Shared pytest fixtures and test configuration for the recount3 package.

This module centralizes fixtures used across the test suite. It intentionally
keeps tests fast, deterministic, and free of real network dependencies. The
fixtures below are designed to:

* Provide a clean runtime configuration (`fake_cfg`) that isolates cache,
  chunk-size, and timeouts to temporary directories.
* Ensure environment variables do not accidentally influence behavior
  (`env_clean`).
* Construct tiny resource descriptions (`sample_descs`) that exercise each
  resource family without requiring the public mirror.
* Prevent network access by default (`no_network`) so unit tests cannot
  inadvertently fetch from remote services.
* Offer small quality-of-life helpers (e.g., `capsys_text`) that simplify
  capturing stdout/stderr when calling the CLI in-process.
"""

from __future__ import annotations

from typing import Any

import pytest

from recount3.config import Config, default_config
from recount3._descriptions import R3ResourceDescription
from recount3.resource import R3Resource


@pytest.fixture(autouse=True)
def env_clean(monkeypatch: pytest.MonkeyPatch) -> None:
    """Clear recount3-relevant environment variables for deterministic tests.

    Pytest automatically applies this fixture to every test (autouse=True).
    Removing configuration environment variables avoids tests being affected by
    a developer's shell environment.

    Args:
      monkeypatch: Pytest monkeypatch fixture.

    Returns:
      None. The environment is modified for the lifetime of each test.
    """
    for key in (
        "RECOUNT3_URL",
        "RECOUNT3_CACHE_DIR",
        "RECOUNT3_CACHE_DISABLE",
        "RECOUNT3_HTTP_TIMEOUT",
        "RECOUNT3_MAX_RETRIES",
        "RECOUNT3_INSECURE_SSL",
        "RECOUNT3_USER_AGENT",
    ):
        monkeypatch.delenv(key, raising=False)


@pytest.fixture()
def fake_cfg(tmp_path) -> Config:
    """Create an isolated :class:`Config` for tests.

    The returned configuration points caching into a temporary directory, uses a
    small chunk-size for fast local writes, and inherits other defaults from
    :func:`default_config`.

    Args:
      tmp_path: Pytest-provided temporary directory path.

    Returns:
      A :class:`Config` instance safe for unit tests.
    """
    base = default_config()
    return Config(
        base_url=base.base_url,
        timeout=base.timeout,
        insecure_ssl=False,
        max_retries=base.max_retries,
        user_agent="recount3-tests",
        cache_dir=tmp_path / ".cache",
        cache_disabled=True,  # Always disabled for unit tests.
        chunk_size=64 * 1024,
    )


@pytest.fixture()
def sample_descs() -> dict[str, R3ResourceDescription]:
    """Return minimal, valid resource descriptions for each family.

    This fixture deliberately covers each resource “mode” accepted by the CLI
    `search` command so that tests can construct :class:`R3Resource` objects
    without hitting remote indices.

    Returns:
      A mapping from a short key to :class:`R3ResourceDescription`.
    """
    return {
        # 1) Annotations (organism + genomic_unit + annotation_file_extension).
        "annotations": R3ResourceDescription(
            resource_type="annotations",
            organism="human",
            genomic_unit="gene",
            annotation_file_extension="G026",
        ),
        # 2) Gene/exon count files.
        "gene_exon": R3ResourceDescription(
            resource_type="count_files_gene_or_exon",
            organism="human",
            data_source="sra",
            genomic_unit="gene",
            project="SRP000001",
            annotation_file_extension="G026",
        ),
        # 3) Junction counts.
        "junctions": R3ResourceDescription(
            resource_type="count_files_junctions",
            organism="human",
            data_source="sra",
            project="SRP000001",
            junction_type="ALL",
            junction_file_extension="MM",
        ),
        # 4) Metadata files (per project + table).
        "metadata": R3ResourceDescription(
            resource_type="metadata_files",
            organism="human",
            data_source="sra",
            project="SRP000002",
            table_name="recount_pred",
        ),
        # 5) BigWig one sample.
        "bigwig": R3ResourceDescription(
            resource_type="bigwig_files",
            organism="mouse",
            data_source="sra",
            project="DRP001299",
            sample="DRR014697",
        ),
        # 6) Data sources (duffel requires organism).
        "sources": R3ResourceDescription(
            resource_type="data_sources",
            organism="human",
        ),
        # 7) Data-source-level metadata (duffel requires organism).
        "source_meta": R3ResourceDescription(
            resource_type="data_source_metadata",
            organism="human",
            data_source="sra",
        ),
    }


@pytest.fixture()
def no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    """Block any attempt to download by default.

    Unit tests should not access the network. This fixture patches
    :meth:`R3Resource.download` so any un-mocked call fails loudly.
    Tests that need to simulate downloads can re-patch at a narrower scope.

    Args:
      monkeypatch: Pytest monkeypatch fixture.

    Returns:
      None. The patch is active for the lifetime of the test.
    """

    def _boom(*_args: Any, **_kwargs: Any) -> str:
        raise RuntimeError("Network disabled in tests")

    monkeypatch.setattr(R3Resource, "download", _boom, raising=True)


@pytest.fixture()
def capsys_text(capsys: pytest.CaptureFixture[str]):
    """Small helper to read captured stdout/stderr as strings.

    This wraps :func:`capsys.readouterr` with a friendlier name so call-sites
    read clearly.

    Args:
      capsys: The pytest capture fixture.

    Returns:
      A function with no arguments that returns a tuple ``(stdout, stderr)``.
    """

    def _read() -> tuple[str, str]:
        out = capsys.readouterr()
        return out.out, out.err

    return _read
