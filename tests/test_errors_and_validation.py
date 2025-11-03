"""Unit tests for parsing and validation helpers inside the CLI module."""

from __future__ import annotations

import json

import recount3.cli as cli
from .helpers import call_cli_inprocess


def test_parse_filters_rejects_bad_tokens() -> None:
    """`key=value` is required; missing '=' should cause a usage error."""
    code, _out, err = call_cli_inprocess(
        ["search", "annotations", "organism:human", "genomic_unit=gene", "annotation_file_extension=G026"]
    )
    # The main dispatcher catches ValueError and returns exit code 2.
    assert code == 2
    assert "Expected key=value" in err


def test_download_inline_malformed_json() -> None:
    """Malformed JSON should return exit code 1 with a clear error."""
    code, _out, err = call_cli_inprocess(
        ["download", "--inline", "{not-json}", "--dest", "."]
    )
    assert code == 1
    assert "Bad --inline JSON" in err
