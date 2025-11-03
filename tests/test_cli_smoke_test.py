"""Tests for `recount3 smoke-test` subcommand behavior.

Mocks:
  * `search_data_source_metadata` to return N tiny resources.
  * Resource downloads to return a path (no real network or files).

Asserts:
  * Exactly `--limit` events are emitted.
  * Exit code is 0 on success.
"""

from __future__ import annotations

import json
from typing import List

from recount3 import search as r3_search
from recount3.resource import R3Resource
from .helpers import call_cli_inprocess, make_resource


def _meta_res():
    return make_resource(resource_type="data_source_metadata", organism="human", data_source="sra")


def test_smoke_test_progress_and_limit(monkeypatch) -> None:
    """Smoke test emits one JSONL progress event per attempted resource."""
    monkeypatch.setattr(r3_search, "search_data_source_metadata", lambda **_: [_meta_res(), _meta_res()], raising=True)
    monkeypatch.setattr(R3Resource, "download", lambda self, **kw: "/dev/null/fake.tsv", raising=True)

    code, out, err = call_cli_inprocess(["smoke-test", "--limit", "2"])
    assert code == 0, err
    events = [json.loads(l) for l in out.strip().splitlines() if l.strip()]
    assert len(events) == 2
    assert all("status" in e and "url" in e for e in events)
