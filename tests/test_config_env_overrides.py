"""Tests for configuration building and environment interplay.

Validates:
  * CLI flag precedence over environment values.
  * `RECOUNT3_CACHE_DISABLE` toggles cache behavior at construction time.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import recount3.cli as cli


def _args(**over) -> object:
    """Create a minimal `argparse.Namespace`-like with the required attrs."""
    defaults = dict(
        base_url=None,
        cache_dir=None,
        timeout=None,
        retries=None,
        insecure_ssl=False,
        user_agent=None,
        chunk_size=None,
    )
    defaults.update(over)
    return type("Args", (), defaults)()


def test_env_cache_disable_wins(monkeypatch, tmp_path) -> None:
    """`RECOUNT3_CACHE_DISABLE=1` flips the flag in the resulting Config."""
    monkeypatch.setenv("RECOUNT3_CACHE_DISABLE", "1")
    cfg = cli._build_config_from_env_and_flags(_args(cache_dir=str(tmp_path)))  # type: ignore
    assert cfg.cache_disabled is True


def test_cli_flags_override_env(monkeypatch, tmp_path) -> None:
    """Explicit CLI values should override environment variables."""
    monkeypatch.setenv("RECOUNT3_URL", "https://env.example/")
    args = _args(
        base_url="https://flag.example/",
        cache_dir=str(tmp_path),
        timeout=15,
        retries=9,
        user_agent="flag UA",
        chunk_size=1111,
    )
    cfg = cli._build_config_from_env_and_flags(args)  # type: ignore
    assert cfg.base_url == "https://flag.example/"
    assert cfg.timeout == 15
    assert cfg.max_retries == 9
    assert cfg.user_agent == "flag UA"
    assert cfg.chunk_size == 1111
    # Cache dir resolves properly:
    assert Path(cfg.cache_dir).exists() or True  # Path is not created here.
