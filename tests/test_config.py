# Copyright (c) 2026, Alexander A. Alsalihi, Robert M. Flight,
# Hunter N.B. Moseley. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * All advertising materials mentioning features or use of this software must
#   display the following acknowledgement: This product includes software
#   developed by the copyright holder.
# * Neither the name of the copyright holder nor the names of its contributors
#   may be used to endorse or promote products derived from this software without
#   specific prior written permission.
# * If the source code is used in a published work, then proper citation of the
#   source code must be included with the published work.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from __future__ import annotations

import os
from pathlib import Path

import pytest

from recount3.config import (
    Config,
    _DEFAULT_CHUNK_SIZE,
    default_config,
    recount3_cache,
    recount3_cache_files,
    recount3_cache_rm,
)


_BASE_URL = "https://example.org/recount3/"


def _minimal_cfg(tmp_path: Path) -> Config:
    """Return a minimal :class:`Config` with a temporary cache directory."""
    return Config(
        base_url=_BASE_URL,
        timeout=30,
        insecure_ssl=False,
        max_retries=3,
        user_agent="test-agent/1.0",
        cache_dir=tmp_path / "cache",
        cache_disabled=False,
    )


class TestConfig:
    """Unit tests for the :class:`Config` frozen dataclass."""

    def test_construction_with_all_fields(self, tmp_path: Path) -> None:
        cfg = Config(
            base_url=_BASE_URL,
            timeout=60,
            insecure_ssl=True,
            max_retries=5,
            user_agent="ua/2.0",
            cache_dir=tmp_path,
            cache_disabled=True,
            chunk_size=512,
        )
        assert cfg.base_url == _BASE_URL
        assert cfg.timeout == 60
        assert cfg.insecure_ssl is True
        assert cfg.max_retries == 5
        assert cfg.user_agent == "ua/2.0"
        assert cfg.cache_dir == tmp_path
        assert cfg.cache_disabled is True
        assert cfg.chunk_size == 512

    def test_chunk_size_default(self, tmp_path: Path) -> None:
        """chunk_size should default to _DEFAULT_CHUNK_SIZE (1 MiB)."""
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua/1.0",
            cache_dir=tmp_path,
            cache_disabled=False,
        )
        assert cfg.chunk_size == _DEFAULT_CHUNK_SIZE
        assert _DEFAULT_CHUNK_SIZE == 1024 * 1024

    def test_frozen_raises_on_mutation(self, tmp_path: Path) -> None:
        """Mutating a frozen dataclass field must raise FrozenInstanceError."""
        cfg = _minimal_cfg(tmp_path)
        with pytest.raises(Exception):
            cfg.timeout = 999  # type: ignore[misc]

    def test_equality(self, tmp_path: Path) -> None:
        """Two Config objects with the same field values must compare equal."""
        cfg_a = _minimal_cfg(tmp_path)
        cfg_b = _minimal_cfg(tmp_path)
        assert cfg_a == cfg_b

    def test_inequality(self, tmp_path: Path) -> None:
        """Two Config objects that differ in any field must compare unequal."""
        cfg_a = _minimal_cfg(tmp_path)
        cfg_b = Config(
            base_url=_BASE_URL,
            timeout=999,
            insecure_ssl=False,
            max_retries=3,
            user_agent="test-agent/1.0",
            cache_dir=tmp_path / "cache",
            cache_disabled=False,
        )
        assert cfg_a != cfg_b

    def test_hashable(self, tmp_path: Path) -> None:
        """Frozen dataclasses are hashable; Config must be usable as a dict key."""
        cfg = _minimal_cfg(tmp_path)
        mapping = {cfg: "value"}
        assert mapping[cfg] == "value"

    def test_repr_contains_field_values(self, tmp_path: Path) -> None:
        """repr() should mention the class name and at least one field value."""
        cfg = _minimal_cfg(tmp_path)
        r = repr(cfg)
        assert "Config" in r
        assert str(cfg.timeout) in r

    def test_slots_no_dict(self, tmp_path: Path) -> None:
        """A slotted dataclass must not have __dict__."""
        cfg = _minimal_cfg(tmp_path)
        assert not hasattr(cfg, "__dict__")


class TestDefaultConfig:
    def test_url_absent_uses_builtin_default(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """When RECOUNT3_URL is absent, the hard-coded duffel URL is used."""
        monkeypatch.delenv("RECOUNT3_URL", raising=False)
        cfg = default_config()
        assert cfg.base_url == "http://duffel.rail.bio/recount3/"

    def test_url_present_without_trailing_slash(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A URL set without a trailing slash must have one appended."""
        monkeypatch.setenv("RECOUNT3_URL", "https://mirror.example.com/r3")
        cfg = default_config()
        assert cfg.base_url == "https://mirror.example.com/r3/"

    def test_url_present_with_trailing_slash(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A URL already ending with '/' must not gain a second slash."""
        monkeypatch.setenv("RECOUNT3_URL", "https://mirror.example.com/r3/")
        cfg = default_config()
        assert cfg.base_url == "https://mirror.example.com/r3/"

    def test_url_with_multiple_trailing_slashes(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Multiple trailing slashes should be stripped to exactly one."""
        monkeypatch.setenv("RECOUNT3_URL", "https://mirror.example.com/r3///")
        cfg = default_config()
        assert cfg.base_url == "https://mirror.example.com/r3/"

    def test_cache_dir_absent_uses_home_default(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """When RECOUNT3_CACHE_DIR is absent, a ~/.cache/recount3/files path is used."""
        monkeypatch.delenv("RECOUNT3_CACHE_DIR", raising=False)
        cfg = default_config()
        expected = Path(os.path.expanduser("~")) / ".cache" / "recount3" / "files"
        assert cfg.cache_dir == expected

    def test_cache_dir_present(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """When RECOUNT3_CACHE_DIR is set, that path is used."""
        monkeypatch.setenv("RECOUNT3_CACHE_DIR", str(tmp_path / "my_cache"))
        cfg = default_config()
        assert cfg.cache_dir == tmp_path / "my_cache"

    def test_timeout_absent_defaults_to_60(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("RECOUNT3_HTTP_TIMEOUT", raising=False)
        cfg = default_config()
        assert cfg.timeout == 60

    def test_timeout_present(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("RECOUNT3_HTTP_TIMEOUT", "120")
        cfg = default_config()
        assert cfg.timeout == 120

    def test_timeout_is_int(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """timeout must be an int, not a string."""
        monkeypatch.setenv("RECOUNT3_HTTP_TIMEOUT", "45")
        cfg = default_config()
        assert isinstance(cfg.timeout, int)
        assert cfg.timeout == 45

    def test_insecure_ssl_absent_is_false(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("RECOUNT3_INSECURE_SSL", raising=False)
        cfg = default_config()
        assert cfg.insecure_ssl is False

    def test_insecure_ssl_set_to_1_is_true(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("RECOUNT3_INSECURE_SSL", "1")
        cfg = default_config()
        assert cfg.insecure_ssl is True

    def test_insecure_ssl_set_to_non_1_is_false(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Any value other than "1" (e.g., "true", "yes", "0") must be False."""
        for value in ("0", "true", "yes", "True", ""):
            monkeypatch.setenv("RECOUNT3_INSECURE_SSL", value)
            cfg = default_config()
            assert cfg.insecure_ssl is False, f"Expected False for {value!r}"

    def test_max_retries_absent_defaults_to_3(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("RECOUNT3_MAX_RETRIES", raising=False)
        cfg = default_config()
        assert cfg.max_retries == 3

    def test_max_retries_present(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("RECOUNT3_MAX_RETRIES", "10")
        cfg = default_config()
        assert cfg.max_retries == 10

    def test_max_retries_is_int(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("RECOUNT3_MAX_RETRIES", "7")
        cfg = default_config()
        assert isinstance(cfg.max_retries, int)

    def test_user_agent_absent_uses_package_default(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv("RECOUNT3_USER_AGENT", raising=False)
        cfg = default_config()
        assert "recount3" in cfg.user_agent.lower()
        assert cfg.user_agent.startswith("recount3-python/")

    def test_user_agent_present_non_empty(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("RECOUNT3_USER_AGENT", "MyApp/3.1")
        cfg = default_config()
        assert cfg.user_agent == "MyApp/3.1"

    def test_user_agent_empty_string_falls_back_to_default(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """An empty-string RECOUNT3_USER_AGENT must fall back to the package default.

        The expression ``os.environ.get("RECOUNT3_USER_AGENT") or <default>``
        treats the empty string as falsy, so the package default is used.
        """
        monkeypatch.setenv("RECOUNT3_USER_AGENT", "")
        cfg = default_config()
        assert cfg.user_agent.startswith("recount3-python/")

    def test_cache_disable_absent_is_false(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("RECOUNT3_CACHE_DISABLE", raising=False)
        cfg = default_config()
        assert cfg.cache_disabled is False

    def test_cache_disable_set_to_1_is_true(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("RECOUNT3_CACHE_DISABLE", "1")
        cfg = default_config()
        assert cfg.cache_disabled is True

    def test_cache_disable_set_to_non_1_is_false(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        for value in ("0", "true", "yes", ""):
            monkeypatch.setenv("RECOUNT3_CACHE_DISABLE", value)
            cfg = default_config()
            assert cfg.cache_disabled is False, f"Expected False for {value!r}"

    def test_chunk_size_absent_uses_constant_default(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv("RECOUNT3_CHUNK_SIZE", raising=False)
        cfg = default_config()
        assert cfg.chunk_size == _DEFAULT_CHUNK_SIZE

    def test_chunk_size_present(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("RECOUNT3_CHUNK_SIZE", "65536")
        cfg = default_config()
        assert cfg.chunk_size == 65536

    def test_chunk_size_is_int(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("RECOUNT3_CHUNK_SIZE", "32768")
        cfg = default_config()
        assert isinstance(cfg.chunk_size, int)

    def test_all_env_vars_set(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """Setting every env var should produce a fully customised Config."""
        monkeypatch.setenv("RECOUNT3_URL", "https://custom.example.com/r3")
        monkeypatch.setenv("RECOUNT3_CACHE_DIR", str(tmp_path / "custom"))
        monkeypatch.setenv("RECOUNT3_HTTP_TIMEOUT", "99")
        monkeypatch.setenv("RECOUNT3_INSECURE_SSL", "1")
        monkeypatch.setenv("RECOUNT3_MAX_RETRIES", "7")
        monkeypatch.setenv("RECOUNT3_USER_AGENT", "CustomBot/9")
        monkeypatch.setenv("RECOUNT3_CACHE_DISABLE", "1")
        monkeypatch.setenv("RECOUNT3_CHUNK_SIZE", "8192")

        cfg = default_config()

        assert cfg.base_url == "https://custom.example.com/r3/"
        assert cfg.cache_dir == tmp_path / "custom"
        assert cfg.timeout == 99
        assert cfg.insecure_ssl is True
        assert cfg.max_retries == 7
        assert cfg.user_agent == "CustomBot/9"
        assert cfg.cache_disabled is True
        assert cfg.chunk_size == 8192

    def test_returns_config_instance(self) -> None:
        cfg = default_config()
        assert isinstance(cfg, Config)


class TestRecount3Cache:
    def test_with_explicit_config_returns_cache_dir(self, tmp_path: Path) -> None:
        cfg = _minimal_cfg(tmp_path)
        result = recount3_cache(config=cfg)
        assert result == cfg.cache_dir
        assert result.is_dir()

    def test_with_no_config_uses_default(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """Passing config=None falls back to default_config()."""
        monkeypatch.setenv("RECOUNT3_CACHE_DIR", str(tmp_path / "default_cache"))
        result = recount3_cache(config=None)
        assert result == tmp_path / "default_cache"
        assert result.is_dir()

    def test_creates_cache_dir_if_missing(self, tmp_path: Path) -> None:
        """The cache directory is created if it does not exist."""
        cache = tmp_path / "brand_new_cache"
        assert not cache.exists()
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        recount3_cache(config=cfg)
        assert cache.is_dir()

    def test_idempotent_if_cache_dir_exists(self, tmp_path: Path) -> None:
        """Calling recount3_cache when the dir already exists must not raise."""
        cfg = _minimal_cfg(tmp_path)
        cfg.cache_dir.mkdir(parents=True, exist_ok=True)
        result = recount3_cache(config=cfg)
        assert result.is_dir()


class TestRecount3CacheFiles:
    """Tests for :func:`recount3_cache_files`."""

    def test_nonexistent_cache_dir_returns_empty(self, tmp_path: Path) -> None:
        """When the cache directory does not exist, an empty list is returned."""
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=tmp_path / "no_such_dir",
            cache_disabled=False,
        )
        result = recount3_cache_files(config=cfg)
        assert result == []

    def test_cache_path_is_file_not_dir_returns_empty(self, tmp_path: Path) -> None:
        """When cache_dir is a regular file (not a directory), return []."""
        file_path = tmp_path / "i_am_a_file"
        file_path.write_text("data")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=file_path,
            cache_disabled=False,
        )
        result = recount3_cache_files(config=cfg)
        assert result == []

    def test_empty_cache_dir_returns_empty(self, tmp_path: Path) -> None:
        """An existing but empty cache directory returns an empty list."""
        cache = tmp_path / "cache"
        cache.mkdir()
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        result = recount3_cache_files(config=cfg)
        assert result == []

    def test_returns_files_in_sorted_order(self, tmp_path: Path) -> None:
        """Files are returned in sorted order regardless of insertion order."""
        cache = tmp_path / "cache"
        cache.mkdir()
        (cache / "zzz.gz").write_bytes(b"z")
        (cache / "aaa.gz").write_bytes(b"a")
        (cache / "mmm.gz").write_bytes(b"m")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        result = recount3_cache_files(config=cfg)
        assert result == sorted(result, key=str)
        assert len(result) == 3

    def test_subdirectory_files_are_included(self, tmp_path: Path) -> None:
        """rglob is used so files in subdirectories are also returned."""
        cache = tmp_path / "cache"
        sub = cache / "sub" / "deep"
        sub.mkdir(parents=True)
        (sub / "nested.txt").write_text("content")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        result = recount3_cache_files(config=cfg)
        assert len(result) == 1
        assert result[0].name == "nested.txt"

    def test_directories_are_not_returned(self, tmp_path: Path) -> None:
        """Directories matched by rglob are excluded; only files are returned."""
        cache = tmp_path / "cache"
        sub = cache / "subdir"
        sub.mkdir(parents=True)
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        result = recount3_cache_files(config=cfg)
        assert result == []

    def test_pattern_filters_matching_files(self, tmp_path: Path) -> None:
        """When pattern is provided, only matching files are returned."""
        cache = tmp_path / "cache"
        cache.mkdir()
        (cache / "data.tsv.gz").write_bytes(b"tsv")
        (cache / "data.parquet").write_bytes(b"parquet")
        (cache / "readme.txt").write_bytes(b"readme")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        result = recount3_cache_files(config=cfg, pattern="*.gz")
        assert len(result) == 1
        assert result[0].name == "data.tsv.gz"

    def test_pattern_none_returns_all_files(self, tmp_path: Path) -> None:
        """When pattern=None, all files are returned (glob_pattern defaults to '*')."""
        cache = tmp_path / "cache"
        cache.mkdir()
        (cache / "a.txt").write_bytes(b"a")
        (cache / "b.txt").write_bytes(b"b")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        result = recount3_cache_files(config=cfg, pattern=None)
        assert len(result) == 2

    def test_pattern_no_match_returns_empty(self, tmp_path: Path) -> None:
        """A pattern that matches nothing returns an empty list."""
        cache = tmp_path / "cache"
        cache.mkdir()
        (cache / "file.txt").write_bytes(b"data")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        result = recount3_cache_files(config=cfg, pattern="*.parquet")
        assert result == []

    def test_with_no_config_uses_default_config(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """config=None must fall back to default_config()."""
        cache = tmp_path / "dc_cache"
        cache.mkdir()
        (cache / "x.txt").write_bytes(b"x")
        monkeypatch.setenv("RECOUNT3_CACHE_DIR", str(cache))
        result = recount3_cache_files(config=None)
        assert len(result) == 1

    def test_returns_list_of_path_objects(self, tmp_path: Path) -> None:
        """Returned entries are :class:`pathlib.Path` instances."""
        cache = tmp_path / "cache"
        cache.mkdir()
        (cache / "file.bin").write_bytes(b"data")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        result = recount3_cache_files(config=cfg)
        assert all(isinstance(p, Path) for p in result)


class TestRecount3CacheRm:
    """Tests for :func:`recount3_cache_rm`."""

    def test_nonexistent_cache_dir_returns_empty(self, tmp_path: Path) -> None:
        """When the cache directory does not exist, return []."""
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=tmp_path / "no_such_dir",
            cache_disabled=False,
        )
        result = recount3_cache_rm(config=cfg)
        assert result == []

    def test_cache_path_is_file_not_dir_returns_empty(self, tmp_path: Path) -> None:
        """When cache_dir is a file (not a directory), return []."""
        file_path = tmp_path / "i_am_a_file"
        file_path.write_text("data")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=file_path,
            cache_disabled=False,
        )
        result = recount3_cache_rm(config=cfg)
        assert result == []

    def test_empty_cache_dir_returns_empty(self, tmp_path: Path) -> None:
        """An empty cache directory returns [] and raises no errors."""
        cache = tmp_path / "cache"
        cache.mkdir()
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        result = recount3_cache_rm(config=cfg)
        assert result == []

    def test_predicate_none_removes_all_files(self, tmp_path: Path) -> None:
        """With predicate=None every file in the cache is removed."""
        cache = tmp_path / "cache"
        cache.mkdir()
        f1 = cache / "alpha.bin"
        f2 = cache / "beta.bin"
        f1.write_bytes(b"a")
        f2.write_bytes(b"b")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        removed = recount3_cache_rm(config=cfg, predicate=None)
        assert sorted(p.name for p in removed) == ["alpha.bin", "beta.bin"]
        assert not f1.exists()
        assert not f2.exists()

    def test_predicate_selects_subset(self, tmp_path: Path) -> None:
        """Only files for which the predicate returns True are removed."""
        cache = tmp_path / "cache"
        cache.mkdir()
        keep = cache / "keep.txt"
        remove = cache / "remove.gz"
        keep.write_bytes(b"keep")
        remove.write_bytes(b"remove")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        removed = recount3_cache_rm(
            config=cfg,
            predicate=lambda p: p.suffix == ".gz",
        )
        assert len(removed) == 1
        assert removed[0].name == "remove.gz"
        assert not remove.exists()
        assert keep.exists()

    def test_predicate_returns_false_for_all_keeps_everything(
        self, tmp_path: Path
    ) -> None:
        """When the predicate always returns False, no files are removed."""
        cache = tmp_path / "cache"
        cache.mkdir()
        f = cache / "important.dat"
        f.write_bytes(b"precious")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        removed = recount3_cache_rm(config=cfg, predicate=lambda p: False)
        assert removed == []
        assert f.exists()

    def test_dry_run_does_not_delete_files(self, tmp_path: Path) -> None:
        """dry_run=True reports candidates but leaves all files on disk."""
        cache = tmp_path / "cache"
        cache.mkdir()
        f1 = cache / "a.gz"
        f2 = cache / "b.gz"
        f1.write_bytes(b"a")
        f2.write_bytes(b"b")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        candidates = recount3_cache_rm(config=cfg, dry_run=True)
        assert sorted(p.name for p in candidates) == ["a.gz", "b.gz"]
        assert f1.exists()
        assert f2.exists()

    def test_dry_run_with_predicate(self, tmp_path: Path) -> None:
        """dry_run + predicate reports the subset that would be removed."""
        cache = tmp_path / "cache"
        cache.mkdir()
        keep = cache / "keep.txt"
        would_remove = cache / "gone.gz"
        keep.write_bytes(b"keep")
        would_remove.write_bytes(b"gone")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        candidates = recount3_cache_rm(
            config=cfg,
            predicate=lambda p: p.suffix == ".gz",
            dry_run=True,
        )
        assert len(candidates) == 1
        assert candidates[0].name == "gone.gz"
        assert would_remove.exists()

    def test_removes_files_in_subdirectories(self, tmp_path: Path) -> None:
        """Files nested inside subdirectories are also candidates for removal."""
        cache = tmp_path / "cache"
        sub = cache / "org" / "proj"
        sub.mkdir(parents=True)
        deep_file = sub / "counts.gz"
        deep_file.write_bytes(b"data")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        removed = recount3_cache_rm(config=cfg)
        assert len(removed) == 1
        assert not deep_file.exists()

    def test_directories_are_not_removed(self, tmp_path: Path) -> None:
        """Only files are removed; directories are left in place."""
        cache = tmp_path / "cache"
        subdir = cache / "subdir"
        subdir.mkdir(parents=True)
        f = subdir / "file.gz"
        f.write_bytes(b"data")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        recount3_cache_rm(config=cfg)
        assert subdir.is_dir()
        assert not f.exists()

    def test_returned_list_is_sorted(self, tmp_path: Path) -> None:
        """The returned list is sorted regardless of filesystem order."""
        cache = tmp_path / "cache"
        cache.mkdir()
        for name in ("zzz.bin", "aaa.bin", "mmm.bin"):
            (cache / name).write_bytes(b"x")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        removed = recount3_cache_rm(config=cfg)
        assert removed == sorted(removed, key=str)

    def test_with_no_config_uses_default_config(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """config=None must fall back to default_config()."""
        cache = tmp_path / "dc_rm_cache"
        cache.mkdir()
        f = cache / "tmp.gz"
        f.write_bytes(b"data")
        monkeypatch.setenv("RECOUNT3_CACHE_DIR", str(cache))
        removed = recount3_cache_rm(config=None)
        assert len(removed) == 1
        assert not f.exists()

    def test_returns_list_of_path_objects(self, tmp_path: Path) -> None:
        """Returned items are :class:`pathlib.Path` instances."""
        cache = tmp_path / "cache"
        cache.mkdir()
        (cache / "file.bin").write_bytes(b"data")
        cfg = Config(
            base_url=_BASE_URL,
            timeout=30,
            insecure_ssl=False,
            max_retries=3,
            user_agent="ua",
            cache_dir=cache,
            cache_disabled=False,
        )
        removed = recount3_cache_rm(config=cfg)
        assert all(isinstance(p, Path) for p in removed)
