# Copyright (c) 2026, Alexander Alsalihi, Robert M. Flight, Hunter N.B. Moseley.
# All rights reserved.
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
#   may be used to endorse or promote products derived from this software
#   without specific prior written permission.
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

import argparse
import dataclasses
import io
import json
import logging
from pathlib import Path
from typing import Any
from unittest import mock

import numpy as np
import pandas as pd
import pytest
import scipy.sparse

from recount3._descriptions import R3ResourceDescription
from recount3.config import Config, default_config
from recount3.errors import (
    CompatibilityError,
    ConfigurationError,
    LoadError,
    Recount3Error,
)
from recount3.resource import R3Resource

from recount3.cli import (
    _build_config_from_env_and_flags,
    _build_parser,
    _cmd_bundle_rse,
    _cmd_bundle_se,
    _cmd_bundle_stack_counts,
    _cmd_download,
    _cmd_ids,
    _cmd_search,
    _cmd_smoke_test,
    _dispatch,
    _init_logging,
    _iter_manifest,
    _parse_filters,
    _reject_unknown_filters,
    _SEARCH_MODE_FILTERS,
    _resource_from_dict,
    _write_jsonl,
    _write_tsv,
    main,
)

_BASE_URL = "http://duffel.rail.bio/recount3/"


def _make_cfg(tmp_path: Path) -> Config:
    return Config(
        base_url=_BASE_URL,
        timeout=30,
        insecure_ssl=False,
        max_retries=3,
        user_agent="test-agent/1.0",
        cache_dir=tmp_path / "cache",
        cache_disabled=True,
        chunk_size=1024,
    )


def _make_annotation_resource(cfg: Config | None = None) -> R3Resource:
    desc = R3ResourceDescription(
        resource_type="annotations",
        organism="human",
        genomic_unit="gene",
        annotation_extension="G026",
    )
    return R3Resource(description=desc, config=cfg)


def _make_metadata_resource(cfg: Config | None = None) -> R3Resource:
    desc = R3ResourceDescription(
        resource_type="metadata_files",
        organism="human",
        data_source="sra",
        project="SRP000001",
        table_name="recount_qc",
    )
    return R3Resource(description=desc, config=cfg)


def _make_namespace(**kwargs: Any) -> argparse.Namespace:
    """Build a minimal Namespace with defaults for global flags."""
    defaults = dict(
        base_url=None,
        cache_dir=None,
        timeout=None,
        retries=None,
        insecure_ssl=False,
        user_agent=None,
        chunk_size=None,
        quiet=False,
        verbose=False,
    )
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


class TestBuildParser:
    def test_returns_parser(self) -> None:
        parser = _build_parser()
        assert isinstance(parser, argparse.ArgumentParser)

    def test_version_flag(self, capsys: pytest.CaptureFixture[str]) -> None:
        parser = _build_parser()
        with pytest.raises(SystemExit) as exc:
            parser.parse_args(["--version"])
        assert exc.value.code == 0

    def test_ids_subcommand_defaults(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["ids"])
        assert args.command == "ids"
        assert args.organism == ""
        assert args.samples_out is None
        assert args.projects_out is None

    def test_search_subcommand(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["search", "annotations", "organism=human"])
        assert args.command == "search"
        assert args.mode == "annotations"
        assert args.filters == ["organism=human"]
        assert args.format == "jsonl"
        assert args.output is None
        assert args.outdir is None

    def test_download_subcommand_manifest(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(
            ["download", "--from=manifest.jsonl", "--dest=/tmp"]
        )
        assert args.command == "download"
        assert args.manifest == "manifest.jsonl"
        assert args.dest == "/tmp"
        assert args.jobs == 8
        assert args.cache == "enable"

    def test_download_subcommand_inline(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["download", "--inline={}", "--dest=/tmp"])
        assert args.inline == "{}"

    def test_bundle_stack_counts(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(
            ["bundle", "stack-counts", "--from=m.jsonl", "--out=out.parquet"]
        )
        assert args.command == "bundle"
        assert args.bundle_cmd == "stack-counts"
        assert args.compat == "family"
        assert args.join == "inner"
        assert args.axis == 1
        assert not args.verify_integrity
        assert not args.densify

    def test_bundle_stack_counts_densify(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(
            [
                "bundle",
                "stack-counts",
                "--from=m.jsonl",
                "--densify",
                "--out=out.parquet",
            ]
        )
        assert args.densify

    def test_bundle_se(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(
            [
                "bundle",
                "se",
                "--from=m.jsonl",
                "--genomic-unit=gene",
                "--out=out.pkl",
            ]
        )
        assert args.bundle_cmd == "se"
        assert args.genomic_unit == "gene"
        assert args.assay_name == "raw_counts"

    def test_bundle_rse(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(
            [
                "bundle",
                "rse",
                "--from=m.jsonl",
                "--genomic-unit=gene",
                "--out=out.pkl",
            ]
        )
        assert args.bundle_cmd == "rse"
        assert not args.allow_fallback_to_se

    def test_smoke_test_defaults(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["smoke-test"])
        assert args.command == "smoke-test"
        assert args.limit == 1
        assert args.dest == "./recount3-smoke"

    def test_smoke_test_custom_dest(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["smoke-test", "--dest=/tmp/my-smoke"])
        assert args.dest == "/tmp/my-smoke"

    def test_global_flags_parsed(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(
            [
                "--base-url=http://example.org/",
                "--cache-dir=/tmp/cache",
                "--timeout=99",
                "--retries=5",
                "--insecure-ssl",
                "--user-agent=myagent",
                "--chunk-size=8192",
                "--quiet",
                "--verbose",
                "smoke-test",
            ]
        )
        assert args.base_url == "http://example.org/"
        assert args.cache_dir == "/tmp/cache"
        assert args.timeout == 99
        assert args.retries == 5
        assert args.insecure_ssl is True
        assert args.user_agent == "myagent"
        assert args.chunk_size == 8192
        assert args.quiet is True
        assert args.verbose is True

    def test_search_output_and_outdir_mutually_exclusive(self) -> None:
        parser = _build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(
                ["search", "annotations", "--output=a.jsonl", "--outdir=/tmp"]
            )


class TestBuildConfigFromEnvAndFlags:
    @pytest.mark.parametrize("backend", ["filesystem", "pybiocfilecache"])
    @pytest.mark.parametrize("directory_source", ["default", "env", "flag"])
    def test_cache_directory_follows_selected_backend(
        self, monkeypatch, tmp_path, backend, directory_source
    ):
        """CLI backend flags select defaults without losing path overrides."""
        other = (
            "filesystem" if backend == "pybiocfilecache" else "pybiocfilecache"
        )
        monkeypatch.setenv("RECOUNT3_CACHE_BACKEND", other)
        monkeypatch.delenv("RECOUNT3_CACHE_DIR", raising=False)
        monkeypatch.setenv("R_USER_CACHE_DIR", str(tmp_path / "shared"))
        flags = ["--cache-backend", backend]
        expected = (
            tmp_path / "shared/R/recount3"
            if backend == "pybiocfilecache"
            else Path.home() / ".cache/recount3/files"
        )
        if directory_source != "default":
            monkeypatch.setenv("RECOUNT3_CACHE_DIR", str(tmp_path / "env"))
            expected = tmp_path / "env"
        if directory_source == "flag":
            flags.extend(["--cache-dir", str(tmp_path / "flag")])
            expected = tmp_path / "flag"
        args = _build_parser().parse_args(
            [*flags, "download", "--inline", "{}"]
        )
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.cache_backend == backend
        assert cfg.cache_dir == expected.resolve()

    def test_all_none_uses_defaults(self) -> None:
        args = _make_namespace()
        cfg = _build_config_from_env_and_flags(args)
        assert isinstance(cfg, Config)
        base = default_config()
        assert cfg.base_url == base.base_url

    def test_base_url_overrides_env(self) -> None:
        args = _make_namespace(base_url="http://custom.example.org/r3/")
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.base_url == "http://custom.example.org/r3/"

    def test_cache_dir_overrides_env(self, tmp_path: Path) -> None:
        args = _make_namespace(cache_dir=str(tmp_path))
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.cache_dir == tmp_path.expanduser().resolve()

    def test_timeout_overrides_default(self) -> None:
        args = _make_namespace(timeout=42)
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.timeout == 42

    def test_retries_overrides_default(self) -> None:
        args = _make_namespace(retries=7)
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.max_retries == 7

    def test_insecure_ssl_true(self) -> None:
        args = _make_namespace(insecure_ssl=True)
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.insecure_ssl is True

    def test_insecure_ssl_false_inherits_env(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("RECOUNT3_INSECURE_SSL", "1")
        args = _make_namespace(insecure_ssl=False)
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.insecure_ssl is True

    def test_user_agent_overrides_default(self) -> None:
        args = _make_namespace(user_agent="mybot/2.0")
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.user_agent == "mybot/2.0"

    def test_chunk_size_overrides_default(self) -> None:
        args = _make_namespace(chunk_size=512)
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.chunk_size == 512

    def test_cache_disable_env_var_one(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("RECOUNT3_CACHE_DISABLE", "1")
        args = _make_namespace()
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.cache_disabled is True

    def test_cache_disable_env_var_not_one(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("RECOUNT3_CACHE_DISABLE", "0")
        args = _make_namespace()
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.cache_disabled is False

    def test_cache_disable_env_var_other(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("RECOUNT3_CACHE_DISABLE", "yes")
        args = _make_namespace()
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.cache_disabled is False

    def test_cache_backend_env_var_used_when_flag_is_absent(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("RECOUNT3_CACHE_BACKEND", "pybiocfilecache")
        args = _build_parser().parse_args(["download", "--inline", "{}"])
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.cache_backend == "pybiocfilecache"

    def test_cache_backend_flag_overrides_env_var(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("RECOUNT3_CACHE_BACKEND", "pybiocfilecache")
        args = _build_parser().parse_args(
            ["--cache-backend", "filesystem", "download", "--inline", "{}"]
        )
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.cache_backend == "filesystem"

    def test_invalid_cache_dir_raises_configuration_error(self) -> None:
        args = _make_namespace(cache_dir="/some/path")
        with mock.patch.object(Path, "resolve", side_effect=OSError("broken")):
            with pytest.raises(ConfigurationError, match="Invalid cache"):
                _build_config_from_env_and_flags(args)


class TestInitLogging:
    def test_quiet_sets_warning(self) -> None:
        args = _make_namespace(quiet=True, verbose=False)
        with mock.patch("logging.basicConfig") as m:
            _init_logging(args)
        m.assert_called_once()
        _, kwargs = m.call_args
        assert kwargs["level"] == logging.WARNING

    def test_verbose_sets_debug(self) -> None:
        args = _make_namespace(quiet=False, verbose=True)
        with mock.patch("logging.basicConfig") as m:
            _init_logging(args)
        _, kwargs = m.call_args
        assert kwargs["level"] == logging.DEBUG

    def test_default_sets_info(self) -> None:
        args = _make_namespace(quiet=False, verbose=False)
        with mock.patch("logging.basicConfig") as m:
            _init_logging(args)
        _, kwargs = m.call_args
        assert kwargs["level"] == logging.INFO

    def test_quiet_wins_over_verbose(self) -> None:
        args = _make_namespace(quiet=True, verbose=True)
        with mock.patch("logging.basicConfig") as m:
            _init_logging(args)
        _, kwargs = m.call_args
        assert kwargs["level"] == logging.DEBUG


class TestParseFilters:
    def test_empty_list_returns_empty_dict(self) -> None:
        assert _parse_filters([]) == {}

    def test_valid_single_token(self) -> None:
        assert _parse_filters(["organism=human"]) == {"organism": "human"}

    def test_key_lowercased(self) -> None:
        result = _parse_filters(["Organism=Human"])
        assert "organism" in result
        assert result["organism"] == "Human"

    def test_multiple_tokens(self) -> None:
        result = _parse_filters(["organism=human", "data_source=sra"])
        assert result == {"organism": "human", "data_source": "sra"}

    def test_value_with_equals(self) -> None:
        # Split on first "=" only.
        result = _parse_filters(["key=a=b"])
        assert result == {"key": "a=b"}

    def test_missing_equals_raises(self) -> None:
        with pytest.raises(ValueError, match="Expected.*key=value"):
            _parse_filters(["badtoken"])

    def test_error_message_mentions_spaces(self) -> None:
        with pytest.raises(ValueError, match="no spaces") as exc_info:
            _parse_filters(["project"])
        assert "ensure no spaces" in str(exc_info.value)

    def test_empty_key_raises(self) -> None:
        with pytest.raises(ValueError, match="Empty filter key"):
            _parse_filters(["=value"])

    def test_whitespace_stripped_from_key(self) -> None:
        result = _parse_filters(["  key  =val"])
        assert "key" in result


class TestRejectUnknownFilters:
    """A selector a mode never reads must not be dropped in silence."""

    def test_accepts_the_keys_a_mode_reads(self) -> None:
        _reject_unknown_filters(
            "project",
            {
                "organism": "human",
                "data_source": "sra",
                "project": "SRP009615",
                "genomic_unit": "gene",
                "annotation_extension": "G026",
            },
        )

    def test_accepts_no_filters(self) -> None:
        _reject_unknown_filters("sources", {})

    def test_rejects_an_unknown_key(self) -> None:
        with pytest.raises(ValueError, match="Unknown filter") as exc_info:
            _reject_unknown_filters("sources", {"organism": "human", "x": "1"})
        assert "'x'" in str(exc_info.value)

    def test_names_the_cli_spelling_for_a_python_keyword(self) -> None:
        """discover() takes genomic_units; the CLI selector is singular."""
        with pytest.raises(ValueError) as exc_info:
            _reject_unknown_filters(
                "project",
                {
                    "organism": "human",
                    "data_source": "sra",
                    "project": "SRP009615",
                    "genomic_units": "gene",
                    "annotations": "G026",
                },
            )
        message = str(exc_info.value)
        assert "'genomic_units' (the Python API spelling; use" in message
        assert "'genomic_unit'" in message
        assert "'annotations'" in message and "'annotation'" in message

    def test_lists_the_keys_the_mode_accepts(self) -> None:
        with pytest.raises(ValueError) as exc_info:
            _reject_unknown_filters("bigwig", {"annotation": "gencode_v26"})
        message = str(exc_info.value)
        assert "This mode reads: data_source, organism, project, sample" in (
            message
        )

    def test_key_valid_for_another_mode_is_still_rejected(self) -> None:
        """'sample' belongs to bigwig mode, not gene-exon."""
        with pytest.raises(ValueError, match="'sample'"):
            _reject_unknown_filters("gene-exon", {"sample": "SRR387777"})

    def test_every_search_mode_declares_its_selectors(self) -> None:
        """The table and the parser's MODE choices must not drift apart.

        A mode missing from the table would silently go unvalidated again;
        a stale entry would describe a subcommand that no longer exists.
        """
        parser = _build_parser()
        subparsers = [
            action
            for action in parser._actions  # pylint: disable=protected-access
            if isinstance(action, argparse._SubParsersAction)
        ][0]
        search = subparsers.choices["search"]
        mode_action = [
            action
            for action in search._actions  # pylint: disable=protected-access
            if action.dest == "mode"
        ][0]

        assert set(mode_action.choices) == set(_SEARCH_MODE_FILTERS)

    def test_unregistered_mode_is_not_second_guessed(self) -> None:
        """An unknown mode is _cmd_search's error to raise, not this one's."""
        _reject_unknown_filters("not-a-mode", {"anything": "goes"})


class TestResourceFromDict:
    def test_creates_resource(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        d = {
            "resource_type": "annotations",
            "organism": "human",
            "genomic_unit": "gene",
            "annotation_extension": "G026",
        }
        res = _resource_from_dict(d, cfg)
        assert isinstance(res, R3Resource)
        assert res.config is cfg

    def test_strips_url_and_arcname(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        d = {
            "resource_type": "annotations",
            "organism": "human",
            "genomic_unit": "gene",
            "annotation_extension": "G026",
            "url": "http://should-be-stripped.example/",
            "arcname": "also-stripped",
        }
        res = _resource_from_dict(d, cfg)
        assert isinstance(res, R3Resource)
        assert res.url != "http://should-be-stripped.example/"

    def test_invalid_resource_type_raises(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        d = {"resource_type": "nonexistent_type"}
        with pytest.raises((ValueError, KeyError)):
            _resource_from_dict(d, cfg)


class TestWriteJsonl:
    def test_writes_to_stdout_when_out_none(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        _write_jsonl([res], None)
        captured = capsys.readouterr()
        obj = json.loads(captured.out.strip())
        assert obj["resource_type"] == "annotations"
        assert "url" in obj
        assert "arcname" in obj

    def test_writes_to_file_when_out_set(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        out = tmp_path / "subdir" / "out.jsonl"
        _write_jsonl([res], out)
        assert out.exists()
        obj = json.loads(out.read_text().strip())
        assert obj["resource_type"] == "annotations"

    def test_dataclass_description_asdict_used(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        assert dataclasses.is_dataclass(res.description)
        _write_jsonl([res], None)
        captured = capsys.readouterr()
        obj = json.loads(captured.out.strip())
        assert obj["organism"] == "human"

    def test_non_dataclass_description_yields_empty_body(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        non_dc_res = mock.MagicMock()
        non_dc_res.description = object()
        non_dc_res.url = "http://example.org/file.gz"
        non_dc_res.arcname = "file.gz"
        with mock.patch(
            "recount3.cli.dataclasses.is_dataclass", return_value=False
        ):
            _write_jsonl([non_dc_res], None)
        captured = capsys.readouterr()
        obj = json.loads(captured.out.strip())
        assert obj["url"] == "http://example.org/file.gz"
        assert "organism" not in obj

    def test_multiple_resources(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        res1 = _make_annotation_resource(cfg)
        res2 = _make_metadata_resource(cfg)
        out = tmp_path / "multi.jsonl"
        _write_jsonl([res1, res2], out)
        lines = [l for l in out.read_text().splitlines() if l.strip()]
        assert len(lines) == 2


class TestWriteTsv:
    def test_writes_to_stdout_when_out_none(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        _write_tsv([res], None)
        captured = capsys.readouterr()
        lines = captured.out.splitlines()
        assert len(lines) == 2
        assert "resource_type" in lines[0]

    def test_writes_to_file(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        out = tmp_path / "sub" / "out.tsv"
        _write_tsv([res], out)
        assert out.exists()
        text = out.read_text()
        assert "resource_type" in text
        assert "human" in text

    def test_non_dataclass_description_empty_row(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        non_dc_res = mock.MagicMock()
        non_dc_res.description = object()
        non_dc_res.url = "http://example.org/file.gz"
        non_dc_res.arcname = "file.gz"
        with mock.patch(
            "recount3.cli.dataclasses.is_dataclass", return_value=False
        ):
            _write_tsv([non_dc_res], None)
        captured = capsys.readouterr()
        lines = captured.out.splitlines()
        assert len(lines) == 2
        data_row = lines[1]
        assert "http://example.org/file.gz" in data_row

    def test_missing_fields_become_empty_string(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        _write_tsv([res], None)
        captured = capsys.readouterr()
        data_line = captured.out.splitlines()[1]
        assert "\t\t" in data_line or data_line.count("\t") >= 2


class TestIterManifest:
    def test_reads_from_file(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        line = json.dumps(
            {
                "resource_type": "annotations",
                "organism": "human",
                "genomic_unit": "gene",
                "annotation_extension": "G026",
            }
        )
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text(line + "\n")
        resources = list(_iter_manifest(str(manifest), cfg))
        assert len(resources) == 1
        assert isinstance(resources[0], R3Resource)

    def test_skips_blank_lines_in_file(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        line = json.dumps(
            {
                "resource_type": "annotations",
                "organism": "human",
                "genomic_unit": "gene",
                "annotation_extension": "G026",
            }
        )
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text("\n" + line + "\n\n")
        resources = list(_iter_manifest(str(manifest), cfg))
        assert len(resources) == 1

    def test_reads_from_stdin(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        line = json.dumps(
            {
                "resource_type": "annotations",
                "organism": "human",
                "genomic_unit": "gene",
                "annotation_extension": "G026",
            }
        )
        fake_stdin = io.StringIO(line + "\n")
        with mock.patch("sys.stdin", fake_stdin):
            resources = list(_iter_manifest("-", cfg))
        assert len(resources) == 1

    def test_skips_blank_lines_in_stdin(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        line = json.dumps(
            {
                "resource_type": "annotations",
                "organism": "human",
                "genomic_unit": "gene",
                "annotation_extension": "G026",
            }
        )
        fake_stdin = io.StringIO("\n" + line + "\n\n")
        with mock.patch("sys.stdin", fake_stdin):
            resources = list(_iter_manifest("-", cfg))
        assert len(resources) == 1

    def test_strips_url_arcname_from_manifest(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        line = json.dumps(
            {
                "resource_type": "annotations",
                "organism": "human",
                "genomic_unit": "gene",
                "annotation_extension": "G026",
                "url": "http://old-url/",
                "arcname": "old-arcname",
            }
        )
        manifest = tmp_path / "m.jsonl"
        manifest.write_text(line + "\n")
        resources = list(_iter_manifest(str(manifest), cfg))
        assert len(resources) == 1
        assert resources[0].url != "http://old-url/"


class TestCmdIds:
    def _make_args(self, **kw: Any) -> argparse.Namespace:
        defaults = dict(organism="", samples_out=None, projects_out=None)
        defaults.update(kw)
        return argparse.Namespace(**defaults)

    def test_writes_to_stdout_when_no_out(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        args = self._make_args()
        with mock.patch(
            "recount3.search.create_sample_project_lists",
            return_value=(["SRR001", "SRR002"], ["SRP001"]),
        ):
            code = _cmd_ids(args, cfg)
        assert code == 0
        captured = capsys.readouterr()
        assert "SRR001" in captured.out
        assert "SRP001" in captured.out

    def test_writes_samples_to_file(self, tmp_path: Path) -> None:
        samples_file = tmp_path / "samples.txt"
        args = self._make_args(samples_out=str(samples_file))
        cfg = _make_cfg(tmp_path)
        with mock.patch(
            "recount3.search.create_sample_project_lists",
            return_value=(["SRR001"], ["SRP001"]),
        ):
            code = _cmd_ids(args, cfg)
        assert code == 0
        assert samples_file.read_text() == "SRR001"

    def test_writes_projects_to_file(self, tmp_path: Path) -> None:
        projects_file = tmp_path / "projects.txt"
        args = self._make_args(projects_out=str(projects_file))
        cfg = _make_cfg(tmp_path)
        with mock.patch(
            "recount3.search.create_sample_project_lists",
            return_value=(["SRR001"], ["SRP001"]),
        ):
            code = _cmd_ids(args, cfg)
        assert code == 0
        assert projects_file.read_text() == "SRP001"

    def test_organism_filter_passed_through(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        args = self._make_args(organism="human")
        with mock.patch(
            "recount3.search.create_sample_project_lists",
            return_value=([], []),
        ) as m:
            _cmd_ids(args, cfg)
        m.assert_called_once_with(organism="human")

    def test_writes_both_to_files(self, tmp_path: Path) -> None:
        samples_file = tmp_path / "s.txt"
        projects_file = tmp_path / "p.txt"
        args = self._make_args(
            samples_out=str(samples_file), projects_out=str(projects_file)
        )
        cfg = _make_cfg(tmp_path)
        with mock.patch(
            "recount3.search.create_sample_project_lists",
            return_value=(["SRR001", "SRR002"], ["SRP001", "SRP002"]),
        ):
            code = _cmd_ids(args, cfg)
        assert code == 0
        assert "SRR001" in samples_file.read_text()
        assert "SRP001" in projects_file.read_text()


def _make_search_args(**kw: Any) -> argparse.Namespace:
    defaults = dict(
        mode="annotations",
        filters=[],
        format="jsonl",
        output=None,
        outdir=None,
    )
    defaults.update(kw)
    return argparse.Namespace(**defaults)


class TestCmdSearchRejectsUnreadFilters:
    """_cmd_search must refuse selectors before it emits a manifest."""

    def test_python_keyword_names_do_not_reach_search_project_all(
        self, tmp_path: Path
    ) -> None:
        """Previously these were dropped and the full default set emitted."""
        cfg = _make_cfg(tmp_path)
        args = _make_search_args(
            mode="project",
            filters=[
                "organism=human",
                "data_source=sra",
                "project=SRP009615",
                "genomic_units=gene",
                "annotations=G026",
            ],
        )
        with mock.patch(
            "recount3.cli.r3_search.search_project_all"
        ) as mock_search:
            with pytest.raises(ValueError, match="Unknown filter"):
                _cmd_search(args, cfg)
        mock_search.assert_not_called()


class TestCmdSearchAnnotations:
    def test_annotations_mode_success(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource()
        args = _make_search_args(
            mode="annotations",
            filters=[
                "organism=human",
                "genomic_unit=gene",
                "annotation_extension=G026",
            ],
        )
        with mock.patch(
            "recount3.cli.r3_search.search_annotations", return_value=[res]
        ):
            code = _cmd_search(args, cfg)
        assert code == 0
        captured = capsys.readouterr()
        obj = json.loads(captured.out.strip())
        assert obj["resource_type"] == "annotations"

    def test_annotations_missing_filter_raises(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        args = _make_search_args(
            mode="annotations",
            filters=["organism=human"],
        )
        with mock.patch(
            "recount3.cli.r3_search.search_annotations", return_value=[]
        ):
            with pytest.raises(ValueError, match="Missing required filters"):
                _cmd_search(args, cfg)


class TestCmdSearchGeneExon:
    def test_gene_exon_mode_with_default_annotation_extension(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        desc = R3ResourceDescription(
            resource_type="count_files_gene_or_exon",
            organism="human",
            data_source="sra",
            genomic_unit="gene",
            project="SRP000001",
            annotation_extension="G026",
        )
        res = R3Resource(description=desc)
        args = _make_search_args(
            mode="gene-exon",
            filters=[
                "organism=human",
                "data_source=sra",
                "genomic_unit=gene",
                "project=SRP000001",
            ],
        )
        with mock.patch(
            "recount3.cli.r3_search.search_count_files_gene_or_exon",
            return_value=[res],
        ) as m:
            code = _cmd_search(args, cfg)
        assert code == 0
        _, kwargs = m.call_args
        assert kwargs["annotation_extension"] == ("G026",)

    def test_gene_exon_mode_with_explicit_annotation_extension(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        desc = R3ResourceDescription(
            resource_type="count_files_gene_or_exon",
            organism="human",
            data_source="sra",
            genomic_unit="gene",
            project="SRP000001",
            annotation_extension="G029",
        )
        res = R3Resource(description=desc)
        args = _make_search_args(
            mode="gene-exon",
            filters=[
                "organism=human",
                "data_source=sra",
                "genomic_unit=gene",
                "project=SRP000001",
                "annotation_extension=G029",
            ],
        )
        with mock.patch(
            "recount3.cli.r3_search.search_count_files_gene_or_exon",
            return_value=[res],
        ) as m:
            code = _cmd_search(args, cfg)
        assert code == 0
        _, kwargs = m.call_args
        assert kwargs["annotation_extension"] == "G029"

    def test_gene_exon_missing_filter_raises(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        args = _make_search_args(
            mode="gene-exon",
            filters=["organism=human"],
        )
        with pytest.raises(ValueError, match="Missing required filters"):
            _cmd_search(args, cfg)


class TestCmdSearchJunctions:
    def test_junctions_mode_defaults(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        desc = R3ResourceDescription(
            resource_type="count_files_junctions",
            organism="human",
            data_source="sra",
            project="SRP000001",
            junction_type="ALL",
            junction_extension="MM",
        )
        res = R3Resource(description=desc)
        args = _make_search_args(
            mode="junctions",
            filters=["organism=human", "data_source=sra", "project=SRP000001"],
        )
        with mock.patch(
            "recount3.cli.r3_search.search_count_files_junctions",
            return_value=[res],
        ) as m:
            code = _cmd_search(args, cfg)
        assert code == 0
        _, kwargs = m.call_args
        assert kwargs["junction_type"] == "ALL"
        assert kwargs["junction_extension"] == "MM"

    def test_junctions_mode_with_explicit_type_and_extension(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        desc = R3ResourceDescription(
            resource_type="count_files_junctions",
            organism="human",
            data_source="sra",
            project="SRP000001",
            junction_type="ALL",
            junction_extension="RR",
        )
        res = R3Resource(description=desc)
        args = _make_search_args(
            mode="junctions",
            filters=[
                "organism=human",
                "data_source=sra",
                "project=SRP000001",
                "junction_type=ALL",
                "junction_extension=RR",
            ],
        )
        with mock.patch(
            "recount3.cli.r3_search.search_count_files_junctions",
            return_value=[res],
        ) as m:
            code = _cmd_search(args, cfg)
        assert code == 0
        _, kwargs = m.call_args
        assert kwargs["junction_extension"] == "RR"

    def test_junctions_missing_filter_raises(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        args = _make_search_args(mode="junctions", filters=["organism=human"])
        with pytest.raises(ValueError, match="Missing required filters"):
            _cmd_search(args, cfg)


class TestCmdSearchMetadata:
    def test_metadata_mode_success(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_metadata_resource()
        args = _make_search_args(
            mode="metadata",
            filters=[
                "organism=human",
                "data_source=sra",
                "table_name=recount_qc",
                "project=SRP000001",
            ],
        )
        with mock.patch(
            "recount3.cli.r3_search.search_metadata_files", return_value=[res]
        ):
            code = _cmd_search(args, cfg)
        assert code == 0

    def test_metadata_missing_filter_raises(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        args = _make_search_args(mode="metadata", filters=["organism=human"])
        with pytest.raises(ValueError, match="Missing required filters"):
            _cmd_search(args, cfg)


class TestCmdSearchBigwig:
    def test_bigwig_mode_success(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        desc = R3ResourceDescription(
            resource_type="bigwig_files",
            organism="human",
            data_source="sra",
            project="SRP000001",
            sample="SRR000001",
        )
        res = R3Resource(description=desc)
        args = _make_search_args(
            mode="bigwig",
            filters=[
                "organism=human",
                "data_source=sra",
                "project=SRP000001",
                "sample=SRR000001",
            ],
        )
        with mock.patch(
            "recount3.cli.r3_search.search_bigwig_files", return_value=[res]
        ):
            code = _cmd_search(args, cfg)
        assert code == 0

    def test_bigwig_missing_filter_raises(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        args = _make_search_args(mode="bigwig", filters=["organism=human"])
        with pytest.raises(ValueError, match="Missing required filters"):
            _cmd_search(args, cfg)


class TestCmdSearchSources:
    def test_sources_mode_success(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        desc = R3ResourceDescription(
            resource_type="data_sources", organism="human"
        )
        res = R3Resource(description=desc)
        args = _make_search_args(mode="sources", filters=["organism=human"])
        with mock.patch(
            "recount3.cli.r3_search.search_data_sources", return_value=[res]
        ):
            code = _cmd_search(args, cfg)
        assert code == 0

    def test_sources_missing_organism_raises(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        args = _make_search_args(mode="sources", filters=[])
        with pytest.raises(ValueError, match="Missing required filters"):
            _cmd_search(args, cfg)


class TestCmdSearchSourceMeta:
    def test_source_meta_success(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        desc = R3ResourceDescription(
            resource_type="data_source_metadata",
            organism="human",
            data_source="sra",
        )
        res = R3Resource(description=desc)
        args = _make_search_args(
            mode="source-meta",
            filters=["organism=human", "data_source=sra"],
        )
        with mock.patch(
            "recount3.cli.r3_search.search_data_source_metadata",
            return_value=[res],
        ):
            code = _cmd_search(args, cfg)
        assert code == 0

    def test_source_meta_missing_filter_raises(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        args = _make_search_args(mode="source-meta", filters=["organism=human"])
        with pytest.raises(ValueError, match="Missing required filters"):
            _cmd_search(args, cfg)


class TestCmdSearchProject:
    def _project_args(
        self, extra_filters: list[str] | None = None
    ) -> argparse.Namespace:
        base_filters = [
            "organism=human",
            "data_source=sra",
            "project=SRP000001",
        ]
        return _make_search_args(
            mode="project",
            filters=base_filters + (extra_filters or []),
        )

    def _run(
        self,
        tmp_path: Path,
        args: argparse.Namespace,
        found: list[R3Resource] | None = None,
    ) -> tuple[int, Any]:
        cfg = _make_cfg(tmp_path)
        if found is None:
            found = []
        with mock.patch(
            "recount3.cli.r3_search.search_project_all", return_value=found
        ) as m:
            code = _cmd_search(args, cfg)
        return code, m

    def test_project_mode_minimal_success(self, tmp_path: Path) -> None:
        args = self._project_args()
        code, m = self._run(tmp_path, args)
        assert code == 0

    def test_project_missing_filter_raises(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        args = _make_search_args(mode="project", filters=["organism=human"])
        with mock.patch(
            "recount3.cli.r3_search.search_project_all", return_value=[]
        ):
            with pytest.raises(ValueError, match="Missing required filters"):
                _cmd_search(args, cfg)

    def test_as_bool_true_values(self, tmp_path: Path) -> None:
        for val in ("1", "true", "t", "yes", "y", "on"):
            args = self._project_args([f"include_bigwig={val}"])
            code, m = self._run(tmp_path, args)
            assert code == 0
            _, kwargs = m.call_args
            assert (
                kwargs["include_bigwig"] is True
            ), f"Expected True for {val!r}"

    def test_as_bool_false_values(self, tmp_path: Path) -> None:
        for val in ("0", "false", "no", "whatever"):
            args = self._project_args([f"include_bigwig={val}"])
            code, m = self._run(tmp_path, args)
            assert code == 0
            _, kwargs = m.call_args
            assert (
                kwargs["include_bigwig"] is False
            ), f"Expected False for {val!r}"

    def test_as_bool_none_uses_default(self, tmp_path: Path) -> None:
        args = self._project_args()
        code, m = self._run(tmp_path, args)
        assert code == 0
        _, kwargs = m.call_args
        assert kwargs["include_bigwig"] is False

    def test_as_bool_include_metadata_default_true(
        self, tmp_path: Path
    ) -> None:
        args = self._project_args()
        code, m = self._run(tmp_path, args)
        _, kwargs = m.call_args
        assert kwargs["include_metadata"] is True

    def test_csv_or_default_none_returns_default(self, tmp_path: Path) -> None:
        args = self._project_args()
        code, m = self._run(tmp_path, args)
        _, kwargs = m.call_args
        assert kwargs["genomic_units"] == ("gene", "exon")

    def test_csv_or_default_empty_string_returns_default(
        self, tmp_path: Path
    ) -> None:
        args = self._project_args(["genomic_unit="])
        code, m = self._run(tmp_path, args)
        _, kwargs = m.call_args
        assert kwargs["genomic_units"] == ("gene", "exon")

    def test_csv_or_default_comma_separated(self, tmp_path: Path) -> None:
        args = self._project_args(["junction_extension=MM,RR"])
        code, m = self._run(tmp_path, args)
        _, kwargs = m.call_args
        assert kwargs["junction_extension"] == ("MM", "RR")

    def test_annotation_filter_used_when_no_extension(
        self, tmp_path: Path
    ) -> None:
        args = self._project_args(["annotation=G026"])
        code, m = self._run(tmp_path, args)
        _, kwargs = m.call_args
        assert kwargs["annotations"] == "G026"

    def test_annotation_extension_overrides_annotation(
        self, tmp_path: Path
    ) -> None:
        args = self._project_args(
            ["annotation=default", "annotation_extension=G029,G030"]
        )
        code, m = self._run(tmp_path, args)
        _, kwargs = m.call_args
        assert kwargs["annotations"] == ("G029", "G030")

    def test_junction_type_default(self, tmp_path: Path) -> None:
        args = self._project_args()
        code, m = self._run(tmp_path, args)
        _, kwargs = m.call_args
        assert kwargs["junction_type"] == "ALL"


class TestCmdSearchUnknownMode:
    def test_unknown_mode_raises_value_error(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        args = _make_search_args(mode="bogus", filters=[])
        with pytest.raises(ValueError, match="Unknown search mode"):
            _cmd_search(args, cfg)


class TestCmdSearchOutputDestination:
    def test_output_flag_writes_to_file(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource()
        out_file = tmp_path / "out.jsonl"
        args = _make_search_args(
            mode="annotations",
            filters=[
                "organism=human",
                "genomic_unit=gene",
                "annotation_extension=G026",
            ],
            output=str(out_file),
        )
        with mock.patch(
            "recount3.cli.r3_search.search_annotations", return_value=[res]
        ):
            code = _cmd_search(args, cfg)
        assert code == 0
        assert out_file.exists()

    def test_outdir_creates_timestamped_file(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource()
        outdir = tmp_path / "outdir"
        outdir.mkdir()
        args = _make_search_args(
            mode="annotations",
            filters=[
                "organism=human",
                "genomic_unit=gene",
                "annotation_extension=G026",
            ],
            outdir=str(outdir),
        )
        with mock.patch(
            "recount3.cli.r3_search.search_annotations", return_value=[res]
        ):
            code = _cmd_search(args, cfg)
        assert code == 0
        files = list(outdir.iterdir())
        assert len(files) == 1
        assert files[0].suffix == ".jsonl"

    def test_tsv_format(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource()
        args = _make_search_args(
            mode="annotations",
            filters=[
                "organism=human",
                "genomic_unit=gene",
                "annotation_extension=G026",
            ],
            format="tsv",
        )
        with mock.patch(
            "recount3.cli.r3_search.search_annotations", return_value=[res]
        ):
            code = _cmd_search(args, cfg)
        assert code == 0
        captured = capsys.readouterr()
        assert "resource_type" in captured.out  # TSV header

    def test_outdir_tsv_creates_tsv_file(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource()
        outdir = tmp_path / "outdir"
        outdir.mkdir()
        args = _make_search_args(
            mode="annotations",
            filters=[
                "organism=human",
                "genomic_unit=gene",
                "annotation_extension=G026",
            ],
            outdir=str(outdir),
            format="tsv",
        )
        with mock.patch(
            "recount3.cli.r3_search.search_annotations", return_value=[res]
        ):
            _cmd_search(args, cfg)
        files = list(outdir.iterdir())
        assert files[0].suffix == ".tsv"

    def test_stdout_logging_when_out_path_none(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource()
        args = _make_search_args(
            mode="annotations",
            filters=[
                "organism=human",
                "genomic_unit=gene",
                "annotation_extension=G026",
            ],
        )
        with mock.patch(
            "recount3.cli.r3_search.search_annotations", return_value=[res]
        ):
            with caplog.at_level(logging.INFO, logger="recount3"):
                _cmd_search(args, cfg)

    def test_file_logging_when_out_path_set(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource()
        out_file = tmp_path / "out.jsonl"
        args = _make_search_args(
            mode="annotations",
            filters=[
                "organism=human",
                "genomic_unit=gene",
                "annotation_extension=G026",
            ],
            output=str(out_file),
        )
        with mock.patch(
            "recount3.cli.r3_search.search_annotations", return_value=[res]
        ):
            code = _cmd_search(args, cfg)
        assert code == 0


from recount3.cli import _download_one


class TestDownloadOne:
    def test_zip_dest_downloads_and_returns_ok(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        dest = tmp_path / "out.zip"
        with mock.patch.object(R3Resource, "download", return_value=None) as m:
            evt = _download_one(res, cfg, dest, "enable", False)
        assert evt["status"] == "ok"
        m.assert_called_once()

    def test_dir_dest_file_not_exists_downloads(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        dest = tmp_path / "downloads"
        dest.mkdir()
        expected_file = dest / Path(res.description.url_path()).name
        assert not expected_file.exists()
        with mock.patch.object(
            R3Resource, "download", return_value=str(expected_file)
        ) as m:
            evt = _download_one(res, cfg, dest, "enable", False)
        assert evt["status"] == "ok"

    def test_dir_dest_file_exists_no_overwrite_cache_enabled_skips(
        self, tmp_path: Path
    ) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        dest = tmp_path / "downloads"
        dest.mkdir()
        dest_file = dest / Path(res.description.url_path()).name
        dest_file.touch()
        with mock.patch.object(R3Resource, "download", return_value=None) as m:
            evt = _download_one(res, cfg, dest, "enable", False)
        assert evt["status"] == "skipped"
        m.assert_called_once_with(path=None, cache_mode="enable")

    def test_dir_dest_file_exists_no_overwrite_cache_disable_skips_no_download(
        self, tmp_path: Path
    ) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        dest = tmp_path / "downloads"
        dest.mkdir()
        dest_file = dest / Path(res.description.url_path()).name
        dest_file.touch()
        with mock.patch.object(R3Resource, "download", return_value=None) as m:
            evt = _download_one(res, cfg, dest, "disable", False)
        assert evt["status"] == "skipped"
        m.assert_not_called()

    def test_exception_returns_error_status(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        dest = tmp_path / "downloads"
        dest.mkdir()
        with mock.patch.object(
            R3Resource, "download", side_effect=OSError("network failure")
        ):
            evt = _download_one(res, cfg, dest, "enable", False)
        assert evt["status"] == "error"
        assert "error" in evt
        assert "OSError" in evt["error"]

    def test_overwrite_true_does_not_skip(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        dest = tmp_path / "downloads"
        dest.mkdir()
        dest_file = dest / Path(res.description.url_path()).name
        dest_file.touch()
        with mock.patch.object(
            R3Resource, "download", return_value=str(dest_file)
        ) as m:
            evt = _download_one(res, cfg, dest, "enable", True)
        assert evt["status"] == "ok"


class TestCmdDownload:
    def _make_args(self, **kw: Any) -> argparse.Namespace:
        defaults = dict(
            inline=None,
            manifest=None,
            dest=".",
            overwrite=False,
            jobs=2,
            cache="enable",
        )
        defaults.update(kw)
        return argparse.Namespace(**defaults)

    def test_inline_valid_json_single_resource(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        inline_obj = {
            "resource_type": "annotations",
            "organism": "human",
            "genomic_unit": "gene",
            "annotation_extension": "G026",
        }
        args = self._make_args(
            inline=json.dumps(inline_obj), dest=str(tmp_path)
        )
        with mock.patch("recount3.cli._download_one") as m:
            m.return_value = {
                "url": "http://x",
                "status": "ok",
                "dest": "/tmp/f",
            }
            code = _cmd_download(args, cfg)
        assert code == 0
        captured = capsys.readouterr()
        evt = json.loads(captured.out.strip())
        assert evt["status"] == "ok"

    def test_inline_bad_json_returns_1(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        args = self._make_args(inline="not valid json!!!")
        code = _cmd_download(args, cfg)
        assert code == 1

    def test_manifest_reads_resources(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        args = self._make_args(manifest="manifest.jsonl", dest=str(tmp_path))
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch("recount3.cli._download_one") as m,
        ):
            m.return_value = {"url": "http://x", "status": "ok", "dest": "/f"}
            code = _cmd_download(args, cfg)
        assert code == 0

    def test_no_resources_returns_0_with_warning(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        args = self._make_args(manifest="m.jsonl", dest=str(tmp_path))
        with mock.patch("recount3.cli._iter_manifest", return_value=[]):
            code = _cmd_download(args, cfg)
        assert code == 0

    def test_zip_dest_creates_parent(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        zip_dest = tmp_path / "subdir" / "out.zip"
        args = self._make_args(manifest="m.jsonl", dest=str(zip_dest))
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch("recount3.cli._download_one") as m,
        ):
            m.return_value = {"url": "http://x", "status": "ok", "dest": None}
            code = _cmd_download(args, cfg)
        assert code == 0
        assert (tmp_path / "subdir").is_dir()

    def test_directory_dest_is_created(self, tmp_path: Path) -> None:
        """A non-existent directory --dest should be created before
        per-resource writes run; otherwise every write races on
        FileNotFoundError (see resource.py case-3 materialization)."""
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        dir_dest = tmp_path / "fresh" / "downloads"
        assert not dir_dest.exists()
        args = self._make_args(manifest="m.jsonl", dest=str(dir_dest))
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch("recount3.cli._download_one") as m,
        ):
            m.return_value = {"url": "http://x", "status": "ok", "dest": "/f"}
            code = _cmd_download(args, cfg)
        assert code == 0
        assert dir_dest.is_dir()

    def test_all_errors_returns_2(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        args = self._make_args(manifest="m.jsonl", dest=str(tmp_path))
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch("recount3.cli._download_one") as m,
        ):
            m.return_value = {
                "url": "http://x",
                "status": "error",
                "dest": None,
                "error": "oops",
            }
            code = _cmd_download(args, cfg)
        assert code == 2

    def test_partial_errors_returns_3(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        res1 = _make_annotation_resource(cfg)
        res2 = _make_annotation_resource(cfg)
        args = self._make_args(manifest="m.jsonl", dest=str(tmp_path), jobs=1)
        results = [
            {"url": "http://a", "status": "ok", "dest": "/f1"},
            {
                "url": "http://b",
                "status": "error",
                "dest": None,
                "error": "oops",
            },
        ]
        with (
            mock.patch(
                "recount3.cli._iter_manifest", return_value=[res1, res2]
            ),
            mock.patch("recount3.cli._download_one", side_effect=results),
        ):
            code = _cmd_download(args, cfg)
        assert code == 3


class TestCmdBundleStackCounts:
    def _make_args(self, out: str, **kw: Any) -> argparse.Namespace:
        defaults = dict(
            manifest="-",
            compat="family",
            join="inner",
            axis=1,
            verify_integrity=False,
            densify=False,
            out=out,
        )
        defaults.update(kw)
        return argparse.Namespace(**defaults)

    def test_success_parquet(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.parquet"
        args = self._make_args(str(out))
        mock_df = mock.MagicMock()
        mock_df.dtypes.items.return_value = []
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.return_value = mock_df
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch(
                "recount3._utils.ensure_parquet_engine", return_value="pyarrow"
            ),
        ):
            code = _cmd_bundle_stack_counts(args, cfg)
        assert code == 0
        mock_df.to_parquet.assert_called_once_with(out)

    def test_sparse_parquet_requires_explicit_densify(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.parquet"
        args = self._make_args(str(out))
        mock_df = mock.MagicMock()
        mock_df.shape = (5, 2)
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.return_value = mock_df
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch(
                "recount3._utils.ensure_parquet_engine",
                return_value="pyarrow",
            ),
            mock.patch(
                "recount3._utils.sparse_column_names",
                return_value=["sample/a", "sample/b"],
            ),
            caplog.at_level(logging.ERROR),
        ):
            code = _cmd_bundle_stack_counts(args, cfg)

        assert code == 2
        mock_df.to_parquet.assert_not_called()
        assert "2 of 2 columns use a pandas sparse dtype" in caplog.text
        assert "--densify" in caplog.text

    def test_sparse_parquet_densifies_when_requested(
        self, tmp_path: Path
    ) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.parquet"
        args = self._make_args(str(out), densify=True)
        sparse_df = mock.MagicMock()
        sparse_df.shape = (5, 2)
        dense_df = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.return_value = sparse_df
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch(
                "recount3._utils.ensure_parquet_engine",
                return_value="pyarrow",
            ),
            mock.patch(
                "recount3._utils.sparse_column_names",
                return_value=["sample/a", "sample/b"],
            ),
            mock.patch(
                "recount3._utils.densify_sparse_columns",
                return_value=dense_df,
            ) as mock_densify,
        ):
            code = _cmd_bundle_stack_counts(args, cfg)

        assert code == 0
        mock_densify.assert_called_once_with(sparse_df)
        dense_df.to_parquet.assert_called_once_with(out)

    def test_success_tsv(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.tsv"
        args = self._make_args(str(out))
        mock_df = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.return_value = mock_df
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
        ):
            code = _cmd_bundle_stack_counts(args, cfg)
        assert code == 0
        mock_df.to_csv.assert_called_once_with(out, sep="\t")

    def test_success_csv_uses_comma_separator(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.csv"
        args = self._make_args(str(out))
        mock_df = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.return_value = mock_df
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
        ):
            code = _cmd_bundle_stack_counts(args, cfg)
        assert code == 0
        mock_df.to_csv.assert_called_once_with(out, sep=",")

    def test_compatibility_error_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.tsv"
        args = self._make_args(str(out))
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.side_effect = CompatibilityError(
            "incompat"
        )
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
        ):
            code = _cmd_bundle_stack_counts(args, cfg)
        assert code == 2

    def test_load_error_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.tsv"
        args = self._make_args(str(out))
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.side_effect = LoadError("bad")
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
        ):
            code = _cmd_bundle_stack_counts(args, cfg)
        assert code == 2

    def test_value_error_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.tsv"
        args = self._make_args(str(out))
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.side_effect = ValueError("bad")
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
        ):
            code = _cmd_bundle_stack_counts(args, cfg)
        assert code == 2

    def test_import_error_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.parquet"
        args = self._make_args(str(out))
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.side_effect = ImportError(
            "missing dep"
        )
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch(
                "recount3._utils.ensure_parquet_engine", return_value="pyarrow"
            ),
        ):
            code = _cmd_bundle_stack_counts(args, cfg)
        assert code == 2
        mock_bundle.stack_count_matrices.assert_called_once()

    def test_write_exception_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.parquet"
        args = self._make_args(str(out))
        mock_df = mock.MagicMock()
        mock_df.dtypes.items.return_value = []
        mock_df.to_parquet.side_effect = OSError("disk full")
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.return_value = mock_df
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch(
                "recount3._utils.ensure_parquet_engine", return_value="pyarrow"
            ),
        ):
            code = _cmd_bundle_stack_counts(args, cfg)
        assert code == 2
        mock_df.to_parquet.assert_called_once_with(out)

    _GENE_IDS = [
        "ENSG00000000003.14",
        "ENSG00000000005.5",
        "ENSG00000000419.12",
    ]
    _EXON_IDS = [
        "ENSG00000000003.14|1",
        "ENSG00000000003.14|2",
        "ENSG00000000419.12|1",
    ]
    _SAMPLES = ["SRR387777", "SRR387778", "SRR387779"]

    def _run_with_frame(
        self,
        frame: pd.DataFrame,
        args: argparse.Namespace,
        cfg: Config,
    ) -> int:
        """Run the command against a real DataFrame instead of a mock."""
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.return_value = frame
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
        ):
            return _cmd_bundle_stack_counts(args, cfg)

    @staticmethod
    def _counts_frame(feature_ids: list[str]) -> pd.DataFrame:
        values = np.array(
            [[10, 0, 7], [0, 3, 0], [5, 5, 5]],
            dtype=np.int64,
        )
        frame = pd.DataFrame(
            values,
            index=pd.Index(feature_ids, name="feature_id"),
            columns=TestCmdBundleStackCounts._SAMPLES,
        )
        return frame

    @pytest.mark.requires_parquet
    @pytest.mark.parametrize(
        ("unit", "feature_ids"),
        [
            ("gene", _GENE_IDS),
            ("exon", _EXON_IDS),
        ],
    )
    def test_parquet_round_trip(
        self, tmp_path: Path, unit: str, feature_ids: list[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / f"{unit}_counts.parquet"
        frame = self._counts_frame(feature_ids)

        code = self._run_with_frame(frame, self._make_args(str(out)), cfg)

        assert code == 0
        assert out.is_file() and out.stat().st_size > 0

        back = pd.read_parquet(out)
        assert list(back.index) == feature_ids
        assert back.index.name == "feature_id"
        assert list(back.columns) == self._SAMPLES
        np.testing.assert_array_equal(back.to_numpy(), frame.to_numpy())
        pd.testing.assert_frame_equal(back, frame)

    @pytest.mark.requires_parquet
    def test_parquet_round_trip_preserves_sample_lookup(
        self, tmp_path: Path
    ) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "counts.parquet"
        frame = self._counts_frame(self._GENE_IDS)

        code = self._run_with_frame(frame, self._make_args(str(out)), cfg)
        assert code == 0

        back = pd.read_parquet(out)
        assert back.loc["ENSG00000000005.5", "SRR387778"] == 3
        assert back["SRR387777"].sum() == 15

    def test_missing_parquet_engine_returns_2_before_any_download(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.parquet"
        args = self._make_args(str(out))
        boom = ImportError(
            "Writing Parquet requires a Parquet engine.\n\n"
            'Install one with:\n\n  pip install "recount3[parquet]"\n'
        )
        with (
            mock.patch(
                "recount3._utils.ensure_parquet_engine", side_effect=boom
            ),
            mock.patch("recount3.cli._iter_manifest") as mock_manifest,
            mock.patch("recount3.cli.R3ResourceBundle") as mock_bundle_cls,
            caplog.at_level(logging.ERROR),
        ):
            code = _cmd_bundle_stack_counts(args, cfg)

        assert code == 2
        # Nothing was read, loaded, or downloaded.
        mock_manifest.assert_not_called()
        mock_bundle_cls.assert_not_called()
        assert not out.exists()
        assert 'pip install "recount3[parquet]"' in caplog.text

    def test_text_output_does_not_require_a_parquet_engine(
        self, tmp_path: Path
    ) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.tsv"
        frame = self._counts_frame(self._GENE_IDS)
        with mock.patch(
            "recount3._utils.ensure_parquet_engine",
            side_effect=AssertionError("engine checked for text output"),
        ):
            code = self._run_with_frame(frame, self._make_args(str(out)), cfg)

        assert code == 0
        back = pd.read_csv(out, sep="\t", index_col=0)
        assert list(back.columns) == self._SAMPLES

    @staticmethod
    def _sparse_junction_frame() -> pd.DataFrame:
        """Mimic a junction MM load: a sparse-backed count matrix."""
        matrix = scipy.sparse.csr_array(
            np.array([[3, 0, 0], [0, 0, 11], [0, 5, 0]], dtype=np.int64)
        )
        frame = pd.DataFrame.sparse.from_spmatrix(matrix)
        frame.columns = TestCmdBundleStackCounts._SAMPLES
        frame.index = pd.Index(["0", "1", "2"], name="junction_id")
        return frame

    @pytest.mark.requires_parquet
    def test_sparse_junctions_to_parquet_without_densify_returns_2(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "junctions.parquet"
        frame = self._sparse_junction_frame()

        with caplog.at_level(logging.ERROR):
            code = self._run_with_frame(frame, self._make_args(str(out)), cfg)

        assert code == 2
        assert not out.exists()
        assert "sparse dtype" in caplog.text
        assert "--densify" in caplog.text

    @pytest.mark.requires_parquet
    def test_sparse_junctions_to_parquet_with_densify_round_trips(
        self, tmp_path: Path
    ) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "junctions.parquet"
        frame = self._sparse_junction_frame()
        args = self._make_args(str(out), densify=True)

        code = self._run_with_frame(frame, args, cfg)

        assert code == 0
        back = pd.read_parquet(out)
        assert list(back.index) == ["0", "1", "2"]
        assert back.index.name == "junction_id"
        assert list(back.columns) == self._SAMPLES
        assert not any(
            isinstance(dtype, pd.SparseDtype) for dtype in back.dtypes
        )
        np.testing.assert_array_equal(
            back.to_numpy(),
            np.array([[3, 0, 0], [0, 0, 11], [0, 5, 0]], dtype=np.int64),
        )

    def test_sparse_junctions_to_tsv_needs_no_densify(
        self, tmp_path: Path
    ) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "junctions.tsv"
        frame = self._sparse_junction_frame()

        code = self._run_with_frame(frame, self._make_args(str(out)), cfg)

        assert code == 0
        back = pd.read_csv(out, sep="\t", index_col=0)
        assert list(back.columns) == self._SAMPLES
        np.testing.assert_array_equal(
            back.to_numpy(),
            np.array([[3, 0, 0], [0, 0, 11], [0, 5, 0]], dtype=np.int64),
        )


def _make_experiment(
    *,
    counts: np.ndarray,
    feature_ids: list[str],
    samples: list[str],
    biocframe: Any,
    extra_row_data: dict[str, list[Any]] | None = None,
    extra_column_data: dict[str, list[Any]] | None = None,
) -> Any:
    """Build a real SummarizedExperiment for AnnData round-trip tests."""
    import summarizedexperiment  # pylint: disable=import-outside-toplevel

    row_data: dict[str, list[Any]] = {"feature_id": list(feature_ids)}
    row_data.update(extra_row_data or {})
    column_data: dict[str, list[Any]] = {"sample": list(samples)}
    column_data.update(extra_column_data or {})

    return summarizedexperiment.SummarizedExperiment(
        assays={"raw_counts": counts},
        row_data=biocframe.BiocFrame(row_data, row_names=feature_ids),
        column_data=biocframe.BiocFrame(column_data, row_names=samples),
    )


class TestCmdBundleSe:
    def _make_args(self, out: str, **kw: Any) -> argparse.Namespace:
        defaults = dict(
            manifest="-",
            genomic_unit="gene",
            annotation=None,
            assay_name="raw_counts",
            join="inner",
            sanitize_columns=False,
            out=out,
        )
        defaults.update(kw)
        return argparse.Namespace(**defaults)

    def test_success_pkl(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.pkl"
        args = self._make_args(str(out))
        mock_se = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.return_value = mock_se
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch("pickle.dump") as mock_pickle,
        ):
            code = _cmd_bundle_se(args, cfg)
        assert code == 0
        mock_pickle.assert_called_once()

    def test_success_h5ad(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.h5ad"
        args = self._make_args(str(out))
        mock_se = mock.MagicMock()
        mock_adata = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.return_value = mock_se
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch("recount3._utils.ensure_anndata_support"),
            mock.patch(
                "recount3._utils.experiment_to_anndata",
                return_value=mock_adata,
            ),
        ):
            code = _cmd_bundle_se(args, cfg)
        assert code == 0
        mock_adata.write_h5ad.assert_called_once_with(out)

    def test_h5ad_unsafe_columns_require_sanitization_flag(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.h5ad"
        args = self._make_args(str(out))
        mock_se = mock.MagicMock()
        mock_adata = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.return_value = mock_se
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch("recount3._utils.ensure_anndata_support"),
            mock.patch(
                "recount3._utils.experiment_to_anndata",
                return_value=mock_adata,
            ),
            mock.patch(
                "recount3._utils.normalize_anndata_for_hdf5",
                return_value=[],
            ),
            mock.patch(
                "recount3._utils.hdf5_unsafe_column_names",
                return_value=["star/metric"],
            ),
            mock.patch(
                "recount3._utils.sanitize_anndata_column_names"
            ) as mock_sanitize,
            caplog.at_level(logging.ERROR),
        ):
            code = _cmd_bundle_se(args, cfg)

        assert code == 2
        mock_sanitize.assert_not_called()
        mock_adata.write_h5ad.assert_not_called()
        assert "--sanitize-columns" in caplog.text

    def test_h5ad_normalizes_and_sanitizes_columns_with_flag(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.h5ad"
        args = self._make_args(str(out), sanitize_columns=True)
        mock_se = mock.MagicMock()
        mock_adata = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.return_value = mock_se
        res = _make_annotation_resource(cfg)
        renames = [("star/metric", "star_metric")]
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch("recount3._utils.ensure_anndata_support"),
            mock.patch(
                "recount3._utils.experiment_to_anndata",
                return_value=mock_adata,
            ),
            mock.patch(
                "recount3._utils.normalize_anndata_for_hdf5",
                return_value=["all_missing"],
            ),
            mock.patch(
                "recount3._utils.hdf5_unsafe_column_names",
                return_value=["star/metric"],
            ),
            mock.patch(
                "recount3._utils.sanitize_anndata_column_names",
                return_value=renames,
            ) as mock_sanitize,
            caplog.at_level(logging.INFO),
        ):
            code = _cmd_bundle_se(args, cfg)

        assert code == 0
        mock_sanitize.assert_called_once_with(mock_adata)
        mock_adata.write_h5ad.assert_called_once_with(out)
        assert "Cast 1 all-missing column" in caplog.text
        assert "star/metric -> star_metric" in caplog.text

    def test_import_error_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.pkl"
        args = self._make_args(str(out))
        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.side_effect = ImportError(
            "no se dep"
        )
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
        ):
            code = _cmd_bundle_se(args, cfg)
        assert code == 2

    def test_general_exception_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.pkl"
        args = self._make_args(str(out))
        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.side_effect = RuntimeError("oops")
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
        ):
            code = _cmd_bundle_se(args, cfg)
        assert code == 2

    def test_write_exception_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.pkl"
        args = self._make_args(str(out))
        mock_se = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.return_value = mock_se
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch("pickle.dump", side_effect=OSError("disk full")),
        ):
            code = _cmd_bundle_se(args, cfg)
        assert code == 2

    def test_h5ad_write_exception_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.h5ad"
        args = self._make_args(str(out))
        mock_se = mock.MagicMock()
        mock_adata = mock.MagicMock()
        mock_adata.write_h5ad.side_effect = OSError("disk full")
        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.return_value = mock_se
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch("recount3._utils.ensure_anndata_support"),
            mock.patch(
                "recount3._utils.experiment_to_anndata",
                return_value=mock_adata,
            ),
        ):
            code = _cmd_bundle_se(args, cfg)
        assert code == 2
        mock_adata.write_h5ad.assert_called_once_with(out)

    def test_h5ad_missing_anndata_returns_2_before_any_build(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.h5ad"
        args = self._make_args(str(out))
        boom = ImportError(
            "Optional dependency 'anndata' is required for this feature.\n\n"
            'Install it with:\n\n  pip install "recount3[anndata]"\n'
        )
        with (
            mock.patch(
                "recount3._utils.ensure_anndata_support", side_effect=boom
            ),
            mock.patch("recount3.cli._iter_manifest") as mock_manifest,
            mock.patch("recount3.cli.R3ResourceBundle") as mock_bundle_cls,
            caplog.at_level(logging.ERROR),
        ):
            code = _cmd_bundle_se(args, cfg)

        assert code == 2
        mock_manifest.assert_not_called()
        mock_bundle_cls.assert_not_called()
        assert not out.exists()
        assert 'pip install "recount3[anndata]"' in caplog.text

    def test_pkl_output_does_not_require_anndata(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.pkl"
        args = self._make_args(str(out))
        mock_obj = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.return_value = mock_obj
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch(
                "recount3._utils.ensure_anndata_support",
                side_effect=AssertionError("anndata checked for .pkl output"),
            ),
            mock.patch("pickle.dump"),
        ):
            code = _cmd_bundle_se(args, cfg)

        assert code == 0

    @pytest.mark.requires_anndata
    def test_h5ad_round_trip(self, tmp_path: Path) -> None:
        anndata = pytest.importorskip("anndata")
        biocframe = pytest.importorskip("biocframe")

        cfg = _make_cfg(tmp_path)
        out = tmp_path / "counts.h5ad"
        args = self._make_args(str(out))

        feature_ids = ["ENSG00000000003.14", "ENSG00000000005.5"]
        samples = ["SRR387777", "SRR387778", "SRR387779"]
        counts = np.array([[10, 0, 7], [0, 3, 0]], dtype=np.int64)
        experiment = _make_experiment(
            counts=counts,
            feature_ids=feature_ids,
            samples=samples,
            biocframe=biocframe,
        )

        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.return_value = experiment
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
        ):
            code = _cmd_bundle_se(args, cfg)

        assert code == 0
        assert out.is_file() and out.stat().st_size > 0

        back = anndata.read_h5ad(out)
        # AnnData is observation-major: samples become obs, features var.
        assert list(back.obs_names) == samples
        assert list(back.var_names) == feature_ids
        np.testing.assert_array_equal(
            np.asarray(back.layers["raw_counts"]), counts.T
        )

    @pytest.mark.requires_anndata
    def test_h5ad_all_missing_column_is_cast_to_nan(
        self, tmp_path: Path
    ) -> None:
        """An absent-for-every-row field would otherwise fail h5py."""
        anndata = pytest.importorskip("anndata")
        biocframe = pytest.importorskip("biocframe")

        cfg = _make_cfg(tmp_path)
        out = tmp_path / "counts.h5ad"
        args = self._make_args(str(out))
        experiment = _make_experiment(
            counts=np.array([[1, 2], [3, 4]], dtype=np.int64),
            feature_ids=["g1", "g2"],
            samples=["s1", "s2"],
            biocframe=biocframe,
            extra_row_data={"phase": [None, None]},
        )

        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.return_value = experiment
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
        ):
            code = _cmd_bundle_se(args, cfg)

        assert code == 0
        back = anndata.read_h5ad(out)
        assert back.var["phase"].isna().all()

    @pytest.mark.requires_anndata
    def test_h5ad_slash_in_column_name_returns_2_without_flag(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """HDF5 reads '/' as a path separator; renaming needs consent."""
        biocframe = pytest.importorskip("biocframe")

        cfg = _make_cfg(tmp_path)
        out = tmp_path / "counts.h5ad"
        args = self._make_args(str(out))
        experiment = _make_experiment(
            counts=np.array([[1, 2], [3, 4]], dtype=np.int64),
            feature_ids=["g1", "g2"],
            samples=["s1", "s2"],
            biocframe=biocframe,
            extra_column_data={
                "recount_qc__star.number_of_splices:_gt/ag": ["1", "2"]
            },
        )

        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.return_value = experiment
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            caplog.at_level(logging.ERROR),
        ):
            code = _cmd_bundle_se(args, cfg)

        assert code == 2
        assert not out.exists()
        assert "--sanitize-columns" in caplog.text
        assert "forward slash" in caplog.text

    @pytest.mark.requires_anndata
    def test_h5ad_slash_in_column_name_round_trips_with_flag(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        anndata = pytest.importorskip("anndata")
        biocframe = pytest.importorskip("biocframe")

        cfg = _make_cfg(tmp_path)
        out = tmp_path / "counts.h5ad"
        args = self._make_args(str(out), sanitize_columns=True)
        experiment = _make_experiment(
            counts=np.array([[1, 2], [3, 4]], dtype=np.int64),
            feature_ids=["g1", "g2"],
            samples=["s1", "s2"],
            biocframe=biocframe,
            extra_column_data={
                "recount_qc__star.number_of_splices:_gt/ag": ["1", "2"]
            },
        )

        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.return_value = experiment
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            caplog.at_level(logging.WARNING),
        ):
            code = _cmd_bundle_se(args, cfg)

        assert code == 0
        assert "_gt/ag -> " in caplog.text

        back = anndata.read_h5ad(out)
        assert "recount_qc__star.number_of_splices:_gt_ag" in back.obs.columns
        assert "recount_qc__star.number_of_splices:_gt/ag" not in back.obs


class TestCmdBundleRse:
    def _make_args(self, out: str, **kw: Any) -> argparse.Namespace:
        defaults = dict(
            manifest="-",
            genomic_unit="gene",
            annotation=None,
            assay_name="raw_counts",
            join="inner",
            allow_fallback_to_se=False,
            sanitize_columns=False,
            out=out,
        )
        defaults.update(kw)
        return argparse.Namespace(**defaults)

    def test_success_pkl(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.pkl"
        args = self._make_args(str(out))
        mock_rse = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.return_value = mock_rse
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch("pickle.dump") as mock_pickle,
        ):
            code = _cmd_bundle_rse(args, cfg)
        assert code == 0
        mock_pickle.assert_called_once()

    def test_success_h5ad(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.h5ad"
        args = self._make_args(str(out))
        mock_rse = mock.MagicMock()
        mock_adata = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.return_value = mock_rse
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch("recount3._utils.ensure_anndata_support"),
            mock.patch(
                "recount3._utils.experiment_to_anndata",
                return_value=mock_adata,
            ),
        ):
            code = _cmd_bundle_rse(args, cfg)
        assert code == 0
        mock_adata.write_h5ad.assert_called_once_with(out)

    def test_h5ad_unsafe_columns_require_sanitization_flag(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.h5ad"
        args = self._make_args(str(out))
        mock_rse = mock.MagicMock()
        mock_adata = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.return_value = mock_rse
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch("recount3._utils.ensure_anndata_support"),
            mock.patch(
                "recount3._utils.experiment_to_anndata",
                return_value=mock_adata,
            ),
            mock.patch(
                "recount3._utils.normalize_anndata_for_hdf5",
                return_value=[],
            ),
            mock.patch(
                "recount3._utils.hdf5_unsafe_column_names",
                return_value=["star/metric"],
            ),
            mock.patch(
                "recount3._utils.sanitize_anndata_column_names"
            ) as mock_sanitize,
            caplog.at_level(logging.ERROR),
        ):
            code = _cmd_bundle_rse(args, cfg)

        assert code == 2
        mock_sanitize.assert_not_called()
        mock_adata.write_h5ad.assert_not_called()
        assert "--sanitize-columns" in caplog.text

    def test_h5ad_normalizes_and_sanitizes_columns_with_flag(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.h5ad"
        args = self._make_args(str(out), sanitize_columns=True)
        mock_rse = mock.MagicMock()
        mock_adata = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.return_value = mock_rse
        res = _make_annotation_resource(cfg)
        renames = [("star/metric", "star_metric")]
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch("recount3._utils.ensure_anndata_support"),
            mock.patch(
                "recount3._utils.experiment_to_anndata",
                return_value=mock_adata,
            ),
            mock.patch(
                "recount3._utils.normalize_anndata_for_hdf5",
                return_value=["all_missing"],
            ),
            mock.patch(
                "recount3._utils.hdf5_unsafe_column_names",
                return_value=["star/metric"],
            ),
            mock.patch(
                "recount3._utils.sanitize_anndata_column_names",
                return_value=renames,
            ) as mock_sanitize,
            caplog.at_level(logging.INFO),
        ):
            code = _cmd_bundle_rse(args, cfg)

        assert code == 0
        mock_sanitize.assert_called_once_with(mock_adata)
        mock_adata.write_h5ad.assert_called_once_with(out)
        assert "Cast 1 all-missing column" in caplog.text
        assert "star/metric -> star_metric" in caplog.text

    def test_allow_fallback_to_se_passed_through(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.pkl"
        args = self._make_args(str(out), allow_fallback_to_se=True)
        mock_rse = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.return_value = mock_rse
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch("pickle.dump"),
        ):
            code = _cmd_bundle_rse(args, cfg)
        assert code == 0
        _, kwargs = mock_bundle.to_ranged_summarized_experiment.call_args
        assert kwargs["allow_fallback_to_se"] is True

    def test_import_error_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.pkl"
        args = self._make_args(str(out))
        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.side_effect = ImportError(
            "missing"
        )
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
        ):
            code = _cmd_bundle_rse(args, cfg)
        assert code == 2

    def test_general_exception_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.pkl"
        args = self._make_args(str(out))
        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.side_effect = RuntimeError(
            "bad"
        )
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
        ):
            code = _cmd_bundle_rse(args, cfg)
        assert code == 2

    def test_write_exception_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.pkl"
        args = self._make_args(str(out))
        mock_rse = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.return_value = mock_rse
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch("pickle.dump", side_effect=OSError("disk full")),
        ):
            code = _cmd_bundle_rse(args, cfg)
        assert code == 2

    def test_h5ad_write_exception_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.h5ad"
        args = self._make_args(str(out))
        mock_rse = mock.MagicMock()
        mock_adata = mock.MagicMock()
        mock_adata.write_h5ad.side_effect = OSError("disk full")
        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.return_value = mock_rse
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch("recount3._utils.ensure_anndata_support"),
            mock.patch(
                "recount3._utils.experiment_to_anndata",
                return_value=mock_adata,
            ),
        ):
            code = _cmd_bundle_rse(args, cfg)
        assert code == 2
        mock_adata.write_h5ad.assert_called_once_with(out)

    def test_h5ad_missing_anndata_returns_2_before_any_build(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.h5ad"
        args = self._make_args(str(out))
        boom = ImportError(
            "Optional dependency 'anndata' is required for this feature.\n\n"
            'Install it with:\n\n  pip install "recount3[anndata]"\n'
        )
        with (
            mock.patch(
                "recount3._utils.ensure_anndata_support", side_effect=boom
            ),
            mock.patch("recount3.cli._iter_manifest") as mock_manifest,
            mock.patch("recount3.cli.R3ResourceBundle") as mock_bundle_cls,
            caplog.at_level(logging.ERROR),
        ):
            code = _cmd_bundle_rse(args, cfg)

        assert code == 2
        mock_manifest.assert_not_called()
        mock_bundle_cls.assert_not_called()
        assert not out.exists()
        assert 'pip install "recount3[anndata]"' in caplog.text

    def test_pkl_output_does_not_require_anndata(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.pkl"
        args = self._make_args(str(out))
        mock_obj = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.return_value = mock_obj
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            mock.patch(
                "recount3._utils.ensure_anndata_support",
                side_effect=AssertionError("anndata checked for .pkl output"),
            ),
            mock.patch("pickle.dump"),
        ):
            code = _cmd_bundle_rse(args, cfg)

        assert code == 0

    @pytest.mark.requires_anndata
    def test_h5ad_round_trip(self, tmp_path: Path) -> None:
        anndata = pytest.importorskip("anndata")
        biocframe = pytest.importorskip("biocframe")

        cfg = _make_cfg(tmp_path)
        out = tmp_path / "counts.h5ad"
        args = self._make_args(str(out))

        feature_ids = ["ENSG00000000003.14", "ENSG00000000005.5"]
        samples = ["SRR387777", "SRR387778", "SRR387779"]
        counts = np.array([[10, 0, 7], [0, 3, 0]], dtype=np.int64)
        experiment = _make_experiment(
            counts=counts,
            feature_ids=feature_ids,
            samples=samples,
            biocframe=biocframe,
        )

        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.return_value = experiment
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
        ):
            code = _cmd_bundle_rse(args, cfg)

        assert code == 0
        assert out.is_file() and out.stat().st_size > 0

        back = anndata.read_h5ad(out)
        # AnnData is observation-major: samples become obs, features var.
        assert list(back.obs_names) == samples
        assert list(back.var_names) == feature_ids
        np.testing.assert_array_equal(
            np.asarray(back.layers["raw_counts"]), counts.T
        )

    @pytest.mark.requires_anndata
    def test_h5ad_all_missing_column_is_cast_to_nan(
        self, tmp_path: Path
    ) -> None:
        """An absent-for-every-row field would otherwise fail h5py."""
        anndata = pytest.importorskip("anndata")
        biocframe = pytest.importorskip("biocframe")

        cfg = _make_cfg(tmp_path)
        out = tmp_path / "counts.h5ad"
        args = self._make_args(str(out))
        experiment = _make_experiment(
            counts=np.array([[1, 2], [3, 4]], dtype=np.int64),
            feature_ids=["g1", "g2"],
            samples=["s1", "s2"],
            biocframe=biocframe,
            extra_row_data={"phase": [None, None]},
        )

        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.return_value = experiment
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
        ):
            code = _cmd_bundle_rse(args, cfg)

        assert code == 0
        back = anndata.read_h5ad(out)
        assert back.var["phase"].isna().all()

    @pytest.mark.requires_anndata
    def test_h5ad_slash_in_column_name_returns_2_without_flag(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """HDF5 reads '/' as a path separator; renaming needs consent."""
        biocframe = pytest.importorskip("biocframe")

        cfg = _make_cfg(tmp_path)
        out = tmp_path / "counts.h5ad"
        args = self._make_args(str(out))
        experiment = _make_experiment(
            counts=np.array([[1, 2], [3, 4]], dtype=np.int64),
            feature_ids=["g1", "g2"],
            samples=["s1", "s2"],
            biocframe=biocframe,
            extra_column_data={
                "recount_qc__star.number_of_splices:_gt/ag": ["1", "2"]
            },
        )

        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.return_value = experiment
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            caplog.at_level(logging.ERROR),
        ):
            code = _cmd_bundle_rse(args, cfg)

        assert code == 2
        assert not out.exists()
        assert "--sanitize-columns" in caplog.text
        assert "forward slash" in caplog.text

    @pytest.mark.requires_anndata
    def test_h5ad_slash_in_column_name_round_trips_with_flag(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        anndata = pytest.importorskip("anndata")
        biocframe = pytest.importorskip("biocframe")

        cfg = _make_cfg(tmp_path)
        out = tmp_path / "counts.h5ad"
        args = self._make_args(str(out), sanitize_columns=True)
        experiment = _make_experiment(
            counts=np.array([[1, 2], [3, 4]], dtype=np.int64),
            feature_ids=["g1", "g2"],
            samples=["s1", "s2"],
            biocframe=biocframe,
            extra_column_data={
                "recount_qc__star.number_of_splices:_gt/ag": ["1", "2"]
            },
        )

        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.return_value = experiment
        res = _make_annotation_resource(cfg)
        with (
            mock.patch("recount3.cli._iter_manifest", return_value=[res]),
            mock.patch(
                "recount3.cli.R3ResourceBundle", return_value=mock_bundle
            ),
            caplog.at_level(logging.WARNING),
        ):
            code = _cmd_bundle_rse(args, cfg)

        assert code == 0
        # Every rename is reported, not applied silently.
        assert "_gt/ag -> " in caplog.text

        back = anndata.read_h5ad(out)
        assert "recount_qc__star.number_of_splices:_gt_ag" in back.obs.columns
        assert "recount_qc__star.number_of_splices:_gt/ag" not in back.obs


class TestCmdSmokeTest:
    def test_smoke_test_downloads_resources(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        args = argparse.Namespace(limit=1, dest="./recount3-smoke")
        desc = R3ResourceDescription(
            resource_type="data_source_metadata",
            organism="human",
            data_source="sra",
        )
        res = R3Resource(description=desc)
        with (
            mock.patch(
                "recount3.cli.r3_search.search_data_source_metadata",
                return_value=[res],
            ),
            mock.patch("recount3.cli._download_one") as m,
            mock.patch("pathlib.Path.mkdir"),
        ):
            m.return_value = {"url": "http://x", "status": "ok", "dest": "/f"}
            code = _cmd_smoke_test(args, cfg)
        assert code == 0
        captured = capsys.readouterr()
        evt = json.loads(captured.out.strip())
        assert evt["status"] == "ok"

    def test_smoke_test_limit_respected(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        args = argparse.Namespace(limit=1, dest="./recount3-smoke")
        desc = R3ResourceDescription(
            resource_type="data_source_metadata",
            organism="human",
            data_source="sra",
        )
        res1 = R3Resource(description=desc)
        res2 = R3Resource(description=desc)
        with (
            mock.patch(
                "recount3.cli.r3_search.search_data_source_metadata",
                return_value=[res1, res2],
            ),
            mock.patch("recount3.cli._download_one") as m,
            mock.patch("pathlib.Path.mkdir"),
        ):
            m.return_value = {"url": "http://x", "status": "ok", "dest": "/f"}
            _cmd_smoke_test(args, cfg)
        assert m.call_count == 1


class TestDispatch:
    def _make_cfg(self, tmp_path: Path) -> Config:
        return _make_cfg(tmp_path)

    def test_dispatch_ids(self, tmp_path: Path) -> None:
        cfg = self._make_cfg(tmp_path)
        args = _make_namespace(
            command="ids", organism="", samples_out=None, projects_out=None
        )
        with mock.patch("recount3.cli._cmd_ids", return_value=0) as m:
            code = _dispatch(args, cfg)
        assert code == 0
        m.assert_called_once_with(args, cfg)

    def test_dispatch_search(self, tmp_path: Path) -> None:
        cfg = self._make_cfg(tmp_path)
        args = _make_namespace(command="search")
        with mock.patch("recount3.cli._cmd_search", return_value=0) as m:
            code = _dispatch(args, cfg)
        assert code == 0
        m.assert_called_once_with(args, cfg)

    def test_dispatch_download(self, tmp_path: Path) -> None:
        cfg = self._make_cfg(tmp_path)
        args = _make_namespace(command="download")
        with mock.patch("recount3.cli._cmd_download", return_value=0) as m:
            code = _dispatch(args, cfg)
        assert code == 0

    def test_dispatch_bundle_stack_counts(self, tmp_path: Path) -> None:
        cfg = self._make_cfg(tmp_path)
        args = _make_namespace(command="bundle", bundle_cmd="stack-counts")
        with mock.patch(
            "recount3.cli._cmd_bundle_stack_counts", return_value=0
        ) as m:
            code = _dispatch(args, cfg)
        assert code == 0

    def test_dispatch_bundle_se(self, tmp_path: Path) -> None:
        cfg = self._make_cfg(tmp_path)
        args = _make_namespace(command="bundle", bundle_cmd="se")
        with mock.patch("recount3.cli._cmd_bundle_se", return_value=0) as m:
            code = _dispatch(args, cfg)
        assert code == 0

    def test_dispatch_bundle_rse(self, tmp_path: Path) -> None:
        cfg = self._make_cfg(tmp_path)
        args = _make_namespace(command="bundle", bundle_cmd="rse")
        with mock.patch("recount3.cli._cmd_bundle_rse", return_value=0) as m:
            code = _dispatch(args, cfg)
        assert code == 0

    def test_dispatch_bundle_unknown_raises_value_error(
        self, tmp_path: Path
    ) -> None:
        cfg = self._make_cfg(tmp_path)
        args = _make_namespace(command="bundle", bundle_cmd="unknown_cmd")
        with pytest.raises(ValueError, match="Unknown bundle subcommand"):
            _dispatch(args, cfg)

    def test_dispatch_smoke_test(self, tmp_path: Path) -> None:
        cfg = self._make_cfg(tmp_path)
        args = _make_namespace(command="smoke-test")
        with mock.patch("recount3.cli._cmd_smoke_test", return_value=0) as m:
            code = _dispatch(args, cfg)
        assert code == 0

    def test_dispatch_unknown_command_raises_value_error(
        self, tmp_path: Path
    ) -> None:
        cfg = self._make_cfg(tmp_path)
        args = _make_namespace(command="nonexistent")
        with pytest.raises(ValueError, match="Unknown command"):
            _dispatch(args, cfg)


class TestMain:
    def _smoke_argv(self) -> list[str]:
        """Minimal argv that runs smoke-test via main()."""
        return ["smoke-test"]

    def test_main_exits_0_on_success(self, tmp_path: Path) -> None:
        with (
            mock.patch("recount3.cli._dispatch", return_value=0),
            mock.patch("recount3.cli._build_config_from_env_and_flags"),
            pytest.raises(SystemExit) as exc,
        ):
            main(["smoke-test"])
        assert exc.value.code == 0

    def test_main_exits_130_on_keyboard_interrupt(self) -> None:
        with (
            mock.patch("recount3.cli._dispatch", side_effect=KeyboardInterrupt),
            mock.patch("recount3.cli._build_config_from_env_and_flags"),
            pytest.raises(SystemExit) as exc,
        ):
            main(["smoke-test"])
        assert exc.value.code == 130

    def test_main_exits_2_on_configuration_error(self) -> None:
        with (
            mock.patch(
                "recount3.cli._dispatch",
                side_effect=ConfigurationError("bad config"),
            ),
            mock.patch("recount3.cli._build_config_from_env_and_flags"),
            pytest.raises(SystemExit) as exc,
        ):
            main(["smoke-test"])
        assert exc.value.code == 2

    def test_main_exits_2_on_recount3_error(self) -> None:
        with (
            mock.patch(
                "recount3.cli._dispatch",
                side_effect=Recount3Error("generic error"),
            ),
            mock.patch("recount3.cli._build_config_from_env_and_flags"),
            pytest.raises(SystemExit) as exc,
        ):
            main(["smoke-test"])
        assert exc.value.code == 2

    def test_main_exits_2_on_value_error(self) -> None:
        with (
            mock.patch(
                "recount3.cli._dispatch",
                side_effect=ValueError("something wrong"),
            ),
            mock.patch("recount3.cli._build_config_from_env_and_flags"),
            pytest.raises(SystemExit) as exc,
        ):
            main(["smoke-test"])
        assert exc.value.code == 2

    def test_main_with_non_zero_code(self, tmp_path: Path) -> None:
        with (
            mock.patch("recount3.cli._dispatch", return_value=3),
            mock.patch("recount3.cli._build_config_from_env_and_flags"),
            pytest.raises(SystemExit) as exc,
        ):
            main(["smoke-test"])
        assert exc.value.code == 3

    def test_main_full_roundtrip_ids(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """End-to-end test using main() with real arg parsing."""
        with (
            mock.patch(
                "recount3.search.create_sample_project_lists",
                return_value=(["SRR001"], ["SRP001"]),
            ),
            pytest.raises(SystemExit) as exc,
        ):
            main(["ids"])
        assert exc.value.code == 0
        captured = capsys.readouterr()
        assert "SRR001" in captured.out
        assert "SRP001" in captured.out


class TestSmokeTestCustomDest:
    def test_uses_custom_dest(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        custom = str(tmp_path / "my-smoke")
        args = argparse.Namespace(limit=1, dest=custom)
        desc = R3ResourceDescription(
            resource_type="data_source_metadata",
            organism="human",
            data_source="sra",
        )
        res = R3Resource(description=desc)
        with (
            mock.patch(
                "recount3.cli.r3_search.search_data_source_metadata",
                return_value=[res],
            ),
            mock.patch("recount3.cli._download_one") as m,
        ):
            m.return_value = {
                "url": "http://x",
                "status": "ok",
                "dest": "/f",
            }
            code = _cmd_smoke_test(args, cfg)
        assert code == 0
        assert Path(custom).exists()


class TestMainBrokenPipe:
    def test_broken_pipe_exits_cleanly(self) -> None:
        with (
            mock.patch(
                "recount3.cli._dispatch",
                side_effect=BrokenPipeError,
            ),
            mock.patch("recount3.cli._build_config_from_env_and_flags"),
            mock.patch("os.open", return_value=99),
            mock.patch("os.dup2"),
            pytest.raises(SystemExit) as exc,
        ):
            main(["smoke-test"])
        assert exc.value.code == 0


class TestMainGuard:
    def test_dunder_main_calls_main(self) -> None:
        import subprocess
        import sys

        result = subprocess.run(
            [sys.executable, "-m", "recount3.cli", "--help"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        assert result.returncode == 0
        assert "recount3" in result.stdout
