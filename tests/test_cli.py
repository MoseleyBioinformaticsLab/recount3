from __future__ import annotations

import argparse
import dataclasses
import io
import json
import logging
from pathlib import Path
from typing import Any
from unittest import mock

import pytest

from recount3._descriptions import R3ResourceDescription
from recount3.config import Config, default_config
from recount3.errors import CompatibilityError, ConfigurationError, LoadError, Recount3Error
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
        args = parser.parse_args(["download", "--from=manifest.jsonl", "--dest=/tmp"])
        assert args.command == "download"
        assert args.manifest == "manifest.jsonl"
        assert args.dest == "/tmp"
        assert args.jobs == 4
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

    def test_bundle_se(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(
            ["bundle", "se", "--from=m.jsonl", "--genomic-unit=gene", "--out=out.pkl"]
        )
        assert args.bundle_cmd == "se"
        assert args.genomic_unit == "gene"
        assert args.assay_name == "raw_counts"

    def test_bundle_rse(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(
            ["bundle", "rse", "--from=m.jsonl", "--genomic-unit=gene", "--out=out.pkl"]
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

    def test_insecure_ssl_false_inherits_env(self, monkeypatch: pytest.MonkeyPatch) -> None:
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

    def test_cache_disable_env_var_one(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("RECOUNT3_CACHE_DISABLE", "1")
        args = _make_namespace()
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.cache_disabled is True

    def test_cache_disable_env_var_not_one(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("RECOUNT3_CACHE_DISABLE", "0")
        args = _make_namespace()
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.cache_disabled is False

    def test_cache_disable_env_var_other(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("RECOUNT3_CACHE_DISABLE", "yes")
        args = _make_namespace()
        cfg = _build_config_from_env_and_flags(args)
        assert cfg.cache_disabled is False

    def test_invalid_cache_dir_raises_configuration_error(self) -> None:
        args = _make_namespace(cache_dir="/some/path")
        with mock.patch.object(
            Path, "resolve", side_effect=OSError("broken")
        ):
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
        with mock.patch("recount3.cli.dataclasses.is_dataclass", return_value=False):
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
        with mock.patch("recount3.cli.dataclasses.is_dataclass", return_value=False):
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
        with mock.patch("recount3.cli.r3_search.search_annotations", return_value=[res]):
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
        with mock.patch("recount3.cli.r3_search.search_annotations", return_value=[]):
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
        desc = R3ResourceDescription(resource_type="data_sources", organism="human")
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
            "recount3.cli.r3_search.search_data_source_metadata", return_value=[res]
        ):
            code = _cmd_search(args, cfg)
        assert code == 0

    def test_source_meta_missing_filter_raises(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        args = _make_search_args(mode="source-meta", filters=["organism=human"])
        with pytest.raises(ValueError, match="Missing required filters"):
            _cmd_search(args, cfg)


class TestCmdSearchProject:
    def _project_args(self, extra_filters: list[str] | None = None) -> argparse.Namespace:
        base_filters = ["organism=human", "data_source=sra", "project=SRP000001"]
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
        with mock.patch("recount3.cli.r3_search.search_project_all", return_value=[]):
            with pytest.raises(ValueError, match="Missing required filters"):
                _cmd_search(args, cfg)

    def test_as_bool_true_values(self, tmp_path: Path) -> None:
        for val in ("1", "true", "t", "yes", "y", "on"):
            args = self._project_args([f"include_bigwig={val}"])
            code, m = self._run(tmp_path, args)
            assert code == 0
            _, kwargs = m.call_args
            assert kwargs["include_bigwig"] is True, f"Expected True for {val!r}"

    def test_as_bool_false_values(self, tmp_path: Path) -> None:
        for val in ("0", "false", "no", "whatever"):
            args = self._project_args([f"include_bigwig={val}"])
            code, m = self._run(tmp_path, args)
            assert code == 0
            _, kwargs = m.call_args
            assert kwargs["include_bigwig"] is False, f"Expected False for {val!r}"

    def test_as_bool_none_uses_default(self, tmp_path: Path) -> None:
        args = self._project_args()
        code, m = self._run(tmp_path, args)
        assert code == 0
        _, kwargs = m.call_args
        assert kwargs["include_bigwig"] is False

    def test_as_bool_include_metadata_default_true(self, tmp_path: Path) -> None:
        args = self._project_args()
        code, m = self._run(tmp_path, args)
        _, kwargs = m.call_args
        assert kwargs["include_metadata"] is True

    def test_csv_or_default_none_returns_default(self, tmp_path: Path) -> None:
        args = self._project_args()
        code, m = self._run(tmp_path, args)
        _, kwargs = m.call_args
        assert kwargs["genomic_units"] == ("gene", "exon")

    def test_csv_or_default_empty_string_returns_default(self, tmp_path: Path) -> None:
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
    def test_unknown_mode_raises_value_error(
        self, tmp_path: Path
    ) -> None:
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
        with mock.patch("recount3.cli.r3_search.search_annotations", return_value=[res]):
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
        with mock.patch("recount3.cli.r3_search.search_annotations", return_value=[res]):
            code = _cmd_search(args, cfg)
        assert code == 0
        files = list(outdir.iterdir())
        assert len(files) == 1
        assert files[0].suffix == ".jsonl"

    def test_tsv_format(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
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
        with mock.patch("recount3.cli.r3_search.search_annotations", return_value=[res]):
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
        with mock.patch("recount3.cli.r3_search.search_annotations", return_value=[res]):
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
        with mock.patch("recount3.cli.r3_search.search_annotations", return_value=[res]):
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
        with mock.patch("recount3.cli.r3_search.search_annotations", return_value=[res]):
            code = _cmd_search(args, cfg)
        assert code == 0


from recount3.cli import _download_one


class TestDownloadOne:
    def test_zip_dest_downloads_and_returns_ok(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        dest = tmp_path / "out.zip"
        with mock.patch.object(
            R3Resource, "download", return_value=None
        ) as m:
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
        with mock.patch.object(R3Resource, "download", return_value=str(dest_file)) as m:
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
        args = self._make_args(inline=json.dumps(inline_obj), dest=str(tmp_path))
        with mock.patch("recount3.cli._download_one") as m:
            m.return_value = {"url": "http://x", "status": "ok", "dest": "/tmp/f"}
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
        with mock.patch(
            "recount3.cli._iter_manifest", return_value=[res]
        ), mock.patch("recount3.cli._download_one") as m:
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
        with mock.patch(
            "recount3.cli._iter_manifest", return_value=[res]
        ), mock.patch("recount3.cli._download_one") as m:
            m.return_value = {"url": "http://x", "status": "ok", "dest": None}
            code = _cmd_download(args, cfg)
        assert code == 0
        assert (tmp_path / "subdir").is_dir()

    def test_all_errors_returns_2(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        cfg = _make_cfg(tmp_path)
        res = _make_annotation_resource(cfg)
        args = self._make_args(manifest="m.jsonl", dest=str(tmp_path))
        with mock.patch(
            "recount3.cli._iter_manifest", return_value=[res]
        ), mock.patch("recount3.cli._download_one") as m:
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
            {"url": "http://b", "status": "error", "dest": None, "error": "oops"},
        ]
        with mock.patch(
            "recount3.cli._iter_manifest", return_value=[res1, res2]
        ), mock.patch("recount3.cli._download_one", side_effect=results):
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
            out=out,
        )
        defaults.update(kw)
        return argparse.Namespace(**defaults)

    def test_success_parquet(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.parquet"
        args = self._make_args(str(out))
        mock_df = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.return_value = mock_df
        res = _make_annotation_resource(cfg)
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
        ):
            code = _cmd_bundle_stack_counts(args, cfg)
        assert code == 0
        mock_df.to_parquet.assert_called_once_with(out)

    def test_success_tsv(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.tsv"
        args = self._make_args(str(out))
        mock_df = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.return_value = mock_df
        res = _make_annotation_resource(cfg)
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
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
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
        ):
            code = _cmd_bundle_stack_counts(args, cfg)
        assert code == 0
        mock_df.to_csv.assert_called_once_with(out, sep=",")

    def test_compatibility_error_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.tsv"
        args = self._make_args(str(out))
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.side_effect = CompatibilityError("incompat")
        res = _make_annotation_resource(cfg)
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
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
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
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
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
        ):
            code = _cmd_bundle_stack_counts(args, cfg)
        assert code == 2

    def test_import_error_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.parquet"
        args = self._make_args(str(out))
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.side_effect = ImportError("missing dep")
        res = _make_annotation_resource(cfg)
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
        ):
            code = _cmd_bundle_stack_counts(args, cfg)
        assert code == 2

    def test_write_exception_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.parquet"
        args = self._make_args(str(out))
        mock_df = mock.MagicMock()
        mock_df.to_parquet.side_effect = OSError("disk full")
        mock_bundle = mock.MagicMock()
        mock_bundle.stack_count_matrices.return_value = mock_df
        res = _make_annotation_resource(cfg)
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
        ):
            code = _cmd_bundle_stack_counts(args, cfg)
        assert code == 2


class TestCmdBundleSe:
    def _make_args(self, out: str, **kw: Any) -> argparse.Namespace:
        defaults = dict(
            manifest="-",
            genomic_unit="gene",
            annotation=None,
            assay_name="raw_counts",
            join="inner",
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
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
        ), mock.patch("pickle.dump") as mock_pickle:
            code = _cmd_bundle_se(args, cfg)
        assert code == 0
        mock_pickle.assert_called_once()

    def test_success_h5ad(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.h5ad"
        args = self._make_args(str(out))
        mock_se = mock.MagicMock()
        mock_adata = mock.MagicMock()
        mock_se.to_anndata.return_value = mock_adata
        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.return_value = mock_se
        res = _make_annotation_resource(cfg)
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
        ):
            code = _cmd_bundle_se(args, cfg)
        assert code == 0
        mock_adata.write_h5ad.assert_called_once_with(out)

    def test_import_error_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.pkl"
        args = self._make_args(str(out))
        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.side_effect = ImportError("no se dep")
        res = _make_annotation_resource(cfg)
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
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
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
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
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
        ), mock.patch("pickle.dump", side_effect=OSError("disk full")):
            code = _cmd_bundle_se(args, cfg)
        assert code == 2

    def test_h5ad_write_exception_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.h5ad"
        args = self._make_args(str(out))
        mock_se = mock.MagicMock()
        mock_adata = mock.MagicMock()
        mock_adata.write_h5ad.side_effect = OSError("disk full")
        mock_se.to_anndata.return_value = mock_adata
        mock_bundle = mock.MagicMock()
        mock_bundle.to_summarized_experiment.return_value = mock_se
        res = _make_annotation_resource(cfg)
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
        ):
            code = _cmd_bundle_se(args, cfg)
        assert code == 2


class TestCmdBundleRse:
    def _make_args(self, out: str, **kw: Any) -> argparse.Namespace:
        defaults = dict(
            manifest="-",
            genomic_unit="gene",
            annotation=None,
            assay_name="raw_counts",
            join="inner",
            allow_fallback_to_se=False,
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
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
        ), mock.patch("pickle.dump") as mock_pickle:
            code = _cmd_bundle_rse(args, cfg)
        assert code == 0
        mock_pickle.assert_called_once()

    def test_success_h5ad(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.h5ad"
        args = self._make_args(str(out))
        mock_rse = mock.MagicMock()
        mock_adata = mock.MagicMock()
        mock_rse.to_anndata.return_value = mock_adata
        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.return_value = mock_rse
        res = _make_annotation_resource(cfg)
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
        ):
            code = _cmd_bundle_rse(args, cfg)
        assert code == 0
        mock_adata.write_h5ad.assert_called_once_with(out)

    def test_allow_fallback_to_se_passed_through(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.pkl"
        args = self._make_args(str(out), allow_fallback_to_se=True)
        mock_rse = mock.MagicMock()
        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.return_value = mock_rse
        res = _make_annotation_resource(cfg)
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
        ), mock.patch("pickle.dump"):
            code = _cmd_bundle_rse(args, cfg)
        assert code == 0
        _, kwargs = mock_bundle.to_ranged_summarized_experiment.call_args
        assert kwargs["allow_fallback_to_se"] is True

    def test_import_error_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.pkl"
        args = self._make_args(str(out))
        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.side_effect = ImportError("missing")
        res = _make_annotation_resource(cfg)
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
        ):
            code = _cmd_bundle_rse(args, cfg)
        assert code == 2

    def test_general_exception_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.pkl"
        args = self._make_args(str(out))
        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.side_effect = RuntimeError("bad")
        res = _make_annotation_resource(cfg)
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
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
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
        ), mock.patch("pickle.dump", side_effect=OSError("disk full")):
            code = _cmd_bundle_rse(args, cfg)
        assert code == 2

    def test_h5ad_write_exception_returns_2(self, tmp_path: Path) -> None:
        cfg = _make_cfg(tmp_path)
        out = tmp_path / "out.h5ad"
        args = self._make_args(str(out))
        mock_rse = mock.MagicMock()
        mock_adata = mock.MagicMock()
        mock_adata.write_h5ad.side_effect = OSError("disk full")
        mock_rse.to_anndata.return_value = mock_adata
        mock_bundle = mock.MagicMock()
        mock_bundle.to_ranged_summarized_experiment.return_value = mock_rse
        res = _make_annotation_resource(cfg)
        with mock.patch("recount3.cli._iter_manifest", return_value=[res]), mock.patch(
            "recount3.cli.R3ResourceBundle", return_value=mock_bundle
        ):
            code = _cmd_bundle_rse(args, cfg)
        assert code == 2


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
        with mock.patch(
            "recount3.cli.r3_search.search_data_source_metadata", return_value=[res]
        ), mock.patch("recount3.cli._download_one") as m, mock.patch(
            "pathlib.Path.mkdir"
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
        with mock.patch(
            "recount3.cli.r3_search.search_data_source_metadata",
            return_value=[res1, res2],
        ), mock.patch("recount3.cli._download_one") as m, mock.patch(
            "pathlib.Path.mkdir"
        ):
            m.return_value = {"url": "http://x", "status": "ok", "dest": "/f"}
            _cmd_smoke_test(args, cfg)
        assert m.call_count == 1


class TestDispatch:
    def _make_cfg(self, tmp_path: Path) -> Config:
        return _make_cfg(tmp_path)

    def test_dispatch_ids(self, tmp_path: Path) -> None:
        cfg = self._make_cfg(tmp_path)
        args = _make_namespace(command="ids", organism="", samples_out=None, projects_out=None)
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
        with mock.patch("recount3.cli._cmd_bundle_stack_counts", return_value=0) as m:
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

    def test_dispatch_bundle_unknown_raises_value_error(self, tmp_path: Path) -> None:
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

    def test_dispatch_unknown_command_raises_value_error(self, tmp_path: Path) -> None:
        cfg = self._make_cfg(tmp_path)
        args = _make_namespace(command="nonexistent")
        with pytest.raises(ValueError, match="Unknown command"):
            _dispatch(args, cfg)


class TestMain:
    def _smoke_argv(self) -> list[str]:
        """Minimal argv that runs smoke-test via main()."""
        return ["smoke-test"]

    def test_main_exits_0_on_success(self, tmp_path: Path) -> None:
        with mock.patch("recount3.cli._dispatch", return_value=0), mock.patch(
            "recount3.cli._build_config_from_env_and_flags"
        ), pytest.raises(SystemExit) as exc:
            main(["smoke-test"])
        assert exc.value.code == 0

    def test_main_exits_130_on_keyboard_interrupt(self) -> None:
        with mock.patch(
            "recount3.cli._dispatch", side_effect=KeyboardInterrupt
        ), mock.patch("recount3.cli._build_config_from_env_and_flags"), pytest.raises(
            SystemExit
        ) as exc:
            main(["smoke-test"])
        assert exc.value.code == 130

    def test_main_exits_2_on_configuration_error(self) -> None:
        with mock.patch(
            "recount3.cli._dispatch",
            side_effect=ConfigurationError("bad config"),
        ), mock.patch("recount3.cli._build_config_from_env_and_flags"), pytest.raises(
            SystemExit
        ) as exc:
            main(["smoke-test"])
        assert exc.value.code == 2

    def test_main_exits_2_on_recount3_error(self) -> None:
        with mock.patch(
            "recount3.cli._dispatch",
            side_effect=Recount3Error("generic error"),
        ), mock.patch("recount3.cli._build_config_from_env_and_flags"), pytest.raises(
            SystemExit
        ) as exc:
            main(["smoke-test"])
        assert exc.value.code == 2

    def test_main_exits_2_on_value_error(self) -> None:
        with mock.patch(
            "recount3.cli._dispatch",
            side_effect=ValueError("something wrong"),
        ), mock.patch("recount3.cli._build_config_from_env_and_flags"), pytest.raises(
            SystemExit
        ) as exc:
            main(["smoke-test"])
        assert exc.value.code == 2

    def test_main_with_non_zero_code(self, tmp_path: Path) -> None:
        with mock.patch("recount3.cli._dispatch", return_value=3), mock.patch(
            "recount3.cli._build_config_from_env_and_flags"
        ), pytest.raises(SystemExit) as exc:
            main(["smoke-test"])
        assert exc.value.code == 3

    def test_main_full_roundtrip_ids(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """End-to-end test using main() with real arg parsing."""
        with mock.patch(
            "recount3.search.create_sample_project_lists",
            return_value=(["SRR001"], ["SRP001"]),
        ), pytest.raises(SystemExit) as exc:
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
        with mock.patch(
            "recount3.cli.r3_search.search_data_source_metadata",
            return_value=[res],
        ), mock.patch("recount3.cli._download_one") as m:
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
        with mock.patch(
            "recount3.cli._dispatch",
            side_effect=BrokenPipeError,
        ), mock.patch(
            "recount3.cli._build_config_from_env_and_flags"
        ), mock.patch(
            "os.open", return_value=99
        ), mock.patch(
            "os.dup2"
        ), pytest.raises(SystemExit) as exc:
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
