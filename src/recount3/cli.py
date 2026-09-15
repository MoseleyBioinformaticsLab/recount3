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
#
"""Command-line interface for recount3.

This module is the adapter between command-line arguments and the package's
Python API. It builds the ``argparse`` command tree, resolves configuration,
encodes and reads manifests, coordinates downloads, writes bundle outputs,
and translates expected failures into process exit codes.

The installed ``recount3`` command and ``python -m recount3`` both call
:func:`main`. Importing this module does not parse arguments or perform I/O.
The end-user command reference and shell recipes live in ``docs/cli.rst``.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import dataclasses
from datetime import datetime
import json
import logging
import os
import pickle
from pathlib import Path
import sys
import types
from typing import Any, Iterable, Iterator, Mapping

from recount3.bundle import R3ResourceBundle
from recount3.config import Config, default_config
from recount3.errors import (
    CompatibilityError,
    ConfigurationError,
    LoadError,
    Recount3Error,
)
from recount3.resource import R3Resource
from recount3.types import CacheMode, CompatibilityMode
from recount3 import _utils
from recount3 import search as r3_search
from recount3.version import __version__


def _build_parser() -> argparse.ArgumentParser:
    """Return the top-level argument parser for the recount3 CLI.

    The parser is split into subparsers for each top-level job. Global flags
    that map to the run-time :class:`Config` are defined here and then read by
    each subcommand through :func:`_build_config_from_env_and_flags`.

    Returns:
      A configured :class:`argparse.ArgumentParser` instance.
    """
    parser = argparse.ArgumentParser(
        prog="recount3",
        description=(
            "Discover, download, and assemble recount3 RNA-seq "
            "resources.\n\n"
            "Typical workflow:\n"
            "  1. recount3 search  — discover resources, "
            "produce a JSONL manifest\n"
            "  2. recount3 download — download resources "
            "from a manifest\n"
            "  3. recount3 bundle  — assemble downloaded "
            "resources into analysis objects\n\n"
            "All flags below are optional and override "
            "environment variables\n"
            "(RECOUNT3_URL, RECOUNT3_CACHE_DIR, "
            "RECOUNT3_HTTP_TIMEOUT, etc.)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--base-url",
        default=None,
        metavar="URL",
        help=(
            "Override the recount3 mirror base URL. "
            "Env: RECOUNT3_URL. "
            "Default: http://duffel.rail.bio/recount3/."
        ),
    )
    parser.add_argument(
        "--cache-backend",
        choices=("filesystem", "pybiocfilecache"),
        default=None,
        help="Cache registry. Env: RECOUNT3_CACHE_BACKEND.",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        metavar="DIR",
        help=(
            "Local directory for caching downloaded files. "
            "Env: RECOUNT3_CACHE_DIR. "
            "Default: ~/.cache/recount3/files for filesystem; "
            "R's recount3 cache directory for pybiocfilecache."
        ),
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=None,
        metavar="SEC",
        help=(
            "HTTP request timeout in seconds. "
            "Env: RECOUNT3_HTTP_TIMEOUT. Default: 60."
        ),
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Maximum retry attempts for transient HTTP "
            "errors (5xx, timeouts). "
            "Env: RECOUNT3_MAX_RETRIES. Default: 3."
        ),
    )
    parser.add_argument(
        "--insecure-ssl",
        action="store_true",
        help=(
            "Disable TLS certificate verification. "
            "Applies only to https base URLs (no-op for the default http "
            "Duffel mirror). NOT recommended outside debugging. "
            "Env: RECOUNT3_INSECURE_SSL=1."
        ),
    )
    parser.add_argument(
        "--user-agent",
        default=None,
        metavar="STRING",
        help=("Custom HTTP User-Agent header. " "Env: RECOUNT3_USER_AGENT."),
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=None,
        metavar="BYTES",
        help=(
            "Streaming read/write chunk size for downloads "
            "and file copies. "
            "Env: RECOUNT3_CHUNK_SIZE. Default: 1048576 (1 MiB)."
        ),
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help=(
            "Suppress informational messages; only show "
            "warnings and errors (WARNING level)."
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help=(
            "Show detailed debug output including HTTP "
            "requests and cache operations (DEBUG level)."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="<command>",
        required=True,
    )

    p_ids = subparsers.add_parser(
        "ids",
        help="Emit unique sample and project ID lists.",
        description=(
            "Query the recount3 mirror for every known sample "
            "and project ID. Useful for scripting loops or "
            "verifying that a project exists before running a "
            "full search."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p_ids.add_argument(
        "--organism",
        default="",
        help=(
            'Restrict to "human" or "mouse". '
            "Omit or leave empty to list IDs for both."
        ),
    )
    p_ids.add_argument(
        "--samples-out",
        default=None,
        metavar="FILE",
        help=(
            "Write the sample ID list to FILE "
            "(one ID per line). Prints to stdout if omitted."
        ),
    )
    p_ids.add_argument(
        "--projects-out",
        default=None,
        metavar="FILE",
        help=(
            "Write the project ID list to FILE "
            "(one ID per line). Prints to stdout if omitted."
        ),
    )

    p_search = subparsers.add_parser(
        "search",
        help="Discover resources and emit a manifest.",
        description=(
            "Discover recount3 resources and emit a "
            "JSONL or TSV manifest.\n\n"
            "Choose a MODE to select the resource family, "
            "then supply the required key=value filters as "
            "positional arguments. The manifest can be piped "
            "directly into 'recount3 download --from=-'."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Required filters per mode:\n"
            "  annotations: organism, genomic_unit,"
            " annotation_extension\n"
            "  gene-exon:   organism, data_source,"
            " genomic_unit, project\n"
            "               (optional:"
            " annotation_extension; default G026)\n"
            "  junctions:   organism, data_source, project\n"
            "               (optional: junction_type=ALL,"
            " junction_extension=MM)\n"
            "  metadata:    organism, data_source,"
            " table_name, project\n"
            "  bigwig:      organism, data_source,"
            " project, sample\n"
            "  project:     organism, data_source, project\n"
            "               (optional: genomic_unit=gene,exon;"
            " annotation=default|all|name;\n"
            "               annotation_extension=G026,G029"
            " [overrides annotation];\n"
            "               junction_extension=MM,RR,ID;"
            " include_metadata=true|false;\n"
            "               include_bigwig=true|false)\n"
            "  sources:     organism\n"
            "  source-meta: organism, data_source\n"
            "\nFilter value reference:\n"
            "  organism:             human | mouse\n"
            "  data_source:          sra | gtex | tcga\n"
            "  genomic_unit:         gene | exon | junction\n"
            "  annotation_extension: G026 | G029 | F006 | "
            "R109 | ERCC | SIRV\n"
            "  table_name:           recount_qc | "
            "recount_project_info | ...\n"
            "\nExamples:\n"
            "  # Discover gene annotations and save as "
            "JSONL:\n"
            "  recount3 search annotations "
            "organism=human genomic_unit=gene \\\n"
            "      annotation_extension=G026 "
            "--format=jsonl > manifest.jsonl\n\n"
            "  # Search gene counts and pipe directly "
            "into download:\n"
            "  recount3 search gene-exon organism=human "
            "data_source=sra \\\n"
            "      genomic_unit=gene project=SRP009615 "
            "--format=jsonl | \\\n"
            "  recount3 download --from=- --dest=./out\n\n"
            "  # Discover all resources for a project:\n"
            "  recount3 search project organism=human "
            "data_source=sra \\\n"
            "      project=SRP009615\n"
        ),
    )
    p_search.add_argument(
        "mode",
        choices=(
            "annotations",
            "gene-exon",
            "junctions",
            "metadata",
            "bigwig",
            "project",
            "sources",
            "source-meta",
        ),
        metavar="MODE",
        help=(
            "Resource family to search. Choices: "
            "annotations, gene-exon, junctions, metadata, "
            "bigwig, project, sources, source-meta."
        ),
    )
    p_search.add_argument(
        "filters",
        nargs="*",
        metavar="key=value",
        help=(
            "One or more key=value filter pairs "
            "(e.g., organism=human data_source=sra "
            "project=SRP009615). Required filters depend "
            "on MODE; see epilog below."
        ),
    )
    p_search.add_argument(
        "--format",
        choices=("jsonl", "tsv"),
        default="jsonl",
        help=(
            "Output format. 'jsonl' emits one JSON object "
            "per resource (ideal for piping into download). "
            "'tsv' emits a tab-separated table."
        ),
    )
    out_dest = p_search.add_mutually_exclusive_group()
    out_dest.add_argument(
        "--output",
        default=None,
        metavar="FILE",
        help=(
            "Write the manifest to FILE instead of stdout. "
            "Mutually exclusive with --outdir."
        ),
    )
    out_dest.add_argument(
        "--outdir",
        default=None,
        metavar="DIR",
        help=(
            "Write the manifest to DIR with an "
            "auto-generated timestamped filename "
            "(e.g., gene-exon-20250315-091200.jsonl). "
            "Mutually exclusive with --output."
        ),
    )

    p_dl = subparsers.add_parser(
        "download",
        help="Download resources from a manifest.",
        description=(
            "Download recount3 resources described by a "
            "JSONL manifest (from 'recount3 search') or a "
            "single inline JSON object.\n\n"
            "Progress is streamed as one JSON event per "
            "resource to stdout (fields: url, status, dest, "
            "error). Status is 'ok', 'skipped', or 'error'."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Download from a saved manifest file:\n"
            "  recount3 download --from=manifest.jsonl "
            "--dest=./data --jobs=8\n\n"
            "  # Pipe search output directly into "
            "download:\n"
            "  recount3 search gene-exon organism=human "
            "data_source=sra \\\n"
            "      genomic_unit=gene project=SRP009615 "
            "--format=jsonl | \\\n"
            "  recount3 download --from=- --dest=./data\n\n"
            "  # Download into a zip archive:\n"
            "  recount3 download --from=manifest.jsonl "
            "--dest=bundle.zip\n\n"
            "Exit codes:\n"
            "  0  All resources downloaded successfully.\n"
            "  3  Some resources failed (partial failure)."
        ),
    )
    src = p_dl.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--from",
        dest="manifest",
        default=None,
        metavar="FILE",
        help=(
            "Path to a JSONL manifest file produced by "
            "'recount3 search'. Use '-' to read from "
            "stdin (enables piping)."
        ),
    )
    src.add_argument(
        "--inline",
        default=None,
        metavar="JSON",
        help=(
            "A single resource as an inline JSON object "
            "on the command line. Useful for downloading "
            "one resource without creating a manifest file."
        ),
    )
    p_dl.add_argument(
        "--dest",
        default=".",
        metavar="PATH",
        help=(
            "Destination directory or .zip file path. "
            "If a directory, files are placed flat inside "
            "it. If a .zip path, all resources are bundled "
            "into the archive."
        ),
    )
    p_dl.add_argument(
        "--overwrite",
        action="store_true",
        help=(
            "Overwrite files that already exist at the "
            "destination. Without this flag, existing files "
            "are skipped (status='skipped')."
        ),
    )
    p_dl.add_argument(
        "--jobs",
        type=int,
        default=8,
        metavar="N",
        help=(
            "Number of parallel download threads. "
            "Increase for high-bandwidth connections or "
            "large manifests."
        ),
    )
    p_dl.add_argument(
        "--cache",
        choices=("enable", "disable", "update"),
        default="enable",
        help=(
            "Cache behavior. 'enable': use local cache, "
            "download only if missing. 'disable': stream "
            "directly to dest, bypass cache entirely. "
            "'update': re-download into cache first, then "
            "copy to dest."
        ),
    )

    p_bundle = subparsers.add_parser(
        "bundle",
        help="Assemble resources into analysis objects.",
        description=(
            "Assemble multiple recount3 resources into "
            "analysis-ready objects: stacked count "
            "matrices, SummarizedExperiments, or "
            "RangedSummarizedExperiments.\n\n"
            "Each bundle subcommand reads a JSONL manifest "
            "(from 'recount3 search') and produces a single "
            "output file."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sp_bundle = p_bundle.add_subparsers(
        dest="bundle_cmd",
        metavar="<bundle-cmd>",
        required=True,
    )

    p_stack = sp_bundle.add_parser(
        "stack-counts",
        help="Stack count matrices into a single table.",
        description=(
            "Load gene/exon or junction count matrices "
            "from a manifest and concatenate them into a "
            "single pandas DataFrame. Outputs TSV, "
            "gzip-compressed TSV, or Parquet."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p_stack.add_argument(
        "--from",
        dest="manifest",
        required=True,
        metavar="FILE",
        help=(
            "Path to a JSONL manifest (or '-' for stdin) "
            "containing count-file resources."
        ),
    )
    p_stack.add_argument(
        "--compat",
        choices=("family", "feature"),
        default="family",
        help=(
            "Compatibility mode for stacking. 'family': "
            "require all resources to be the same type "
            "(gene, exon, or junction). 'feature': allow "
            "mixing types that share the same feature index."
        ),
    )
    p_stack.add_argument(
        "--join",
        choices=("inner", "outer"),
        default="inner",
        help=(
            "How to handle mismatched row/column indices "
            "across matrices. 'inner': keep only shared "
            "features/samples. 'outer': keep all, fill "
            "missing values with NaN."
        ),
    )
    p_stack.add_argument(
        "--axis",
        type=int,
        choices=(0, 1),
        default=1,
        help=(
            "Concatenation axis. 0: stack rows (add more "
            "features). 1: stack columns (add more samples, "
            "the typical case)."
        ),
    )
    p_stack.add_argument(
        "--verify-integrity",
        action="store_true",
        help=(
            "Raise an error if the resulting index "
            "contains duplicates. Useful for catching "
            "accidentally overlapping manifests."
        ),
    )
    p_stack.add_argument(
        "--densify",
        action="store_true",
        help=(
            "Convert sparse columns to dense before "
            "writing Parquet. Junction count matrices "
            "load sparse and no Parquet engine accepts "
            "pandas sparse dtypes, so .parquet output "
            "fails without this flag. Densifying "
            "materializes every implicit zero and can "
            "need far more memory. Ignored for text "
            "formats, which already write dense rows."
        ),
    )
    p_stack.add_argument(
        "--out",
        required=True,
        metavar="FILE",
        help=(
            "Output file path. Extension determines "
            "format: .csv, .tsv (tab-separated), .tsv.gz "
            "(gzip-compressed TSV), or .parquet "
            "(Apache Parquet; requires a Parquet engine: "
            "pip install 'recount3[parquet]', or any "
            "pyarrow/fastparquet install pandas accepts). "
            "Engine availability is checked before "
            "downloading anything."
        ),
    )

    p_se = sp_bundle.add_parser(
        "se",
        help="Build a SummarizedExperiment.",
        description=(
            "Assemble resources from a manifest into a "
            "BiocPy SummarizedExperiment. Requires the "
            "'recount3[biocpy]' extra. Output as a Python "
            "pickle (.pkl), or as an AnnData HDF5 file "
            "(.h5ad) with the 'recount3[anndata]' extra. "
            "Availability is checked before downloading "
            "anything."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p_se.add_argument(
        "--from",
        dest="manifest",
        required=True,
        metavar="FILE",
        help=(
            "Path to a JSONL manifest (or '-' for stdin) "
            "produced by 'recount3 search'."
        ),
    )
    p_se.add_argument(
        "--genomic-unit",
        choices=("gene", "exon", "junction"),
        required=True,
        help=(
            "Which feature type to assemble. Only "
            "resources matching this unit are included."
        ),
    )
    p_se.add_argument(
        "--annotation",
        default=None,
        metavar="KEY",
        help=(
            "Annotation extension to select for gene/exon "
            "counts (e.g., 'G026', 'G029', 'R109'). "
            "If omitted, uses the first available."
        ),
    )
    p_se.add_argument(
        "--assay-name",
        default="raw_counts",
        help=(
            "Name for the assay slot in the SE. Change "
            "this if you plan to merge multiple assays."
        ),
    )
    p_se.add_argument(
        "--join",
        choices=("inner", "outer"),
        default="inner",
        help=(
            "Join policy when stacking matrices across "
            "projects. 'inner': only shared features. "
            "'outer': union of all features (NaN fill)."
        ),
    )
    p_se.add_argument(
        "--sanitize-columns",
        action="store_true",
        help=(
            "Replace '/' with '_' in obs/var column "
            "names before writing .h5ad. HDF5 reads a "
            "forward slash as a path separator, and "
            "recount3 STAR QC fields are named after "
            "splice motifs (for example, "
            "'recount_qc__star.number_of_splices:_gt/ag'), "
            "so .h5ad output fails without this flag. "
            "It renames the columns your analysis indexes "
            "by. Ignored for .pkl output, which keeps the "
            "names verbatim."
        ),
    )
    p_se.add_argument(
        "--out",
        required=True,
        metavar="FILE",
        help=(
            "Output file path. Use .pkl for a Python "
            "pickle or .h5ad for an AnnData HDF5 file "
            "(requires the 'recount3[anndata]' extra)."
        ),
    )

    p_rse = sp_bundle.add_parser(
        "rse",
        help="Build a RangedSummarizedExperiment.",
        description=(
            "Assemble resources from a manifest into a "
            "BiocPy RangedSummarizedExperiment, which adds "
            "genomic coordinates (GenomicRanges) to each "
            "feature row. Requires the 'recount3[biocpy]' "
            "extra. Output as .pkl, or as .h5ad with the "
            "'recount3[anndata]' extra."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p_rse.add_argument(
        "--from",
        dest="manifest",
        required=True,
        metavar="FILE",
        help=(
            "Path to a JSONL manifest (or '-' for stdin) "
            "produced by 'recount3 search'."
        ),
    )
    p_rse.add_argument(
        "--genomic-unit",
        choices=("gene", "exon", "junction"),
        required=True,
        help=(
            "Which feature type to assemble. Only "
            "resources matching this unit are included."
        ),
    )
    p_rse.add_argument(
        "--annotation",
        default=None,
        metavar="KEY",
        help=(
            "Annotation extension to select for gene/exon "
            "counts (e.g., 'G026', 'G029', 'R109'). "
            "If omitted, uses the first available."
        ),
    )
    p_rse.add_argument(
        "--assay-name",
        default="raw_counts",
        help=(
            "Name for the assay slot in the RSE. Change "
            "this if you plan to merge multiple assays."
        ),
    )
    p_rse.add_argument(
        "--join",
        choices=("inner", "outer"),
        default="inner",
        help=(
            "Join policy when stacking matrices across "
            "projects. 'inner': only shared features. "
            "'outer': union of all features (NaN fill)."
        ),
    )
    p_rse.add_argument(
        "--allow-fallback-to-se",
        action="store_true",
        help=(
            "If genomic ranges cannot be derived for the "
            "selected annotation, fall back to emitting a "
            "plain SummarizedExperiment instead of raising "
            "an error."
        ),
    )
    p_rse.add_argument(
        "--sanitize-columns",
        action="store_true",
        help=(
            "Replace '/' with '_' in obs/var column "
            "names before writing .h5ad. HDF5 reads a "
            "forward slash as a path separator, and "
            "recount3 STAR QC fields are named after "
            "splice motifs (for example, "
            "'recount_qc__star.number_of_splices:_gt/ag'), "
            "so .h5ad output fails without this flag. "
            "It renames the columns your analysis indexes "
            "by. Ignored for .pkl output, which keeps the "
            "names verbatim."
        ),
    )
    p_rse.add_argument(
        "--out",
        required=True,
        metavar="FILE",
        help=(
            "Output file path. Use .pkl for a Python "
            "pickle or .h5ad for an AnnData HDF5 file "
            "(requires the 'recount3[anndata]' extra)."
        ),
    )

    p_smoke = subparsers.add_parser(
        "smoke-test",
        help="Quick connectivity smoke test.",
        description=(
            "Verify that the recount3 mirror is reachable "
            "by discovering and downloading a small number "
            "of resources. Intended for CI pipelines, dev "
            "setup verification, and network diagnostics."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p_smoke.add_argument(
        "--dest",
        default="./recount3-smoke",
        metavar="DIR",
        help=(
            "Destination directory for smoke-test "
            "downloads (default: ./recount3-smoke)."
        ),
    )
    p_smoke.add_argument(
        "--limit",
        type=int,
        default=1,
        metavar="N",
        help=(
            "Number of resources to discover and attempt "
            "to download. Increase to exercise more of "
            "the mirror."
        ),
    )

    return parser


def _build_config_from_env_and_flags(args: argparse.Namespace) -> Config:
    """Construct a :class:`Config` from environment and CLI flags.

    Flag values override environment; unspecified flags keep defaults from
    :func:`default_config`.

    Args:
      args: Parsed CLI arguments.

    Returns:
      A :class:`Config` instance suitable for network and caching ops.

    Raises:
      ConfigurationError: If a provided path is invalid or incompatible.
    """
    base = default_config(cache_backend=getattr(args, "cache_backend", None))

    base_url = args.base_url or base.base_url
    cache_dir = Path(args.cache_dir) if args.cache_dir else base.cache_dir
    timeout = args.timeout if args.timeout is not None else base.timeout
    retries = args.retries if args.retries is not None else base.max_retries
    insecure_ssl = bool(args.insecure_ssl) or base.insecure_ssl
    user_agent = args.user_agent or base.user_agent
    chunk_size = (
        args.chunk_size if args.chunk_size is not None else base.chunk_size
    )

    try:
        cache_dir = cache_dir.expanduser().resolve()
    except Exception as exc:
        raise ConfigurationError(
            f"Invalid cache directory: {cache_dir!s}"
        ) from exc

    return Config(
        base_url=base_url,
        timeout=int(timeout),
        insecure_ssl=bool(insecure_ssl),
        max_retries=int(retries),
        user_agent=str(user_agent),
        cache_dir=cache_dir,
        cache_backend=base.cache_backend,
        cache_disabled=os.environ.get("RECOUNT3_CACHE_DISABLE", "0") == "1",
        chunk_size=int(chunk_size),
    )


def _init_logging(args: argparse.Namespace) -> None:
    """Configure logging verbosity per user preference.

    Args:
      args: Parsed CLI arguments with ``--quiet`` or ``--verbose`` flags.
    """
    level = logging.INFO
    if args.quiet:
        level = logging.WARNING
    if args.verbose:
        level = logging.DEBUG
    logging.basicConfig(
        level=level, format="%(levelname)s:%(name)s:%(message)s"
    )


def _parse_filters(tokens: Iterable[str]) -> dict[str, str]:
    """Parse ``key=value`` tokens into a dictionary.

    Args:
      tokens: Iterable of strings, each of the form ``k=v``.

    Returns:
      A mapping of lower-cased keys to their values (unmodified).

    Raises:
      ValueError: If a token is not of the form ``k=v``.
    """
    result: dict[str, str] = {}
    for t in tokens:
        if "=" not in t:
            raise ValueError(
                f"Expected 'key=value' filter, got: {t!r} "
                f"(ensure no spaces around '=' or after commas)"
            )
        k, v = t.split("=", 1)
        k = k.strip().lower()
        v = v.strip()
        if not k:
            raise ValueError(f"Empty filter key in token: {t!r}")
        result[k] = v
    return result


_PROJECT_SELECTORS = frozenset({"organism", "data_source", "project"})
_JUNCTION_SELECTORS = frozenset({"junction_type", "junction_extension"})

_SEARCH_MODE_FILTERS: Mapping[str, frozenset[str]] = types.MappingProxyType(
    {
        "annotations": frozenset(
            {"organism", "genomic_unit", "annotation_extension"}
        ),
        "gene-exon": _PROJECT_SELECTORS
        | {"genomic_unit", "annotation_extension"},
        "junctions": _PROJECT_SELECTORS | _JUNCTION_SELECTORS,
        "metadata": _PROJECT_SELECTORS | {"table_name"},
        "bigwig": _PROJECT_SELECTORS | {"sample"},
        "project": _PROJECT_SELECTORS
        | _JUNCTION_SELECTORS
        | {
            "genomic_unit",
            "annotation",
            "annotation_extension",
            "include_metadata",
            "include_bigwig",
        },
        "sources": frozenset({"organism"}),
        "source-meta": frozenset({"organism", "data_source"}),
    }
)

_PYTHON_FILTER_ALIASES: Mapping[str, str] = types.MappingProxyType(
    {
        "annotations": "annotation",
        "data_sources": "data_source",
        "genomic_units": "genomic_unit",
        "junction_exts": "junction_extension",
        "junction_extensions": "junction_extension",
        "junction_types": "junction_type",
        "organisms": "organism",
        "projects": "project",
        "samples": "sample",
        "table_names": "table_name",
    }
)


def _reject_unknown_filters(mode: str, filters: Mapping[str, str]) -> None:
    """Raise if a filter key is one the given search mode never reads.

    ``key=value`` selectors are free-form, so a misspelling or a Python
    keyword name used at the shell would otherwise be dropped without a word
    and the command would emit its defaults instead. For ``search project``
    that means a full ten-resource manifest in place of the narrowed one the
    user asked for, which is wrong rather than merely unhelpful.

    Args:
      mode: The ``search`` subcommand mode, for example ``"project"``.
      filters: Parsed ``key=value`` selectors.

    Raises:
      ValueError: If any key is not read by ``mode``. The message names each
        offending key, the CLI spelling when the key is a known Python-API
        alias, and the keys the mode does accept.
    """
    allowed = _SEARCH_MODE_FILTERS.get(mode)
    if allowed is None:
        return

    unknown = sorted(set(filters) - allowed)
    if not unknown:
        return

    details = []
    for key in unknown:
        alias = _PYTHON_FILTER_ALIASES.get(key)
        if alias is not None and alias in allowed:
            details.append(
                f"{key!r} (the Python API spelling; use {alias!r} here)"
            )
        else:
            details.append(repr(key))

    raise ValueError(
        f"Unknown filter(s) for mode={mode!r}: {', '.join(details)}. "
        f"This mode reads: {', '.join(sorted(allowed))}. "
        "Filters it does not read are ignored, which would emit a manifest "
        "that does not match the request."
    )


def _resource_from_dict(mapping: Mapping[str, Any], cfg: Config) -> R3Resource:
    """Create an :class:`R3Resource` from a manifest mapping.

    The mapping may contain convenience keys like ``url`` and ``arcname``.
    These are ignored during rehydration; the canonical URL and archive name
    are derived from the description and the current :class:`Config`.

    Args:
      mapping: JSON-like mapping for a single resource
        (often a line from JSONL).
      cfg: :class:`Config` for URL construction and caching behavior.

    Returns:
      A configured :class:`R3Resource` instance.

    Raises:
      KeyError: If ``resource_type`` is missing.
      ValueError: If the resource fields are invalid for that type.
    """
    return R3Resource.from_mapping(mapping, config=cfg)


def _write_jsonl(
    resources: Iterable[R3Resource], output_path: Path | None
) -> None:
    """Write one JSON object per line describing each resource.

    The JSON schema contains all dataclass fields from the description plus
    derived ``url`` and ``arcname`` fields so the manifest is self-contained.

    Args:
      resources: Iterable of resources to serialize.
      output_path: Optional file path; when ``None``, writes to stdout.
    """
    sink: Any
    close = False
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        sink = output_path.open("w", encoding="utf-8")
        close = True
    else:
        sink = sys.stdout

    try:
        for res in resources:
            desc = res.description
            body = (
                dataclasses.asdict(desc)
                if dataclasses.is_dataclass(desc)
                else {}
            )
            body["url"] = res.url
            body["arcname"] = res.arcname  # type: ignore[attr-defined]
            sink.write(json.dumps(body, ensure_ascii=False) + "\n")
    finally:
        if close:
            sink.close()


def _write_tsv(
    resources: Iterable[R3Resource], output_path: Path | None
) -> None:
    """Write a TSV manifest with a stable column ordering.

    Columns:
      resource_type, organism, data_source, genomic_unit, project, sample,
      table_name, junction_type, annotation_extension,
      junction_extension, url, arcname

    Missing fields are left blank.

    Args:
      resources: Iterable of resources to serialize.
      output_path: Optional file path; when ``None``, writes to stdout.
    """
    order = [
        "resource_type",
        "organism",
        "data_source",
        "genomic_unit",
        "project",
        "sample",
        "table_name",
        "junction_type",
        "annotation_extension",
        "junction_extension",
        "url",
        "arcname",
    ]

    sink: Any
    close = False
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        sink = output_path.open("w", encoding="utf-8", newline="")
        close = True
    else:
        sink = sys.stdout

    try:
        sink.write("\t".join(order) + "\n")
        for res in resources:
            desc = res.description
            row = (
                dataclasses.asdict(desc)
                if dataclasses.is_dataclass(desc)
                else {}
            )
            row["url"] = res.url
            row["arcname"] = res.arcname  # type: ignore[attr-defined]
            sink.write(
                "\t".join(str(row.get(k, "") or "") for k in order) + "\n"
            )
    finally:
        if close:
            sink.close()


def _iter_manifest(path_or_dash: str, cfg: Config) -> Iterator[R3Resource]:
    """Yield :class:`R3Resource` objects from a JSONL manifest.

    Args:
      path_or_dash: Path to a JSONL file or '-' to read from stdin.
      cfg: :class:`Config` to attach to each resource.

    Yields:
      :class:`R3Resource` objects created from each JSON line.

    Raises:
      FileNotFoundError: If the manifest path does not exist.
      json.JSONDecodeError: For malformed JSON input.
      KeyError/ValueError: For invalid manifest objects.
    """
    if path_or_dash == "-":
        stream = sys.stdin
    else:
        p = Path(path_or_dash)
        with p.open("r", encoding="utf-8") as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                obj = json.loads(line)
                yield _resource_from_dict(obj, cfg)
        return

    for line in stream:
        line = line.strip()
        if not line:
            continue
        obj = json.loads(line)
        yield _resource_from_dict(obj, cfg)


def _cmd_ids(args: argparse.Namespace, cfg: Config) -> int:
    """Implement the ``ids`` subcommand.

    Writes two lists (samples and projects) to the specified files or stdout.

    Args:
      args: Parsed CLI arguments for the ``ids`` subcommand.
      cfg: Unused. Present for handler uniformity.

    Returns:
      Process exit code (0 for success).
    """
    del cfg  # Unused in this handler.

    from recount3.search import (  # pylint: disable=import-outside-toplevel
        create_sample_project_lists,
    )

    samples, projects = create_sample_project_lists(organism=args.organism)

    def _emit(items: list[str], path: str | None, label: str) -> None:
        if path:
            p = Path(path)
            p.parent.mkdir(parents=True, exist_ok=True)
            with p.open("w", encoding="utf-8") as fh:
                fh.write("\n".join(items))
            logging.info("Wrote %s list to: %s", label, p)
        else:
            for item in items:
                sys.stdout.write(item + "\n")

    _emit(samples, args.samples_out, "samples")
    _emit(projects, args.projects_out, "projects")
    return 0


def _cmd_search(args: argparse.Namespace, cfg: Config) -> int:
    """Execute the 'search' subcommand.

    Dispatches to one of the supported modes and writes a manifest to
    stdout or to a file when --output is provided.

    Supported modes:
      annotations   : List supported annotation file extensions.
      gene-exon     : Gene/exon count files for a filter set.
      junctions     : Junction artifacts (MM/RR/ID) for a filter set.
      metadata      : Project-level metadata tables.
      bigwig        : Per-sample coverage files.
      project       : Validate a project and enumerate all related files.
      sources       : List available data sources for an organism.
      source-meta   : Data-source-level metadata tables.

    Project mode filters:
      Required: organism, data_source, project
      Optional:
        genomic_unit=gene,exon
        annotation=default|all|gencode_v26,gencode_v29
          Human-readable name or convenience alias that is
          resolved to extension codes internally.  Accepts
          'default' (organism primary), 'all' (every
          available), a comma list of annotation names
          (e.g. gencode_v26), or raw codes (e.g. G026).
        annotation_extension=G026,G029
          Raw annotation file-extension codes.  When set,
          overrides 'annotation' completely.
        junction_type=ALL
        junction_extension=MM,RR,ID
        include_metadata=true|false
        include_bigwig=true|false

    Notes:
      * The values for comma-separated filters (for example, genomic_unit
        and junction_extension) are parsed as CSV.
      * Boolean filters accept 1/0, true/false, t/f, yes/no, y/n, on/off.
      * The emitted manifest can be consumed by 'recount3 download'.

    Args:
      args: Parsed command-line arguments for the 'search' subcommand.
      cfg: Global configuration (I/O settings, concurrency, cache).

    Returns:
      Zero on success; non-zero on error.

    Raises:
      ValueError: If required filters are missing or invalid for a mode.
    """
    filters = _parse_filters(args.filters)
    mode = args.mode
    _reject_unknown_filters(mode, filters)

    def _require(*keys: str) -> None:
        missing = [k for k in keys if k not in filters]
        if missing:
            missing_str = ", ".join(missing)
            raise ValueError(
                f"Missing required filters for mode={mode!r}: {missing_str}"
            )

    found: list[R3Resource]
    match mode:
        case "annotations":
            _require("organism", "genomic_unit", "annotation_extension")
            found = r3_search.search_annotations(
                organism=filters["organism"],
                genomic_unit=filters["genomic_unit"],
                annotation_extension=filters["annotation_extension"],
            )

        case "gene-exon":
            _require("organism", "data_source", "genomic_unit", "project")
            ann_ext = filters.get("annotation_extension", ("G026",))
            found = r3_search.search_count_files_gene_or_exon(
                organism=filters["organism"],
                data_source=filters["data_source"],
                genomic_unit=filters["genomic_unit"],
                project=filters["project"],
                annotation_extension=ann_ext,
            )

        case "junctions":
            _require("organism", "data_source", "project")
            found = r3_search.search_count_files_junctions(
                organism=filters["organism"],
                data_source=filters["data_source"],
                project=filters["project"],
                junction_type=filters.get("junction_type", "ALL"),
                junction_extension=filters.get("junction_extension", "MM"),
            )

        case "metadata":
            _require("organism", "data_source", "table_name", "project")
            found = r3_search.search_metadata_files(
                organism=filters["organism"],
                data_source=filters["data_source"],
                table_name=filters["table_name"],
                project=filters["project"],
            )

        case "bigwig":
            _require("organism", "data_source", "project", "sample")
            found = r3_search.search_bigwig_files(
                organism=filters["organism"],
                data_source=filters["data_source"],
                project=filters["project"],
                sample=filters["sample"],
            )

        case "project":
            _require("organism", "data_source", "project")

            def _as_bool(s: str | None, default: bool = False) -> bool:
                if s is None:
                    return default
                return s.lower() in ("1", "true", "t", "yes", "y", "on")

            def _csv_or_default(
                s: str | None, default: tuple[str, ...]
            ) -> tuple[str, ...]:
                if not s:
                    return default
                return tuple(p.strip() for p in s.split(",") if p.strip())

            gu = _csv_or_default(filters.get("genomic_unit"), ("gene", "exon"))
            # 'annotation' is a human-readable name or alias ("default",
            # "all", "gencode_v26", or a raw code like "G026").
            # 'annotation_extension' is the raw file-extension code
            # (e.g. "G026").  When both are provided the explicit
            # extension wins, mirroring the library's two-tier API
            # (see search._resolve_annotation_exts).
            annotations = filters.get("annotation", "default")
            ann_ext = _csv_or_default(
                filters.get("annotation_extension"), tuple()
            )
            jext = _csv_or_default(filters.get("junction_extension"), ("MM",))
            jtype = filters.get("junction_type", "ALL")
            inc_meta = _as_bool(filters.get("include_metadata"), True)
            inc_bw = _as_bool(filters.get("include_bigwig"), False)

            found = r3_search.search_project_all(
                organism=filters["organism"],
                data_source=filters["data_source"],
                project=filters["project"],
                genomic_units=gu,
                annotations=ann_ext if ann_ext else annotations,
                junction_type=jtype,
                junction_extension=jext,
                include_metadata=inc_meta,
                include_bigwig=inc_bw,
            )

        case "sources":
            _require("organism")
            found = r3_search.search_data_sources(organism=filters["organism"])

        case "source-meta":
            _require("organism", "data_source")
            found = r3_search.search_data_source_metadata(
                organism=filters["organism"],
                data_source=filters["data_source"],
            )

        case _:
            raise ValueError(f"Unknown search mode: {mode!r}")

    configured = [dataclasses.replace(r, config=cfg) for r in found]

    out_path: Path | None
    if args.output:
        out_path = Path(args.output)
    elif args.outdir:
        ext = "jsonl" if args.format == "jsonl" else "tsv"
        ts = datetime.now().strftime("%Y%m%d-%H%M%S")
        fname = f"{args.mode}-{ts}.{ext}"
        out_path = Path(args.outdir) / fname
    else:
        out_path = None

    if args.format == "jsonl":
        _write_jsonl(configured, out_path)
    else:
        _write_tsv(configured, out_path)

    if out_path is None:
        logging.info("Emitted %d resources to stdout.", len(configured))
    else:
        logging.info("Emitted %d resources to: %s", len(configured), out_path)
    return 0


def _download_one(
    res: R3Resource,
    cfg: Config,
    dest: Path,
    cache_mode: CacheMode,
    overwrite: bool,
) -> dict[str, Any]:
    """Download one resource and return a JSON-serializable event.

    Args:
      res: The resource to download.
      cfg: The resolved :class:`Config`.
      dest: Destination directory or .zip file path.
      cache_mode: Cache behavior to use.
      overwrite: Whether to overwrite an existing file in directory mode.

    Returns:
      A dict event with keys: ``url``, ``status``, ``dest`` and optionally
      ``error``. ``status`` is one of ``ok``, ``skipped``, or ``error``.
    """
    res = dataclasses.replace(res, config=cfg)

    is_zip = dest.suffix.lower() == ".zip"
    will_skip = False
    dest_file: Path | None = None
    if not is_zip:
        # Same naming rule as resource.download(): basename of url_path().
        dest_file = dest / Path(res.description.url_path()).name
        if dest_file.exists() and not overwrite:
            will_skip = True

    try:
        if will_skip:
            # Ensure cache (if desired) but avoid copy.
            if cache_mode != "disable":
                res.download(path=None, cache_mode=cache_mode)
            status = "skipped"
            out_path = str(dest_file) if dest_file else None
        else:
            out_path = res.download(
                path=dest,
                cache_mode=cache_mode,
                overwrite=overwrite,
                chunk_size=cfg.chunk_size,
            )
            status = "ok"
        return {"url": res.url, "status": status, "dest": out_path}
    except Exception as exc:  # pylint: disable=broad-exception-caught
        return {
            "url": res.url,
            "status": "error",
            "dest": None,
            "error": repr(exc),
        }


def _cmd_download(args: argparse.Namespace, cfg: Config) -> int:
    """Implement the ``download`` subcommand.

    Prints one JSONL progress event per attempted download to stdout. The exit
    code is:
      0 - all OK,
      1 - malformed ``--inline`` JSON,
      2 - fatal I/O or setup problems,
      3 - partial failures.

    Args:
      args: Parsed CLI arguments for the ``download`` subcommand.
      cfg: :class:`Config` for networking and caching.

    Returns:
      Process exit code.
    """
    if args.inline:
        try:
            one = json.loads(args.inline)
        except json.JSONDecodeError as exc:
            logging.error("Bad --inline JSON: %r", exc)
            return 1
        resources = [_resource_from_dict(one, cfg)]
    else:
        resources = list(_iter_manifest(args.manifest, cfg))

    if not resources:
        logging.warning("No resources to download.")
        return 0

    dest = Path(args.dest)
    cache_mode: CacheMode = args.cache  # type: ignore[assignment]
    overwrite: bool = bool(args.overwrite)

    if dest.suffix.lower() == ".zip":
        dest.parent.mkdir(parents=True, exist_ok=True)
    else:
        dest.mkdir(parents=True, exist_ok=True)

    total = len(resources)
    errors = 0

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = [
            pool.submit(_download_one, r, cfg, dest, cache_mode, overwrite)
            for r in resources
        ]
        for fut in concurrent.futures.as_completed(futs):
            evt = fut.result()
            sys.stdout.write(json.dumps(evt, ensure_ascii=False) + "\n")
            if evt.get("status") == "error":
                errors += 1

    if errors == 0:
        logging.info(
            "Downloaded %d/%d successfully to: %s",
            total,
            total,
            dest,
        )
        return 0
    if errors < total:
        logging.warning(
            "Partial success: %d/%d failed. Destination: %s",
            errors,
            total,
            dest,
        )
        return 3
    logging.error(
        "All downloads failed: %d/%d errors. Destination: %s",
        errors,
        total,
        dest,
    )
    return 2


def _cmd_bundle_stack_counts(args: argparse.Namespace, cfg: Config) -> int:
    """Implement ``bundle stack-counts`` subcommand.

    Loads resources from a manifest and concatenates compatible count matrices.

    Parquet output is validated before any download happens, so a missing
    engine fails in a second rather than after assembling the matrix.

    Args:
      args: Parsed CLI arguments for ``bundle stack-counts``.
      cfg: :class:`Config` used for resource loading.

    Returns:
      Process exit code (0 for success, non-zero for failure).
    """
    out = Path(args.out)
    to_parquet = out.suffix.lower() == ".parquet"

    if to_parquet:
        try:
            engine = _utils.ensure_parquet_engine()
        except ImportError as exc:
            logging.error("%s", exc)
            return 2
        logging.debug("Parquet engine resolved by pandas: %s", engine)

    resources = list(_iter_manifest(args.manifest, cfg))
    bundle = R3ResourceBundle()
    bundle.extend(resources)

    compat: CompatibilityMode = args.compat  # type: ignore[assignment]
    try:
        df = bundle.stack_count_matrices(
            join_policy=args.join,
            axis=int(args.axis),
            verify_integrity=bool(args.verify_integrity),
            autoload=True,
            compat=compat,
        )
    except (CompatibilityError, LoadError, ValueError, ImportError) as exc:
        logging.error("Failed to stack count matrices (reason: %r)", exc)
        return 2

    if to_parquet:
        sparse_columns = _utils.sparse_column_names(df)
        if sparse_columns:
            if not args.densify:
                logging.error(
                    "Cannot write Parquet: %d of %d columns use a pandas "
                    "sparse dtype (for example, %s), which no Parquet engine "
                    "supports. Junction count matrices are sparse-backed. "
                    "Re-run with --densify to materialize every implicit "
                    "zero (this can need far more memory than the sparse "
                    "matrix), or write .tsv/.tsv.gz/.csv instead, which "
                    "handles sparse columns without densifying the whole "
                    "frame in memory.",
                    len(sparse_columns),
                    df.shape[1],
                    ", ".join(sparse_columns[:3]),
                )
                return 2
            logging.warning(
                "Densifying %d sparse column(s) for Parquet output; "
                "memory use scales with the full %d x %d matrix.",
                len(sparse_columns),
                df.shape[0],
                df.shape[1],
            )
            df = _utils.densify_sparse_columns(df)

    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        if to_parquet:
            df.to_parquet(out)
        else:
            sep = "," if out.suffix.lower() == ".csv" else "\t"
            df.to_csv(out, sep=sep)
    except Exception as exc:  # pylint: disable=broad-exception-caught
        logging.error("Failed to write output (reason: %r): %s", exc, out)
        return 2

    logging.info("Wrote stacked table to: %s", out)
    return 0


def _cmd_bundle_se(args: argparse.Namespace, cfg: Config) -> int:
    """Implement ``bundle se`` subcommand.

    Loads resources from a manifest and assembles a
    :class:`~summarizedexperiment.SummarizedExperiment`.

    Args:
      args: Parsed CLI arguments for ``bundle se``.
      cfg: :class:`Config` used for resource loading.

    Returns:
      Process exit code:
        0 - success,
        2 - fatal error (missing dependency, build failure, or write failure).
    """
    out = Path(args.out)
    to_h5ad = out.suffix.lower() == ".h5ad"

    if to_h5ad:
        try:
            _utils.ensure_anndata_support()
        except ImportError as exc:
            logging.error("Cannot write .h5ad: %s", exc)
            return 2

    resources = list(_iter_manifest(args.manifest, cfg))
    bundle = R3ResourceBundle()
    bundle.extend(resources)

    try:
        se = bundle.to_summarized_experiment(
            genomic_unit=args.genomic_unit,
            annotation_extension=args.annotation,
            assay_name=args.assay_name,
            join_policy=args.join,
            autoload=True,
        )
    except ImportError as exc:
        logging.error("Missing optional dependency: %s", exc)
        return 2
    except Exception as exc:  # pylint: disable=broad-exception-caught
        logging.error("Failed to build SE (reason: %r).", exc)
        return 2

    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        if to_h5ad:
            adata = _utils.experiment_to_anndata(se)
            cast = _utils.normalize_anndata_for_hdf5(adata)
            if cast:
                logging.info(
                    "Cast %d all-missing column(s) to NaN for HDF5: %s",
                    len(cast),
                    ", ".join(cast),
                )
            unsafe = _utils.hdf5_unsafe_column_names(adata)
            unsafe += _utils.hdf5_unsafe_uns_keys(adata)
            if unsafe:
                if not args.sanitize_columns:
                    logging.error(
                        "Cannot write .h5ad: %d name(s) contain a "
                        "forward slash, which HDF5 reads as a path separator "
                        "(for example, %s). recount3 STAR QC fields are named "
                        "after splice motifs, so most projects with SRA "
                        "metadata hit this, in the sample columns and in the "
                        "uns provenance map keyed by them. Re-run with "
                        "--sanitize-columns to replace each '/' with '_' "
                        "(this renames the columns your analysis indexes by), "
                        "or write a .pkl instead, which keeps the names "
                        "verbatim.",
                        len(unsafe),
                        ", ".join(unsafe[:3]),
                    )
                    return 2
                renames = _utils.sanitize_anndata_column_names(adata)
                renames += _utils.sanitize_anndata_uns_keys(adata)
                logging.warning(
                    "Renamed %d name(s) for HDF5: %s",
                    len(renames),
                    ", ".join(
                        f"{before} -> {after}" for before, after in renames
                    ),
                )
            adata.write_h5ad(out)
        else:
            with open(out, "wb") as fh:
                pickle.dump(se, fh, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as exc:  # pylint: disable=broad-exception-caught
        logging.error("Failed to write output (reason: %r): %s", exc, out)
        return 2

    logging.info("Wrote SummarizedExperiment to: %s", out)
    return 0


def _cmd_bundle_rse(args: argparse.Namespace, cfg: Config) -> int:
    """Implement ``bundle rse`` subcommand.

    Loads resources from a manifest and assembles a
    :class:`~summarizedexperiment.RangedSummarizedExperiment`.

    Args:
      args: Parsed CLI arguments for ``bundle rse``.
      cfg: :class:`Config` used for resource loading.

    Returns:
      Process exit code:
        0 - success,
        2 - fatal error (missing dependency, build failure, or write failure).
    """
    out = Path(args.out)
    to_h5ad = out.suffix.lower() == ".h5ad"

    if to_h5ad:
        try:
            _utils.ensure_anndata_support()
        except ImportError as exc:
            logging.error("Cannot write .h5ad: %s", exc)
            return 2

    resources = list(_iter_manifest(args.manifest, cfg))
    bundle = R3ResourceBundle()
    bundle.extend(resources)

    try:
        rse = bundle.to_ranged_summarized_experiment(
            genomic_unit=args.genomic_unit,
            annotation_extension=args.annotation,
            assay_name=args.assay_name,
            join_policy=args.join,
            autoload=True,
            allow_fallback_to_se=bool(args.allow_fallback_to_se),
        )
    except ImportError as exc:
        logging.error("Missing optional dependency: %s", exc)
        return 2
    except Exception as exc:  # pylint: disable=broad-exception-caught
        logging.error("Failed to build RSE (reason: %r).", exc)
        return 2

    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        if to_h5ad:
            adata = _utils.experiment_to_anndata(rse)
            cast = _utils.normalize_anndata_for_hdf5(adata)
            if cast:
                logging.info(
                    "Cast %d all-missing column(s) to NaN for HDF5: %s",
                    len(cast),
                    ", ".join(cast),
                )
            unsafe = _utils.hdf5_unsafe_column_names(adata)
            unsafe += _utils.hdf5_unsafe_uns_keys(adata)
            if unsafe:
                if not args.sanitize_columns:
                    logging.error(
                        "Cannot write .h5ad: %d name(s) contain a "
                        "forward slash, which HDF5 reads as a path separator "
                        "(for example, %s). recount3 STAR QC fields are named "
                        "after splice motifs, so most projects with SRA "
                        "metadata hit this, in the sample columns and in the "
                        "uns provenance map keyed by them. Re-run with "
                        "--sanitize-columns to replace each '/' with '_' "
                        "(this renames the columns your analysis indexes by), "
                        "or write a .pkl instead, which keeps the names "
                        "verbatim.",
                        len(unsafe),
                        ", ".join(unsafe[:3]),
                    )
                    return 2
                renames = _utils.sanitize_anndata_column_names(adata)
                renames += _utils.sanitize_anndata_uns_keys(adata)
                logging.warning(
                    "Renamed %d name(s) for HDF5: %s",
                    len(renames),
                    ", ".join(
                        f"{before} -> {after}" for before, after in renames
                    ),
                )
            adata.write_h5ad(out)
        else:
            with open(out, "wb") as fh:
                pickle.dump(rse, fh, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as exc:  # pylint: disable=broad-exception-caught
        logging.error("Failed to write output (reason: %r): %s", exc, out)
        return 2

    logging.info("Wrote (Ranged)SummarizedExperiment to: %s", out)
    return 0


def _cmd_smoke_test(args: argparse.Namespace, cfg: Config) -> int:
    """Implement the ``smoke-test`` subcommand (small connectivity test).

    This downloads a handful of sources metadata files for 'human' (which are
    tiny), printing JSONL progress events to stdout like the main download
    command. Intended for CI or quick validation.

    Args:
      args: Parsed CLI arguments for the ``smoke-test`` subcommand.
      cfg: :class:`Config` for networking and caching.

    Returns:
      Process exit code.
    """
    limit = max(1, int(args.limit))
    resources = r3_search.search_data_source_metadata(
        organism="human", data_source="sra"
    )
    resources = [dataclasses.replace(r, config=cfg) for r in resources[:limit]]

    dest = Path(args.dest)
    dest.mkdir(parents=True, exist_ok=True)

    for res in resources:
        evt = _download_one(
            res, cfg, dest, cache_mode="enable", overwrite=False
        )
        sys.stdout.write(json.dumps(evt, ensure_ascii=False) + "\n")

    return 0


def _dispatch(args: argparse.Namespace, cfg: Config) -> int:
    """Dispatch parsed arguments to the appropriate subcommand handler.

    Args:
      args: Parsed CLI arguments.
      cfg: Resolved :class:`Config`.

    Returns:
      Process exit code.

    Raises:
      ValueError: If the command is unknown (should not happen with argparse).
    """
    match args.command:
        case "ids":
            return _cmd_ids(args, cfg)
        case "search":
            return _cmd_search(args, cfg)
        case "download":
            return _cmd_download(args, cfg)
        case "bundle":
            match args.bundle_cmd:
                case "stack-counts":
                    return _cmd_bundle_stack_counts(args, cfg)
                case "se":
                    return _cmd_bundle_se(args, cfg)
                case "rse":
                    return _cmd_bundle_rse(args, cfg)
                case _:
                    raise ValueError(
                        f"Unknown bundle subcommand: {args.bundle_cmd!r}"
                    )
        case "smoke-test":
            return _cmd_smoke_test(args, cfg)
        case _:
            raise ValueError(f"Unknown command: {args.command!r}")


def main(argv: list[str] | None = None) -> None:
    """Program entry point for the recount3 CLI.

    This function is deliberately small and import-safe. It sets up logging,
    builds the configuration, dispatches to a subcommand handler, and finally
    calls :func:`sys.exit` with the handler's return code.

    Args:
      argv: Optional list of arguments. If not provided, :data:`sys.argv[1:]`
        is used. Accepting an explicit list makes this function easy to test.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)

    _init_logging(args)
    cfg = _build_config_from_env_and_flags(args)

    try:
        code = _dispatch(args, cfg)
    except KeyboardInterrupt:
        logging.error("Interrupted.")
        code = 130
    except BrokenPipeError:
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        code = 0
    except (ConfigurationError, Recount3Error, ValueError) as exc:
        logging.error("Fatal error: %r", exc)
        code = 2

    sys.exit(code)


if __name__ == "__main__":
    main()
