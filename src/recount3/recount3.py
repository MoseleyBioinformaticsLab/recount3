#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Utilities for describing, downloading, and loading recount3 resources.

This module provides a small, typed API for interacting with the recount3 data
repository. It focuses on three core tasks:

1) **Describe** a resource with :class:`R3ResourceDescription`
2) **Fetch/Materialize** it with :class:`R3Resource.download`
3) **Load** tabular resources with :class:`R3Resource.load`

It also includes convenience search helpers (e.g., :func:`search_annotations`)
and a lightweight :class:`Recount3Bundle` to batch-download and combine results.

Overview
--------
- **Typed URL building:** `R3ResourceDescription` validates the fields required
  for each data type and constructs a stable URL path (duffel layout).
- **Deterministic materialization:** `R3Resource` derives a canonical ZIP
  `arcname` from the URL path and uses safe linking/copying for local files.
- **Explicit caching:** `download(cache_mode=...)` controls whether the cache is
  used, bypassed, or refreshed.
- **Simple loading:** `load()` parses TSV/TSV.GZ/MD/MD.GZ into a
  `pandas.DataFrame` (if pandas is available) or a list of dicts.
- **Discovery helpers:** `search_*` functions generate `R3Resource` objects for
  common queries; `create_sample_project_lists()` scans metadata to derive
  sample and project IDs.
- **Batch workflows:** `Recount3Bundle` downloads, caches loaded objects, and
  can concatenate multiple count matrices.

Supported resource types
------------------------
The `resource_type` field of :class:`R3ResourceDescription` accepts:

- `"annotations"`: GTF-like gene/exon annotation files.
- `"count_files_gene_or_exon"`: per-project gene/exon count tables.
- `"count_files_junctions"`: per-project junction count tables.
- `"metadata_files"`: per-project metadata tables.
- `"bigwig_files"`: per-sample BigWig files (download only; `load()` is
  deprecated for BigWig).
- `"data_sources"`: organism-level data source index (duffel `homes_index`).
- `"data_source_metadata"`: top-level metadata table per data source.

Caching and materialization
---------------------------
`R3Resource.download(path, cache_mode)` controls how bytes are obtained
and where they go:

- `cache_mode="enable"`: Use cache if present; otherwise download into cache.
- `cache_mode="disable"`: Bypass cache and stream directly to the destination.
- `cache_mode="update"`: Re-download into cache (atomic replace), then behave
  like `"enable"`.

`path` semantics:
- `None`: cache-only (ensure bytes exist in cache; no filesystem output).
- Directory path: materialize by hardlink (or copy) using the URL basename.
- `*.zip`: add an entry whose name matches the URL path (without the leading
  slash) for deterministic ZIP contents.

Environment variables
---------------------
- ``RECOUNT3_URL``: Base URL (default: `http://duffel.rail.bio/recount3/`).
- ``RECOUNT3_CACHE_DIR``: Cache directory (default: `~/.cache/recount3/files`).
- ``RECOUNT3_CACHE_DISABLE``: If `"1"`, disable cache behavior globally.
- ``RECOUNT3_HTTP_TIMEOUT``: Network timeout in seconds (default: `60`).
- ``RECOUNT3_MAX_RETRIES``: Max retry attempts for HTTP downloads (default: `3`).
- ``RECOUNT3_INSECURE_SSL``: If `"1"`, disable TLS verification (not recommended).
- ``RECOUNT3_USER_AGENT``: Custom HTTP User-Agent string.

Exceptions
----------
- :class:`DeprecationError`: Raised by `R3Resource.load()` for `"bigwig_files"`.
- :class:`ValueError`: Invalid parameter combinations or path semantics.
- Network errors from `urllib` (e.g., :class:`URLError`, :class:`HTTPError`).

Examples
--------
Basic usage:

    >>> desc = R3ResourceDescription(
    ...     resource_type="annotations",
    ...     organism="human",
    ...     genomic_unit="gene",
    ...     annotation_file_extension="G026",
    ... )
    >>> res = R3Resource(desc)
    >>> # Ensure present in cache (no filesystem output):
    >>> res.download(path=None, cache_mode="enable")
    >>> # Materialize into a directory:
    >>> out_path = res.download(path="./data", cache_mode="enable")
    >>> # Stream directly into a ZIP (bypass cache):
    >>> res.download(path="./bundle.zip", cache_mode="disable")
    >>> # Load a tabular file:
    >>> df = res.load()

Batch discovery and stacking:

    >>> resources = search_count_files_gene_or_exon(
    ...     organism="human",
    ...     data_source="sra",
    ...     genomic_unit="gene",
    ...     project=["SRP000001", "SRP000002"],
    ...     annotation_file_extension="G026",
    ... )
    >>> bundle = materialize_bundle(resources, destination_directory="counts", load=True)
    >>> combined = bundle.stack_count_matrices(join="inner")

Command-line behavior
---------------------
When executed as a script, this module:
1) Derives human `samplist.txt` and `projlist.txt` via
   :func:`create_sample_project_lists`.
2) Runs a small smoke test that prints example URLs and attempts downloads
   into the current directory.

Dependencies
------------
- `pandas` for DataFrame loading
"""

from __future__ import annotations

import contextlib
import csv
import errno
import gzip
import hashlib
import itertools
import io
import os
import shutil
import socket
import threading
import time
import traceback
import typing as _t
import urllib.error
import urllib.parse
import urllib.request
import ssl
import zipfile
import dataclasses
from pathlib import Path
from typing import Literal, Optional

# ---------------------------------------------------------------------------
# Required dependencies
# ---------------------------------------------------------------------------

import pandas as pd

# ---------------------------------------------------------------------------
# Public constants
# ---------------------------------------------------------------------------

BASE_URL = os.environ.get("RECOUNT3_URL", "http://duffel.rail.bio/recount3/").rstrip("/") + "/"
DEFAULT_TIMEOUT = int(os.environ.get("RECOUNT3_HTTP_TIMEOUT", "60"))
INSECURE_SSL = os.environ.get("RECOUNT3_INSECURE_SSL", "0") == "1"
MAX_RETRIES = int(os.environ.get("RECOUNT3_MAX_RETRIES", "3"))
USER_AGENT = os.environ.get("RECOUNT3_USER_AGENT") or "recount3-python/0.2 (+https://github.com/MoseleyBioinformaticsLab/recount3)"

# Cache root may be overridden by env var. Nest files under "files/" to
# avoid collisions if additional metadata directories are added later.
DEFAULT_CACHE_DIR = os.environ.get(
    "RECOUNT3_CACHE_DIR",
    os.path.join(os.path.expanduser("~"), ".cache", "recount3", "files"),
)
CACHE_DISABLED = os.environ.get("RECOUNT3_CACHE_DISABLE", "0") == "1"

# Internal validation enums.
_VALID_ORGANISMS = {"human", "mouse"}
_VALID_DATA_SOURCES = {"sra", "gtex", "tcga"}
_VALID_GENOMIC_UNITS = {"gene", "exon"}

# Serialize writes to the same cache path.
_file_lock = threading.Lock()

# Chunk size for streaming network/file copies.
_DEFAULT_CHUNK_SIZE = 1024 * 1024

# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class DeprecationError(RuntimeError):
    """Raised when an old feature is invoked.
    
    Temporary, using to help with refactoring process.
    """


# ---------------------------------------------------------------------------
# Helpers: filesystem, HTTP, and cache utilities
# ---------------------------------------------------------------------------


def _ensure_dir(path: str | Path) -> None:
    """Ensure a directory exists, creating parents as needed.

    Args:
      path: Directory path to create.

    Raises:
      NotADirectoryError: If the path exists as a non-directory.
    """
    p = Path(path)
    if not p.exists():
        p.mkdir(parents=True, exist_ok=True)
    elif not p.is_dir():
        raise NotADirectoryError(str(p))


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _cache_key_for_url(url: str) -> str:
    parsed = urllib.parse.urlparse(url)
    key = parsed.path.lstrip("/")
    digest = _sha256(url)[:16]
    return f"{digest}__{Path(key).name}"


def _cache_path(url: str, cache_root: str | Path) -> Path:
    _ensure_dir(cache_root)
    return Path(cache_root) / _cache_key_for_url(url)


def _ssl_insecure_context():
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def _http_open(
    url: str,
    *,
    timeout: float = DEFAULT_TIMEOUT,
    headers: dict[str, _t.Any] | None = None,
) -> _t.BinaryIO:
    """Open a URL with urllib, ensuring all headers are valid strings.

    Args:
      url: Absolute URL to open.
      timeout: Socket timeout (seconds) passed to ``opener.open(...)``.
      headers: Optional extra headers to send. Non-string values (including
        ``None``) are sanitized: ``None`` is dropped; others are ``str(...)``-cast.

    Returns:
      A file-like binary object for the HTTP response (e.g., an
      ``http.client.HTTPResponse``), suitable for use with ``contextlib.closing``.

    Raises:
      URLError, HTTPError, socket.timeout: Network/HTTP errors from urllib.
      ValueError: If ``url`` is empty.

    Environment:
      RECOUNT3_USER_AGENT: If set, overrides the default User-Agent string.
        You may also define a module-level ``USER_AGENT``; if present,
        it takes precedence over the env var.
    """
    if not url:
        raise ValueError("Empty URL passed to _http_open().")

    ua = (
        globals().get("USER_AGENT")
        or os.environ.get("RECOUNT3_USER_AGENT")
        or "recount3-python/0.3 (+https://github.com/MoseleyBioinformaticsLab/recount3)"
    )

    merged_headers: dict[str, str] = {"User-Agent": str(ua)}

    if headers:
        for k, v in headers.items():
            if v is None:
                continue
            merged_headers[str(k)] = str(v)

    opener = urllib.request.build_opener()
    if INSECURE_SSL:
        opener.add_handler(urllib.request.HTTPSHandler(context=_ssl_insecure_context()))  # pragma: no cover

    req = urllib.request.Request(url, headers=merged_headers)
    return opener.open(req, timeout=timeout)


def _with_retries(func, *, attempts: int = MAX_RETRIES, base_sleep: float = 0.5):
    last_exc: BaseException | None = None
    for i in range(max(1, attempts)):
        try:
            return func()
        except (urllib.error.URLError, socket.timeout) as exc:  # pragma: no cover - network
            last_exc = exc
            if i < attempts - 1:
                time.sleep(base_sleep * (2**i))
    if last_exc:
        raise last_exc


def _stream_copy(src_fh, dst_fh, *, chunk_size: int = _DEFAULT_CHUNK_SIZE) -> int:
    total = 0
    while True:
        chunk = src_fh.read(chunk_size)
        if not chunk:
            break
        dst_fh.write(chunk)
        total += len(chunk)
    return total


def _hardlink_or_copy(src: Path, dst: Path) -> None:
    """Materialize by hardlink, fallback to copy on cross-device/perm errors."""
    _ensure_dir(dst.parent)  # ensure destination directory exists
    try:
        if dst.exists():
            dst.unlink()
        os.link(src, dst)
    except OSError as e:
        if e.errno in (errno.EXDEV, errno.EPERM, errno.EACCES, errno.EMLINK):
            _ensure_dir(dst.parent)
            shutil.copy2(src, dst)
        else:  # pragma: no cover - unexpected
            raise


def _atomic_replace(src_tmp: Path, final_path: Path) -> None:
    _ensure_dir(final_path.parent)
    os.replace(src_tmp, final_path)


def _download_to_file(url: str, out_path: Path, *, chunk_size: int) -> None:
    def _do():
        with contextlib.closing(_http_open(url)) as resp:
            _ensure_dir(out_path.parent)
            tmp = out_path.with_suffix(out_path.suffix + ".downloading")
            with open(tmp, "wb") as fh:
                _stream_copy(resp, fh, chunk_size=chunk_size)
            _atomic_replace(tmp, out_path)

    _with_retries(_do)


def _download_stream_to_zip(
    url: str,
    zip_path: Path,
    arcname: str,
    *,
    chunk_size: int,
    overwrite: bool,
) -> None:
    """Stream URL content directly into a .zip file (no cache involvement)."""
    _ensure_dir(zip_path.parent)
    mode = "a" if zip_path.exists() else "w"
    with zipfile.ZipFile(zip_path, mode=mode, compression=zipfile.ZIP_DEFLATED) as zf:
        if not overwrite:
            try:
                zf.getinfo(arcname)
                return  # already present
            except KeyError:
                pass
        with contextlib.closing(_http_open(url)) as resp:
            with zf.open(arcname, "w") as zfh:
                _stream_copy(resp, zfh, chunk_size=chunk_size)


def _write_cached_file_to_zip(
    cached_file: Path,
    zip_path: Path,
    arcname: str,
    *,
    overwrite: bool,
) -> None:
    _ensure_dir(zip_path.parent)
    mode = "a" if zip_path.exists() else "w"
    with zipfile.ZipFile(zip_path, mode=mode, compression=zipfile.ZIP_DEFLATED) as zf:
        if not overwrite:
            try:
                zf.getinfo(arcname)
                return
            except KeyError:
                pass
        zf.write(cached_file, arcname)


# ---------------------------------------------------------------------------
# Multi-factory resource descriptions
# ---------------------------------------------------------------------------

CacheMode = Literal["enable", "disable", "update"]


def _p2(s: _t.Optional[str]) -> str:
    """Return the last two characters of ``s`` for duffel shard directories.

    - ``None`` or ``""``  -> ``""``
    - ``"A"``              -> ``"A"``
    - ``"SRP107565"``      -> ``"65"``
    """
    if not s:
        return ""
    return s[-2:]


@dataclasses.dataclass(slots=True)
class _R3CommonFields:
    """Uniform field schema shared by all description subclasses.

    The *multi-factory* base returns concrete subclasses that still expose all
    of these attributes so existing code paths remain source-compatible.
    """
    resource_type: str
    organism: Optional[str] = None
    data_source: Optional[str] = None
    genomic_unit: Optional[str] = None
    project: Optional[str] = None
    sample: Optional[str] = None
    annotation_file_extension: Optional[str] = None
    junction_type: Optional[str] = None
    junction_file_extension: Optional[str] = None
    table_name: Optional[str] = None


class R3ResourceDescription:
    """Abstract base class and multi-factory for recount3 resource descriptors.

    Usage:
      >>> desc = R3ResourceDescription(resource_type="annotations", organism="human", ...)

    The base class overrides :py:meth:`__new__` to return an instance of the
    registered subclass for the given ``resource_type``. Subclasses should
    inherit from :class:`_R3CommonFields` and :class:`R3ResourceDescription`,
    be declared as ``@dataclasses.dataclass(slots=True)``, and implement
    :py:meth:`url_path()` plus any needed validation in ``__post_init__``.
    """

    # Registry mapping resource_type -> subclass
    _TYPE_REGISTRY: dict[str, type["R3ResourceDescription"]] = {}

    # ----- multi-factory -------------------------------------------------

    def __new__(cls, *args, **kwargs):
        """Return an instance of the appropriate concrete subclass.

        If ``cls`` is already a subclass (direct instantiation), default object
        creation is used. When invoked as ``R3ResourceDescription(...)``,
        the ``resource_type`` argument is examined and a suitable subclass is
        instantiated instead.

        Raises:
          KeyError: Missing ``resource_type``.
          ValueError: Unknown ``resource_type``.
        """
        if cls is not R3ResourceDescription:
            return super().__new__(cls)

        # Accept either positional or keyword 'resource_type', but prefer kw.
        rtype = kwargs.get("resource_type")
        if rtype is None and args:
            rtype = args[0]  # positional resource_type in _R3CommonFields

        if not rtype:
            raise KeyError("resource_type is required to construct R3ResourceDescription.")

        subcls = cls._TYPE_REGISTRY.get(str(rtype))
        if subcls is None:
            raise ValueError(f"Unsupported resource_type: {rtype!r}")

        return super().__new__(subcls)  # type: ignore

    # ----- registration helper ------------------------------------------

    @classmethod
    def register_type(cls, resource_type: str):
        """Decorator to register a concrete subclass in the factory.

        Example:
          >>> @R3ResourceDescription.register_type("annotations")
          ... @dataclasses.dataclass(slots=True)
          ... class R3Annotations(_R3CommonFields, R3ResourceDescription):
          ...     def url_path(self) -> str: ...
        """
        def _decorator(subcls: type["R3ResourceDescription"]):
            cls._TYPE_REGISTRY[resource_type] = subcls
            setattr(subcls, "_RESOURCE_TYPE", resource_type)
            return subcls
        return _decorator

    # ----- shared validation helpers ------------------------------------

    def _require(self, *names: str) -> None:
        """Raise ``KeyError`` if any named attribute is missing or empty."""
        for n in names:
            if getattr(self, n, None) in (None, ""):
                raise KeyError(f"Missing required field: {n}")

    def _check_organism(self) -> None:
        if getattr(self, "organism", None) not in _VALID_ORGANISMS:
            raise ValueError(f"Invalid organism: {getattr(self, 'organism', None)!r}")

    def _check_data_source(self) -> None:
        if getattr(self, "data_source", None) not in _VALID_DATA_SOURCES:
            raise ValueError(f"Invalid data_source: {getattr(self, 'data_source', None)!r}")

    def _check_genomic_unit(self) -> None:
        if getattr(self, "genomic_unit", None) not in _VALID_GENOMIC_UNITS:
            raise ValueError(f"Invalid genomic_unit: {getattr(self, 'genomic_unit', None)!r}")

    # ----- abstract API --------------------------------------------------

    def url_path(self) -> str:  # pragma: no cover - implemented by subclasses
        """Return the URL path (including leading slash)."""
        raise NotImplementedError


# ----- concrete description subclasses -------------------------------------

@R3ResourceDescription.register_type("annotations")
@dataclasses.dataclass(slots=True)
class R3Annotations(_R3CommonFields, R3ResourceDescription):
    """Descriptor for gene/exon annotation tables."""

    def __post_init__(self) -> None:
        self._require("organism", "genomic_unit", "annotation_file_extension")
        self._check_organism()
        self._check_genomic_unit()

    def url_path(self) -> str:
        return (
            f"{self.organism}/annotations/{self.genomic_unit}_sums/"
            f"{self.organism}.{self.genomic_unit}_sums.{self.annotation_file_extension}.gtf.gz"
        )


@R3ResourceDescription.register_type("count_files_gene_or_exon")
@dataclasses.dataclass(slots=True)
class R3GeneOrExonCounts(_R3CommonFields, R3ResourceDescription):
    """Descriptor for per-project gene/exon count matrices."""

    def __post_init__(self) -> None:
        # Annotation extension is required on duffel mirrors.
        self._require("organism", "data_source", "genomic_unit", "project", "annotation_file_extension")
        self._check_organism()
        self._check_data_source()
        self._check_genomic_unit()

    def url_path(self) -> str:
        return (
            f"{self.organism}/data_sources/{self.data_source}/{self.genomic_unit}_sums/"
            f"{_p2(self.project)}/{self.project}/"
            f"{self.data_source}.{self.genomic_unit}_sums.{self.project}.{self.annotation_file_extension}.gz"
        )


@R3ResourceDescription.register_type("count_files_junctions")
@dataclasses.dataclass(slots=True)
class R3JunctionCounts(_R3CommonFields, R3ResourceDescription):
    """Descriptor for per-project junction count files."""

    def __post_init__(self) -> None:
        self._require("organism", "data_source", "project", "junction_type", "junction_file_extension")
        self._check_organism()
        self._check_data_source()

    def url_path(self) -> str:
        return (
            f"{self.organism}/data_sources/{self.data_source}/junctions/"
            f"{_p2(self.project)}/{self.project}/"
            f"{self.data_source}.junctions.{self.project}.{self.junction_type}.{self.junction_file_extension}.gz"
        )


@R3ResourceDescription.register_type("metadata_files")
@dataclasses.dataclass(slots=True)
class R3ProjectMetadata(_R3CommonFields, R3ResourceDescription):
    """Descriptor for per-project metadata tables."""

    def __post_init__(self) -> None:
        self._require("organism", "data_source", "project", "table_name")
        self._check_organism()
        self._check_data_source()

    def url_path(self) -> str:
        return (
            f"{self.organism}/data_sources/{self.data_source}/metadata/"
            f"{_p2(self.project)}/{self.project}/"
            f"{self.data_source}.{self.table_name}.{self.project}.MD.gz"
        )


@R3ResourceDescription.register_type("bigwig_files")
@dataclasses.dataclass(slots=True)
class R3BigWig(_R3CommonFields, R3ResourceDescription):
    """Descriptor for per-sample BigWig coverage files (download-only)."""

    def __post_init__(self) -> None:
        self._require("organism", "data_source", "project", "sample")
        self._check_organism()
        self._check_data_source()

    def url_path(self) -> str:
        return (
            f"{self.organism}/data_sources/{self.data_source}/base_sums/"
            f"{_p2(self.project)}/{self.project}/{_p2(self.sample)}/"
            f"{self.data_source}.base_sums.{self.project}_{self.sample}.ALL.bw"
        )


@R3ResourceDescription.register_type("data_sources")
@dataclasses.dataclass(slots=True)
class R3DataSources(_R3CommonFields, R3ResourceDescription):
    """Descriptor for the organism-level data-source index (homes_index)."""

    def __post_init__(self) -> None:
        self._require("organism")
        self._check_organism()

    def url_path(self) -> str:
        return f"{self.organism}/homes_index"


@R3ResourceDescription.register_type("data_source_metadata")
@dataclasses.dataclass(slots=True)
class R3DataSourceMetadata(_R3CommonFields, R3ResourceDescription):
    """Descriptor for source-level metadata listing tables."""

    def __post_init__(self) -> None:
        self._require("organism", "data_source")
        self._check_organism()
        self._check_data_source()

    def url_path(self) -> str:
        return (
            f"{self.organism}/data_sources/{self.data_source}/metadata/"
            f"{self.data_source}.recount_project.MD.gz"
        )


# ---------------------------------------------------------------------------
# R3Resource: orchestrates URL, caching, downloading, and loading
# ---------------------------------------------------------------------------


@dataclasses.dataclass(slots=True)
class R3Resource:
    """Resource that manages URL, caching, materialization, and loading.

    Attributes:
      description: :class:`R3ResourceDescription` used to build the path.
      url: Full URL derived from :data:`BASE_URL` + ``description.url_path()``.
      filepath: Optional local path when materialized to a directory (set by
        :meth:`download` when ``path`` is a directory).

    Methods:
      download(): Ensure bytes exist (cache-or-network) and optionally
        materialize to a directory or a zip archive.
      load(): Load the resource into a Python object based on its type.

    Cache modes:
      - ``"enable"``: use cache (create if missing).
      - ``"disable"``: bypass cache; stream directly to destination.
      - ``"update"``: re-download into cache under lock, then act like "enable".

    Path semantics:
      - ``None``: cache-only; do not materialize to filesystem.
      - Directory path: link/copy cached file, or download direct if cache
        is disabled.
      - ``.zip`` path: add entry to zip (arcname is URL path sans leading slash).
    """

    description: R3ResourceDescription
    url: str | None = None
    filepath: Optional[str] = None

    # ---- Derived properties ---------------------------------------------

    def __post_init__(self) -> None:
        if self.url is None:
            self.url = urllib.parse.urljoin(BASE_URL, self.description.url_path())

    @property
    def arcname(self) -> str:
        """Deterministic ZIP arcname derived from the URL path."""
        return self.description.url_path().lstrip("/")

    # ---- Cache helpers ---------------------------------------------------

    def _cache_root(self) -> Path:
        return Path(DEFAULT_CACHE_DIR)

    def _cached_path(self) -> Path:
        return _cache_path(self.url or "", self._cache_root())

    def _ensure_cached(self, *, mode: CacheMode, chunk_size: int) -> Path:
        """Ensure bytes are present in cache per mode; return cache path."""
        cache_path = self._cached_path()
        if mode == "disable":
            raise ValueError("Cache is not used in 'disable' mode")

        _ensure_dir(cache_path.parent)

        if mode == "enable":
            if cache_path.exists():
                return cache_path
            with _file_lock:
                if cache_path.exists():
                    return cache_path
                _download_to_file(self.url or "", cache_path, chunk_size=chunk_size)
            return cache_path

        if mode == "update":
            with _file_lock:
                _download_to_file(self.url or "", cache_path, chunk_size=chunk_size)
            return cache_path

        raise ValueError(f"Unknown cache mode: {mode!r}")

    # ---- Public API ------------------------------------------------------

    def download(
        self,
        path: Optional[str] = None,
        *,
        cache_mode: CacheMode = "enable",
        overwrite: bool = False,
        chunk_size: int = _DEFAULT_CHUNK_SIZE,
    ) -> Optional[str]:
        """Ensure bytes are available and optionally materialize.

        Args:
          path: ``None`` for cache-only; a directory to link/copy into; or a
            ``.zip`` archive to add the resource under its URL path.
          cache_mode: ``"enable"`` uses cache; ``"disable"`` streams to ``path``;
            ``"update"`` refreshes cache first then behaves like ``"enable"``.
          overwrite: Overwrite existing directory materialization when True.
          chunk_size: Chunk size in bytes for streaming/copying.

        Returns:
          Optional[str]: Destination file path if a directory materialization
          occurred; ``None`` for cache-only or zip additions.

        Raises:
          ValueError: For invalid combinations (e.g., ``path=None`` with
            ``cache_mode='disable'``) or invalid ``path`` semantics.
          urllib.error.URLError, socket.timeout: Network errors.
        """
        if cache_mode not in ("enable", "disable", "update"):
            raise ValueError(f"Invalid cache_mode: {cache_mode!r}")

        # Case 1: cache-only
        if path is None:
            if cache_mode == "disable":
                raise ValueError(
                    "path=None with cache_mode='disable' is invalid; "
                    "cache-only requires cache enabled."
                )
            self._ensure_cached(mode=cache_mode, chunk_size=chunk_size)
            return None

        # Determine whether path is a directory or a .zip file.
        path_p = Path(path)
        is_zip = path_p.suffix.lower() == ".zip"
        is_dir = path_p.suffix == "" or path_p.is_dir()

        if not (is_zip or is_dir):
            raise ValueError(f"path must be a directory or a .zip file, got: {path!r}")

        # Case 2: .zip archive materialization
        if is_zip:
            arcname = self.arcname
            if cache_mode == "disable":
                _download_stream_to_zip(
                    self.url or "",
                    path_p,
                    arcname,
                    chunk_size=chunk_size,
                    overwrite=overwrite,
                )
                return None
            cached = self._ensure_cached(mode=cache_mode, chunk_size=chunk_size)
            _write_cached_file_to_zip(cached, path_p, arcname, overwrite=overwrite)
            return None

        # Case 3: directory materialization
        # Decide final destination filename from URL path basename.
        dest_filename = Path(self.description.url_path()).name
        dest_path = path_p / dest_filename
        if cache_mode == "disable":
            # Stream directly to dest (no cache writes).
            _download_to_file(self.url or "", dest_path, chunk_size=chunk_size)
            self.filepath = str(dest_path)
            return str(dest_path)

        cached = self._ensure_cached(mode=cache_mode, chunk_size=chunk_size)
        if dest_path.exists():
            if not overwrite:
                self.filepath = str(dest_path)
                return str(dest_path)
            # Overwrite via hardlink/copy.
        _hardlink_or_copy(cached, dest_path)
        self.filepath = str(dest_path)
        return str(dest_path)

    # ------------------------------------------------------------------

    def load(self):
        """Load the resource into memory based on its type.

        * Tabular (``.tsv``/``.tsv.gz``) resources: return a ``pandas.DataFrame``
          when pandas is installed; otherwise return a list of dict rows.
        * BigWig resources: loading is **deprecated**; a :class:`DeprecationError`
          will be raised with guidance.

        Returns:
            Any: Loaded object depending on resource type.

        Raises:
            DeprecationError: For BigWig resource loading.
            ValueError: If the resource type is unsupported for loading.
            FileNotFoundError: If the cached file is missing after download.
        """
        if getattr(self.description, "resource_type", None) == "bigwig_files":
            raise DeprecationError(
                "Loading BigWig data via pyBigWig is deprecated. "
                "Please migrate to a BiocPy-compatible workflow or process "
                "BigWig files externally. You can still download or add them "
                "to zip archives with R3Resource.download()."
            )

        # Ensure bytes are present locally (cache enabled by default).
        self.download(path=None, cache_mode="enable")
        cached = self._cached_path()
        if not cached.exists():  # pragma: no cover - defensive
            raise FileNotFoundError(str(cached))

        # Load TSV/TSV.GZ into DataFrame or list-of-dicts.
        name = cached.name.lower()
        if name.endswith(".tsv") or name.endswith(".tsv.gz") or name.endswith(".md.gz"):
            if pd is not None:
                # pandas will infer compression and delimiter.
                return pd.read_table(cached, compression="infer")  # type: ignore[no-any-return]

            # Parse as TSV using Python.
            if name.endswith(".gz"):
                fh: io.TextIOBase
                fh = io.TextIOWrapper(gzip.open(cached, "rb"), encoding="utf-8")
            else:
                fh = open(cached, "rt", encoding="utf-8")
            with fh as tsv:
                reader = csv.DictReader(tsv, delimiter="\t")
                return list(reader)

        raise ValueError(
            f"Unsupported load() for resource type {getattr(self.description, 'resource_type', None)!r} "
            f"with file {cached.name!r}"
        )

    # ------------------------------------------------------------------

    def __repr__(self) -> str:  # pragma: no cover - representation only
        cls = self.__class__.__name__
        return f"{cls}(url={self.url!r}, arcname={self.arcname!r}, filepath={self.filepath!r})"


# ---------------------------------------------------------------------------
# Discovery/search helpers
# ---------------------------------------------------------------------------


_StringOrIterable = _t.Union[str, _t.Iterable[str]]


def _as_tuple(value: _StringOrIterable) -> _t.Tuple[str, ...]:
    """Normalize a string or iterable of strings into a tuple of strings."""
    if isinstance(value, str):
        return (value,)
    return tuple(value)


def _build_param_grid(resource_type: str, **required_values: _StringOrIterable) -> list[dict[str, str]]:
    """Cartesian product over provided values."""
    keys = list(required_values.keys())
    vals = [_as_tuple(required_values[k]) for k in keys]
    grid: list[dict[str, str]] = []
    for combo in itertools.product(*vals):
        d: dict[str, str] = {"resource_type": resource_type}
        d.update({k: v for k, v in zip(keys, combo)})
        grid.append(d)
    return grid


def _make_resources(
    param_dicts: _t.Iterable[dict[str, str]],
    *,
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Instantiate R3Resource for each param dict, with optional de-duplication."""
    out: list[R3Resource] = []
    seen: set[str] = set()
    for params in param_dicts:
        try:
            desc = R3ResourceDescription(**params)
            res = R3Resource(desc)
        except Exception:
            if strict:
                raise
            # skip invalid combos in non-strict mode
            continue
        if deduplicate:
            if res.url in seen:
                continue
            seen.add(res.url)  # type: ignore
        out.append(res)
    return out


# ---- Type-specific search wrappers -----------------------------------

def search_annotations(
    *,
    organism: _StringOrIterable,
    genomic_unit: _StringOrIterable,
    annotation_file_extension: _StringOrIterable,
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Return R3Resource objects for GTF/annotation files."""
    grid = _build_param_grid(
        "annotations",
        organism=organism,
        genomic_unit=genomic_unit,
        annotation_file_extension=annotation_file_extension,
    )
    return _make_resources(grid, strict=strict, deduplicate=deduplicate)


def search_count_files_gene_or_exon(
    *,
    organism: _StringOrIterable,
    data_source: _StringOrIterable,
    genomic_unit: _StringOrIterable,
    project: _StringOrIterable,
    annotation_file_extension: _StringOrIterable = ("G026",),
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Return R3Resource objects for per-project gene/exon counts.

    Note:
      ``annotation_file_extension`` is required on duffel mirrors. A default
      can be provided here if you routinely target a particular catalog (e.g.
      ``G026`` for Gencode v26). Override per call as needed.
    """
    grid = _build_param_grid(
        "count_files_gene_or_exon",
        organism=organism,
        data_source=data_source,
        genomic_unit=genomic_unit,
        project=project,
        annotation_file_extension=annotation_file_extension,
    )
    return _make_resources(grid, strict=strict, deduplicate=deduplicate)


def search_count_files_junctions(
    *,
    organism: _StringOrIterable,
    data_source: _StringOrIterable,
    project: _StringOrIterable,
    junction_type: _StringOrIterable = "ALL",
    junction_file_extension: _StringOrIterable = "MM",
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Return R3Resource objects for per-project junction count files."""
    grid = _build_param_grid(
        "count_files_junctions",
        organism=organism,
        data_source=data_source,
        project=project,
        junction_type=junction_type,
        junction_file_extension=junction_file_extension,
    )
    return _make_resources(grid, strict=strict, deduplicate=deduplicate)


def search_metadata_files(
    *,
    organism: _StringOrIterable,
    data_source: _StringOrIterable,
    table_name: _StringOrIterable,
    project: _StringOrIterable,
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Return R3Resource objects for top-level metadata tables."""
    grid = _build_param_grid(
        "metadata_files",
        organism=organism,
        data_source=data_source,
        table_name=table_name,
        project=project,
    )
    return _make_resources(grid, strict=strict, deduplicate=deduplicate)


def search_bigwig_files(
    *,
    organism: _StringOrIterable,
    data_source: _StringOrIterable,
    project: _StringOrIterable,
    sample: _StringOrIterable,
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Return R3Resource objects for BigWig resources (download-only; load() is deprecated)."""
    grid = _build_param_grid(
        "bigwig_files",
        organism=organism,
        data_source=data_source,
        project=project,
        sample=sample,
    )
    return _make_resources(grid, strict=strict, deduplicate=deduplicate)


def search_data_sources(
    *,
    organism: _StringOrIterable,
    strict: bool = True,
) -> list[R3Resource]:
    """Return an R3Resource representing the data-source index.

    The duffel layout requires an organism to select the correct homes index.
    """
    grid = _build_param_grid("data_sources", organism=organism)
    return _make_resources(grid, strict=strict, deduplicate=True)


def search_data_source_metadata(
    *,
    organism: _StringOrIterable,
    data_source: _StringOrIterable,
    strict: bool = True,
    deduplicate: bool = True,
) -> list[R3Resource]:
    """Return R3Resource objects for data-source level metadata listings."""
    grid = _build_param_grid(
        "data_source_metadata",
        organism=organism,
        data_source=data_source,
    )
    return _make_resources(grid, strict=strict, deduplicate=deduplicate)


# ---------------------------------------------------------------------------
# Aggregation / materialization conveniences
# ---------------------------------------------------------------------------

@dataclasses.dataclass(slots=True)
class Recount3Bundle:
    """Container for a set of R3Resource objects and (optionally) loaded data.

    - Holds the resources themselves.
    - Optionally caches the loaded tabular objects (DataFrames if pandas is
      available, else list-of-dicts).
    """
    resources: list[R3Resource] = dataclasses.field(default_factory=list)
    loaded: dict[str, list[object]] = dataclasses.field(default_factory=dict)

    def add(self, resource: R3Resource) -> None:
        """Add a resource to the bundle."""
        self.resources.append(resource)

    def load_all(self, *, strict: bool = True) -> "Recount3Bundle":
        """Load all resources and cache by resource_type."""
        for res in self.resources:
            try:
                obj = res.load()
            except Exception:
                if strict:
                    raise
                continue
            key = getattr(res.description, "resource_type", "unknown")
            self.loaded.setdefault(key, []).append(obj)
        return self

    def get_loaded(self, resource_type: str) -> list[object]:
        """Return previously loaded objects for a given resource type."""
        return self.loaded.get(resource_type, [])

    def stack_count_matrices(
        self,
        *,
        join: str = "inner",
        axis: int = 1,
        verify_integrity: bool = False,
        autoload: bool = True,
    ):
        """Concatenate count matrices (gene/exon or junction) as DataFrames.

        Args:
          join: 'inner' or 'outer' join on index.
          axis: 1 to stack by columns (default), 0 to stack by rows.
          verify_integrity: Passed to pandas.concat.
          autoload: If True, call load() for any count resources that aren't
            already cached in self.loaded.

        Returns:
          pandas.DataFrame

        Raises:
          ImportError: If pandas is unavailable.
          ValueError: If there are no count resources in the bundle.
        """
        if pd is None:
            raise ImportError("pandas is required for stack_count_matrices().")

        # Collect appropriate resources
        wanted = {
            "count_files_gene_or_exon",
            "count_files_junctions",
        }
        count_res = [r for r in self.resources if getattr(r.description, "resource_type", "") in wanted]
        if not count_res:
            raise ValueError("No count-file resources available to stack.")

        # Ensure loaded objects are present
        dfs: list[pd.DataFrame] = []
        for r in count_res:
            key = getattr(r.description, "resource_type", "unknown")
            cache = self.loaded.setdefault(key, [])
            df = None
            # try cache
            if cache:
                # consume one matching object if it's a DataFrame
                for i, obj in enumerate(cache):
                    if isinstance(obj, pd.DataFrame):
                        df = cache.pop(i)
                        break
            # load if needed
            if df is None:
                if not autoload:
                    raise ValueError(f"Resource {r.url} not loaded and autoload=False.")
                obj = r.load()
                if not isinstance(obj, pd.DataFrame):
                    raise TypeError(f"Loaded object for {r.url} is not a DataFrame.")
                df = obj
            dfs.append(df)  # type: ignore

        return pd.concat(  # type: ignore
            dfs,
            axis=axis,  # type: ignore
            join=join,  # type: ignore
            verify_integrity=verify_integrity,
        )


def materialize_resource(
    resource: R3Resource,
    *,
    destination_directory: str = ".",
    cache_mode: CacheMode = "enable",
    overwrite: bool = False,
) -> str:
    """Download a single resource into a directory and return its filepath."""
    out = resource.download(path=destination_directory, cache_mode=cache_mode, overwrite=overwrite)
    # R3Resource.download returns the on-disk path for directory targets.
    if out is None:
        raise RuntimeError("materialize_resource expected a directory path; got zip or cache-only.")
    return out


def materialize_bundle(
    resources: _t.Iterable[R3Resource],
    *,
    destination_directory: str = ".",
    cache_mode: CacheMode = "enable",
    strict: bool = True,
    load: bool = False,
    overwrite: bool = False,
) -> Recount3Bundle:
    """Download (and optionally load) a list of resources into a bundle.

    Args:
      resources: Iterable of R3Resource.
      destination_directory: Directory to materialize files via link/copy.
      cache_mode: Passed through to R3Resource.download.
      strict: If True, abort on first failure; else skip bad entries.
      load: If True, call R3Resource.load() and cache in the bundle.
      overwrite: Overwrite existing files when materializing.

    Returns:
      Recount3Bundle
    """
    bundle = Recount3Bundle()
    for r in resources:
        try:
            r.download(path=destination_directory, cache_mode=cache_mode, overwrite=overwrite)
            bundle.add(r)
            if load:
                obj = r.load()
                key = getattr(r.description, "resource_type", "unknown")
                bundle.loaded.setdefault(key, []).append(obj)
        except Exception:
            if strict:
                raise
            continue
    return bundle


# ---------------------------------------------------------------------------
# Discovery convenience: create sample & project lists
# ---------------------------------------------------------------------------

def create_sample_project_lists(organism: str = "") -> tuple[list[str], list[str]]:
    """Return (samples, projects) discovered from per-data-source metadata.

    This variant avoids the fragile global “sources index” endpoint and instead
    iterates the known data sources (sra, gtex, tcga). It loads each source's
    metadata table and extracts sample/project identifiers, optionally filtering
    by organism if the column is present.

    Args:
      organism: Optional organism filter (\"human\" or \"mouse\", case-insensitive).

    Returns:
      (samples, projects): Sorted unique IDs.
    """
    data_sources = sorted(_VALID_DATA_SOURCES)

    samples: set[str] = set()
    projects: set[str] = set()

    # Helper to normalize a loaded object (DataFrame or list-of-dicts) into rows of dict.
    def _rows(obj) -> list[dict]:
        if pd is not None and isinstance(obj, pd.DataFrame):
            return obj.to_dict("records")
        return list(obj)  # already list-of-dicts

    # Columns to look for (loose matching across mirrors)
    proj_keys = ("project", "project_id", "study", "study_accession")
    samp_keys = ("sample", "sample_id", "run", "run_accession", "external_id")
    org_keys = ("organism", "species")

    for ds in data_sources:
        try:
            res = R3Resource(R3ResourceDescription(resource_type="data_source_metadata", organism="human", data_source=ds))
            meta = res.load()
        except Exception:
            # If a particular data source isn’t present on the mirror, skip it.
            continue

        for row in _rows(meta):
            # Optional organism filter if a recognizable column exists.
            if organism:
                ov = None
                for k in org_keys:
                    if k in row and row[k]:
                        ov = str(row[k]).lower()
                        break
                if ov and ov != organism.lower():
                    continue

            # Project
            pv = None
            for k in proj_keys:
                if k in row and row[k]:
                    pv = str(row[k])
                    break
            if pv:
                projects.add(pv)

            # Sample
            sv = None
            for k in samp_keys:
                if k in row and row[k]:
                    sv = str(row[k])
                    break
            if sv:
                samples.add(sv)

    return sorted(samples), sorted(projects)


# =============================================================================
# Script entry point
# =============================================================================

def _example_descriptions() -> list[R3ResourceDescription]:
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
            annotation_file_extension="G026",   # <-- REQUIRED
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
        # 5) BigWig (download-only)
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
    """Smoke-test connectivity and download logic for a few exemplar resources.

    Prints each URL and attempts a download (to the current directory), catching
    and printing exceptions rather than aborting the run—matching the old behavior.
    """
    for desc in _example_descriptions():
        try:
            res = R3Resource(desc)
            print(res.url)
            folder_name = f"./downloads/{desc.resource_type}"  # type: ignore
            res.download(path=folder_name, cache_mode="enable")
        except Exception:
            print(traceback.format_exc())


def main() -> None:
    """Recreate the original script's behavior using the refactored API:

    - Build human sample/project lists and write them to samplist.txt/projlist.txt.
    - Run the example download smoke tests.
    """
    # Generate sample and project lists for human data
    samplist, projlist = create_sample_project_lists("human")

    # Write lists to text files (UTF-8, one ID per line)
    with open("samplist.txt", "w", encoding="utf-8") as samples_file:
        samples_file.write("\n".join(samplist))

    with open("projlist.txt", "w", encoding="utf-8") as projects_file:
        projects_file.write("\n".join(projlist))

    # Run download tests for the example resource set
    test_download()


if __name__ == "__main__":
    # Optional: very lightweight logging config for quick visibility
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s:%(name)s:%(message)s",
    )
    main()
