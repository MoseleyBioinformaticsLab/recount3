"""Resource orchestration: URL, caching, downloading, and loading.

This module mirrors the original `R3Resource` behavior, with configuration
threaded through (see :mod:`recount3.config`). No logic changes are made.
"""

from __future__ import annotations

import dataclasses
import threading
import urllib.parse
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd  # Required dependency in the original script.

from .config import Config, default_config
from .descriptions import R3ResourceDescription
from .errors import LoadError
from .types import CacheMode
from .bigwig import BigWigFile
from ._internal.cache import _cache_path
from ._internal.fs import _ensure_dir, _hardlink_or_copy
from ._internal.http import (
    download_stream_to_zip,
    download_to_file,
    write_cached_file_to_zip,
)

# Serialize writes to the same cache path (as in the original code).
_FILE_LOCK = threading.Lock()


@dataclass(slots=True)
class R3Resource:
    """Resource that manages URL, caching, materialization, and loading.

    Attributes:
      description: A :class:`R3ResourceDescription` used to build the path.
      url: Full URL derived from config.base_url + ``description.url_path()``.
      filepath: Optional local path when materialized to a directory.
      config: Optional :class:`Config`; defaults to :func:`default_config`.
    """

    description: R3ResourceDescription
    url: str | None = None
    filepath: Optional[str] = None
    config: Config | None = None

    _cached_data: object | None = dataclasses.field(default=None, init=False, repr=False)

    # ---- Derived properties --------------------------------------------

    def __post_init__(self) -> None:
        cfg = self.config or default_config()
        if self.url is None:
            self.url = urllib.parse.urljoin(cfg.base_url, self.description.url_path())

    @property
    def arcname(self) -> str:
        """Deterministic ZIP arcname derived from the URL path."""
        return self.description.url_path().lstrip("/")

    # ---- Cache helpers --------------------------------------------------

    def _cache_root(self) -> Path:
        cfg = self.config or default_config()
        return Path(cfg.cache_dir)

    def _cached_path(self) -> Path:
        return _cache_path(self.url or "", self._cache_root())

    def _ensure_cached(self, *, mode: CacheMode, chunk_size: int) -> Path:
        """Ensure bytes are present in cache per mode; return cache path."""
        cfg = self.config or default_config()
        cache_path = self._cached_path()
        if mode == "disable":
            raise ValueError("Cache is not used in 'disable' mode")

        _ensure_dir(cache_path.parent)

        if mode == "enable":
            if cache_path.exists():
                return cache_path
            with _FILE_LOCK:
                if cache_path.exists():
                    return cache_path
                download_to_file(
                    self.url or "",
                    cache_path,
                    chunk_size=chunk_size,
                    timeout=cfg.timeout,
                    insecure_ssl=cfg.insecure_ssl,
                    user_agent=cfg.user_agent,
                    attempts=cfg.max_retries,
                )
            return cache_path

        if mode == "update":
            with _FILE_LOCK:
                download_to_file(
                    self.url or "",
                    cache_path,
                    chunk_size=chunk_size,
                    timeout=cfg.timeout,
                    insecure_ssl=cfg.insecure_ssl,
                    user_agent=cfg.user_agent,
                    attempts=cfg.max_retries,
                )
            return cache_path

        raise ValueError(f"Unknown cache mode: {mode!r}")

    # ---- Public API -----------------------------------------------------

    def download(
        self,
        path: Optional[str] = None,
        *,
        cache_mode: CacheMode = "enable",
        overwrite: bool = False,
        chunk_size: int | None = None,
    ) -> Optional[str]:
        """Ensure bytes are available and optionally materialize.

        Args:
          path: ``None`` for cache-only; a directory to link/copy into; or a
            ``.zip`` archive to add the resource under its URL path.
          cache_mode: ``"enable"`` uses cache; ``"disable"`` streams to ``path``;
            ``"update"`` refreshes cache first then behaves like ``"enable"``.
          overwrite: Overwrite existing directory materialization when True.
          chunk_size: Chunk size in bytes for streaming/copying. Defaults to the
            configured value.

        Returns:
          Destination file path if a directory materialization occurred; ``None``
          for cache-only or zip additions.

        Raises:
          ValueError: For invalid combinations (e.g., ``path=None`` with
            ``cache_mode='disable'``) or invalid ``path`` semantics.
        """
        cfg = self.config or default_config()
        cs = chunk_size or cfg.chunk_size

        if cache_mode not in ("enable", "disable", "update"):
            raise ValueError(f"Invalid cache_mode: {cache_mode!r}")

        # Case 1: cache-only
        if path is None:
            if cache_mode == "disable":
                raise ValueError(
                    "path=None with cache_mode='disable' is invalid; "
                    "cache-only requires cache enabled."
                )
            self._ensure_cached(mode=cache_mode, chunk_size=cs)
            return None

        # Determine directory vs zip path.
        path_p = Path(path)
        is_zip = path_p.suffix.lower() == ".zip"
        is_dir = path_p.suffix == "" or path_p.is_dir()
        if not (is_zip or is_dir):
            raise ValueError(f"path must be a directory or a .zip file, got: {path!r}")

        # Case 2: .zip archive materialization
        if is_zip:
            arcname = self.arcname
            if cache_mode == "disable":
                download_stream_to_zip(
                    self.url or "",
                    path_p,
                    arcname,
                    chunk_size=cs,
                    overwrite=overwrite,
                    timeout=cfg.timeout,
                    insecure_ssl=cfg.insecure_ssl,
                    user_agent=cfg.user_agent,
                    attempts=cfg.max_retries,
                )
                return None
            cached = self._ensure_cached(mode=cache_mode, chunk_size=cs)
            write_cached_file_to_zip(cached, path_p, arcname, overwrite=overwrite)
            return None

        # Case 3: directory materialization
        dest_filename = Path(self.description.url_path()).name
        dest_path = path_p / dest_filename
        if cache_mode == "disable":
            download_to_file(
                self.url or "",
                dest_path,
                chunk_size=cs,
                timeout=cfg.timeout,
                insecure_ssl=cfg.insecure_ssl,
                user_agent=cfg.user_agent,
                attempts=cfg.max_retries,
            )
            self.filepath = str(dest_path)
            return str(dest_path)

        cached = self._ensure_cached(mode=cache_mode, chunk_size=cs)
        if dest_path.exists():
            if not overwrite:
                self.filepath = str(dest_path)
                return str(dest_path)
        _hardlink_or_copy(cached, dest_path)
        self.filepath = str(dest_path)
        return str(dest_path)

    # ------------------------------------------------------------------

    def load(self, *, force: bool = False) -> object:
        """Load the resource and cache the result on the instance.

        Returns:
          The loaded object (typically a pandas.DataFrame for tabular resources,
          or a :class:`BigWigFile` for BigWig resources).

        Raises:
          FileNotFoundError: If bytes are unexpectedly missing after download.
          LoadError: If the resource type is unsupported for loading.
        """
        if not force and self._cached_data is not None:
            return self._cached_data

        # BigWig: open with pyBigWig and return a managed reader.
        if getattr(self.description, "resource_type", None) == "bigwig_files":
            self.download(path=None, cache_mode="enable")
            cached = self._cached_path()
            if not cached.exists():
                raise FileNotFoundError(str(cached))
            reader = BigWigFile(cached)
            reader._ensure_open()
            self._cached_data = reader
            return reader

        # Ensure bytes are present locally (cache).
        self.download(path=None, cache_mode="enable")
        cached = self._cached_path()
        if not cached.exists():
            raise FileNotFoundError(str(cached))

        name = cached.name.lower()
        if name.endswith(".tsv") or name.endswith(".tsv.gz") or name.endswith(".md.gz"):
            if pd is not None:
                obj = pd.read_table(cached, compression="infer")
                self._cached_data = obj
                return obj

            # Fallback without pandas (legacy parity).
            import csv
            import gzip
            import io

            if name.endswith(".gz"):
                fh = io.TextIOWrapper(gzip.open(cached, "rb"), encoding="utf-8")
            else:
                fh = open(cached, "rt", encoding="utf-8")
            with fh as tsv:
                reader = csv.DictReader(tsv, delimiter="\t")
                obj = list(reader)
                self._cached_data = obj
                return obj

        raise LoadError(
            "Unsupported load() for resource type "
            f"{getattr(self.description, 'resource_type', None)!r} "
            f"with file {cached.name!r}"
        )

    # ------------------------------------------------------------------

    def is_loaded(self) -> bool:
        """Return whether this resource has an in-memory loaded object."""
        return self._cached_data is not None

    def get_loaded(self) -> object | None:
        """Return the in-memory loaded object without triggering I/O."""
        return self._cached_data

    def clear_loaded(self) -> None:
        """Clear only the in-memory loaded object (does not touch on-disk cache)."""
        obj = self._cached_data
        try:
            if isinstance(obj, BigWigFile):
                obj.close()
        finally:
            self._cached_data = None

    def __repr__(self) -> str:  # pragma: no cover
        cls = self.__class__.__name__
        return f"{cls}(url={self.url!r}, arcname={self.arcname!r}, filepath={self.filepath!r})"
