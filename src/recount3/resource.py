"""Resource orchestration: URL resolution, caching, downloading, and loading.

This module defines :class:`R3Resource`, the central class of the recount3
package. Every downloadable file (count matrices, metadata tables, annotation
GTFs, BigWig coverage files, junction counts) is represented as an
``R3Resource``. The class manages the full lifecycle of a resource:

1. Description -> URL: a :class:`~recount3._descriptions.R3ResourceDescription`
   provides the structured parameters (organism, project, genomic unit, …) that
   are used to construct the deterministic duffel-mirror URL.
2. URL -> cache: :meth:`~R3Resource.download` fetches the file over HTTP
   and stores it in a persistent on-disk cache keyed by URL hash.
3. Cache -> materialization: the cached file can be hard-linked or copied
   to a user-supplied directory, or appended to a ZIP archive.
4. Cache -> in-memory object: :meth:`~R3Resource.load` parses the cached
   file into an appropriate Python object (see Notes below).

Typical usage example::

    from recount3 import R3Resource, R3GeneOrExonCounts

    desc = R3GeneOrExonCounts(
        organism="human",
        data_source="sra",
        genomic_unit="gene",
        project="SRP009615",
        annotation_extension="G026",
    )
    res = R3Resource(desc)

    # Cache the file (no local copy):
    res.download(path=None, cache_mode="enable")

    # Or copy into a directory:
    dest = res.download(path="/data/recount3")

    # Parse into a DataFrame (for count resources):
    df = res.load()

Note:
    Downloads are protected by a module-level :class:`threading.Lock`, so
    multiple threads can safely call :meth:`~R3Resource.download` on
    resources sharing a common cache path without corrupting the cache.

Note:
    :meth:`~R3Resource.load` returns different types depending on the
    resource type:

    * Gene/exon count resources -> :class:`pandas.DataFrame` (features x samples).
    * Junction MM resources -> :class:`scipy.sparse.csr_matrix`.
    * Junction ID/RR resources -> :class:`pandas.DataFrame`.
    * BigWig resources -> :class:`~recount3._bigwig.BigWigFile`.
"""

from __future__ import annotations

import dataclasses
import gzip
import threading
import urllib.parse
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import scipy.io
import scipy.sparse

from recount3._bigwig import BigWigFile
from recount3._descriptions import R3ResourceDescription
from recount3._utils import (
    _cache_path,
    _derive_junction_sidecar_url,
    _ensure_dir,
    _hardlink_or_copy,
    download_stream_to_zip,
    download_to_file,
    write_cached_file_to_zip,
)
from recount3.config import Config, default_config
from recount3.errors import LoadError
from recount3.types import CacheMode

_FILE_LOCK = threading.Lock()


def _ensure_cached_url(
    *,
    url: str,
    cache_root: Path,
    cfg: Config,
    chunk_size: int,
) -> Path:
    """Ensure a given URL is downloaded and cached locally.

    Checks the specified cache root directory for the presence of the file
    corresponding to the provided URL. If the file is missing, it acquires a
    thread lock and downloads the file in chunks.

    Args:
        url: The full, absolute remote URL of the file to download.
        cache_root: The file system directory path where cached data should
            be persistently stored.
        cfg: The global configuration object, which provides network parameters
            such as connection timeouts, maximum retries, SSL verification
            preferences, and the user agent string.
        chunk_size: The specific number of bytes to read and write per chunk
            during the download process.

    Returns:
        The absolute path pointing to the locally cached file.
    """
    cache_path = _cache_path(url, cache_root)
    _ensure_dir(cache_path.parent)
    if cache_path.exists():
        return cache_path

    with _FILE_LOCK:
        if cache_path.exists():
            return cache_path
        download_to_file(
            url,
            cache_path,
            chunk_size=chunk_size,
            timeout=cfg.timeout,
            insecure_ssl=cfg.insecure_ssl,
            user_agent=cfg.user_agent,
            attempts=cfg.max_retries,
        )
    return cache_path


def _read_id_rail_ids(id_path: Path) -> list[str]:
    """Read a junction ``.ID`` or ``.ID.gz`` file and return rail IDs in order.

    The file is a tab-separated table with a header row. The function looks for
    a column named ``rail_id`` (case-insensitive); if absent, it falls back to
    the first column. A second parse attempt with ``sep=None`` (Python engine)
    is made if the initial tab-delimited parse fails.

    Args:
        id_path: Path to the junction ID file (``.ID`` or ``.ID.gz``).

    Returns:
        An ordered list of rail ID strings, one per row.

    Raises:
        LoadError: If the file parses to an empty DataFrame or if the resulting
            list of IDs is empty.
    """
    try:
        df = pd.read_csv(
            id_path,
            compression="infer",
            sep="\t",
            header=0,
            engine="c",
            low_memory=False,
        )
    except Exception:
        df = pd.read_csv(
            id_path,
            compression="infer",
            sep=None,
            header=0,
            engine="python",
        )

    if df.empty:
        raise LoadError(f"Junction ID file {id_path.name!r} parsed empty.")

    cols_lower = [str(c).strip().lower() for c in df.columns]
    rail_col = None
    if "rail_id" in cols_lower:
        rail_col = df.columns[cols_lower.index("rail_id")]
    else:
        rail_col = df.columns[0]

    rail_ids = df[rail_col].astype(str).tolist()
    if not rail_ids:
        raise LoadError(f"Junction ID file {id_path.name!r} has no rail IDs.")
    return rail_ids


def _read_mm_matrix(mm_path: Path) -> scipy.sparse.csr_matrix:
    """Read a MatrixMarket file into a SciPy CSR sparse matrix.

    Args:
        mm_path: The local file system path pointing to the MatrixMarket
            file (typically ending in .MM or .MM.gz).

    Returns:
        A 2-dimensional sparse matrix accurately representing the data parsed
        from the file, formatted as a SciPy CSR matrix.

    Raises:
        LoadError: If the file cannot be opened, fails to parse as a valid
            MatrixMarket format, or if the resulting parsed matrix does not
            contain exactly two dimensions.
    """
    try:
        if mm_path.suffix.lower() == ".gz":
            with gzip.open(mm_path, "rb") as fh:
                mat = scipy.io.mmread(fh)
        else:
            mat = scipy.io.mmread(str(mm_path))
    except Exception as exc:
        raise LoadError(f"Failed to read MatrixMarket from {mm_path.name!r}.") from exc

    if getattr(mat, "ndim", None) != 2:
        raise LoadError(f"MatrixMarket {mm_path.name!r} is not 2-dimensional.")
    return scipy.sparse.csr_matrix(mat)


@dataclass(slots=True)
class R3Resource:
    """Resource that manages URL resolution, caching, materialization, and loading.

    Attributes:
        description: An instance of `R3ResourceDescription` that specifies the
            metadata and the hierarchical path used to correctly locate and
            define the resource.
        url: The full, absolute network URL pointing to the remote resource. If
            not explicitly provided during initialization, it is derived by
            joining the configured base URL with the description's relative URL.
        filepath: An optional string representing the absolute local path where
            the resource was successfully materialized (either copied or linked).
        config: An optional `Config` instance dictating strict network and cache
            behaviors. If omitted, the global default configuration is dynamically
            applied.
    """

    description: R3ResourceDescription
    url: str | None = None
    filepath: str | None = None
    config: Config | None = None

    _cached_data: object | None = dataclasses.field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        """Initialize derived fields after dataclass instantiation."""
        cfg = self.config or default_config()
        if self.url is None:
            self.url = urllib.parse.urljoin(cfg.base_url, self.description.url_path())

    @property
    def arcname(self) -> str:
        """Deterministic ZIP arcname derived from the URL path."""
        return self.description.url_path().lstrip("/")

    def _cache_root(self) -> Path:
        """Resolve the root directory path used for local file caching.

        Returns:
            The absolute directory path where safely downloaded files are
            globally cached, determined entirely by the active configuration.
        """
        cfg = self.config or default_config()
        return Path(cfg.cache_dir)

    def _cached_path(self) -> Path:
        """Compute the local cache file path for this specific resource.

        Returns:
            The fully resolved absolute path indicating exactly where this
            resource's data persistently resides within the local cache hierarchy.
        """
        return _cache_path(self.url or "", self._cache_root())

    def _ensure_cached(self, *, mode: CacheMode, chunk_size: int) -> Path:
        """Ensure the resource data is available in the local cache.

        Args:
            mode: The primary caching strategy to employ. Must be 'enable' to
                use pre-existing cached files (downloading only if absent),
                or 'update' to force a fresh network download, overriding 
                any previously existing files.
            chunk_size: The fixed size in bytes utilized for reading and
            writing active network streams.

        Returns:
            The absolute path object directing to the successfully locally
            cached file.

        Raises:
            ValueError: If the caching mode is strictly set to 'disable' (as
                local caching operations are fundamentally invalid in this
                mode), or if the universally provided mode string is completely
                unrecognized.
        """
        cfg = self.config or default_config()
        cache_path = self._cached_path()
        if mode == "disable":
            raise ValueError("Cache is not used in 'disable' mode")

        if mode == "enable":
            return _ensure_cached_url(
                url=self.url or "",
                cache_root=self._cache_root(),
                cfg=cfg,
                chunk_size=chunk_size,
            )

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

    def download(
        self,
        path: str | None = None,
        *,
        cache_mode: CacheMode = "enable",
        overwrite: bool = False,
        chunk_size: int | None = None,
    ) -> str | None:
        """Ensure resource availability and optionally materialize it.

        Transitions the remote resource to the local system. Caches the file,
        writes it to a specific directory, or appends it to a ZIP archive
        depending on the arguments provided.

        Args:
            path: Target destination. If None, performs a cache-only download.
                If a directory path, links or copies the file there. If a
                '.zip' path, injects the file into the archive using `arcname`.
            cache_mode: Caching behavior. 'enable' uses existing cache, 'disable'
                streams directly to `path` without caching, 'update' forces a
                cache refresh before materialization.
            overwrite: If True, replaces existing files at the destination.
            chunk_size: Byte size for streaming operations. Defaults to the
                configured chunk size.

        Returns:
            The final file path if materialized to a directory. None if performing
            a cache-only download or appending to a ZIP archive.

        Raises:
            ValueError: Combinations are invalid (e.g., path=None with
                cache_mode='disable') or path has an unsupported format.

        Examples:
            Cache the file without copying it anywhere::

                res.download(path=None, cache_mode="enable")

            Copy the cached file into a local directory::

                dest = res.download(path="/data/recount3")

            Append the file to a ZIP archive::

                res.download(path="/data/recount3.zip")

            Force a cache refresh before copying::

                dest = res.download(path="/data/recount3", cache_mode="update")
        """
        cfg = self.config or default_config()
        if cfg.cache_disabled:
            cache_mode = "disable"
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

    def load(self, *, force: bool = False) -> object:
        """Parse the resource into an appropriate in-memory Python object.

        Downloads and caches the resource if missing. Uses the resource description
        to determine the parsing strategy (e.g., BigWig, tabular counts, junctions).
        Caches the resulting object internally to prevent redundant disk I/O.

        Args:
            force: If True, bypasses the in-memory object cache and re-parses
                data directly from disk.

        Returns:
            The parsed object. Tabular and junction counts return a
            `pandas.DataFrame`. BigWig files return a
            `recount3._bigwig.BigWigFile` instance.

        Raises:
            FileNotFoundError: The file is missing from the cache post-download.
            LoadError: Parsing fails, matrix shapes mismatch, or the resource
                type is currently unsupported.

        Examples:
            Load a gene/exon count matrix as a DataFrame::

                counts_df = gene_count_res.load()  # -> pd.DataFrame

            Load a BigWig coverage file (close it when done)::

                bw = bigwig_res.load()  # -> BigWigFile
                vals = bw.values("chr1", 0, 1_000_000)
                bw.close()

            Re-parse from disk, bypassing the in-memory cache::

                counts_df = res.load(force=True)
        """
        if not force and self._cached_data is not None:
            return self._cached_data

        if getattr(self.description, "resource_type", None) == "bigwig_files":
            self.download(path=None, cache_mode="enable")
            cached = self._cached_path()
            if not cached.exists():
                raise FileNotFoundError(str(cached))
            reader = BigWigFile(cached)
            reader._ensure_open()
            self._cached_data = reader
            return reader

        rtype = getattr(self.description, "resource_type", None)
        if rtype == "count_files_gene_or_exon":
            self.download(path=None, cache_mode="enable")
            cached = self._cached_path()
            if not cached.exists():
                raise FileNotFoundError(str(cached))
            try:
                df = pd.read_csv(
                    cached,
                    compression="infer",
                    sep="\t",
                    header=0,
                    comment="#",
                    index_col=None,
                    engine="c",
                    low_memory=False,
                )
            except Exception:
                df = pd.read_csv(
                    cached,
                    compression="infer",
                    sep=None,
                    header=0,
                    comment="#",
                    index_col=None,
                    engine="python",
                )
            if df.empty or df.shape[1] < 2:
                raise LoadError(
                    f"Parsed an empty or 1-column matrix from {cached.name!r}."
                )
            cols_lower = [str(c).strip().lower() for c in df.columns]
            try_index = None
            for candidate in ("gene_id", "exon_id", "feature_id"):
                if candidate in cols_lower:
                    try_index = df.columns[cols_lower.index(candidate)]
                    break
            if try_index is None:
                try_index = df.columns[0]
            df = df.set_index(try_index)
            df.columns = [str(c) for c in df.columns]
            df.index = df.index.astype(str)
            self._cached_data = df
            return df

        if rtype == "count_files_junctions":
            jxn_ext = str(getattr(self.description, "junction_extension", "MM")).upper()

            cfg = self.config or default_config()
            self.download(path=None, cache_mode="enable")
            mm_cached = self._cached_path()
            if not mm_cached.exists():
                raise FileNotFoundError(str(mm_cached))

            if jxn_ext == "MM":
                try:
                    id_url = _derive_junction_sidecar_url(self.url or "", "ID")
                except Exception as exc:
                    raise LoadError(
                        f"Cannot derive junction ID URL from {self.url!r}."
                    ) from exc

                id_cached = _ensure_cached_url(
                    url=id_url,
                    cache_root=self._cache_root(),
                    cfg=cfg,
                    chunk_size=cfg.chunk_size,
                )

                rail_ids = _read_id_rail_ids(id_cached)
                mat = _read_mm_matrix(mm_cached)

                mat_shape = mat.shape
                if mat_shape is None:
                    raise LoadError(f"MatrixMarket {mm_cached.name!r} has undefined shape.")
                if mat_shape[1] != len(rail_ids):
                    raise LoadError(
                        "Junction MM column count does not match ID length: "
                        f"{mat_shape[1]} != {len(rail_ids)} "
                        f"({mm_cached.name!r} vs {id_cached.name!r})."
                    )

                df = pd.DataFrame.sparse.from_spmatrix(mat)
                df.columns = [str(x) for x in rail_ids]
                df.index = [str(i) for i in range(df.shape[0])]

                self._cached_data = df
                return df

            if jxn_ext in {"ID", "RR"}:
                obj = pd.read_table(mm_cached, compression="infer")
                self._cached_data = obj
                return obj

            raise LoadError(
                f"Unsupported junction_extension {jxn_ext!r} for {mm_cached.name!r}."
            )

        self.download(path=None, cache_mode="enable")
        cached = self._cached_path()
        if not cached.exists():
            raise FileNotFoundError(str(cached))

        name = cached.name.lower()
        if name.endswith(".tsv") or name.endswith(".tsv.gz") or name.endswith(".md.gz"):
            obj = pd.read_table(cached, compression="infer")
            self._cached_data = obj
            return obj

        raise LoadError(
            "Unsupported load() for resource type "
            f"{getattr(self.description, 'resource_type', None)!r} "
            f"with file {cached.name!r}"
        )

    def is_loaded(self) -> bool:
        """Check if the resource currently holds a parsed in-memory object.

        Returns:
            True if an object is cached in memory, False otherwise.
        """
        return self._cached_data is not None

    def get_loaded(self) -> object | None:
        """Retrieve the parsed in-memory object without triggering disk I/O.

        Returns:
            The loaded object if present, otherwise None.
        """
        return self._cached_data

    def clear_loaded(self) -> None:
        """Evict the in-memory object cache and close file handles if applicable.

        Does not delete or modify the on-disk file cache.
        """
        obj = self._cached_data
        try:
            if isinstance(obj, BigWigFile):
                obj.close()
        finally:
            self._cached_data = None

    def __repr__(self) -> str:
        """Return a string representation of the instance."""
        cls = self.__class__.__name__
        return f"{cls}(url={self.url!r}, arcname={self.arcname!r}, filepath={self.filepath!r})"
