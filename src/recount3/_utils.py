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
"""Internal utilities for HTTP, filesystem, and caching operations.

The module is organized into four sections:

1. Caching Utilities: URL-based cache key generation and path management
2. Filesystem Utilities: Directory management and atomic file operations
3. HTTP Utilities: Network requests with retries, streaming, and ZIP handling
4. Optional Dependency Management: lazy import helpers for BiocPy packages
   (``biocframe``, ``summarizedexperiment``, ``genomicranges``), ``pyBigWig``,
   ``anndata``, and Parquet engines; raises a standardized :exc:`ImportError`
   with an install command when a required optional package is missing.

Note:
    This module is considered internal implementation detail and may change
    without notice. Users should prefer the public API in the main package.

Examples:
    Typical usage is through the public API, but direct usage might look like::

        from pathlib import Path
        from recount3._utils import _cache_path

        path = _cache_path("https://example.com/data.tsv", Path("/tmp/cache"))

Attributes:
    _PATH_LOCKS_GUARD: Threading lock to safely interact with the weakref
        dictionary.
    _PATH_LOCKS: A weak reference dictionary mapping canonical paths to locks.
        This synchronizes cache and archive writes per file, preventing
        race conditions without leaking memory for inactive cache paths.
"""

from __future__ import annotations

import contextlib
import datetime
import errno
import hashlib
import http.client
import logging
import os
import re
import shutil
import socket
import ssl
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import functools
import importlib
import types
import weakref
from collections.abc import Iterator, Mapping
from pathlib import Path
from typing import BinaryIO, Any, cast, TYPE_CHECKING

import pandas as pd

from recount3 import errors

if TYPE_CHECKING:  # pragma: no cover
    import biocframe  # type: ignore[import-not-found]
    import genomicranges  # type: ignore[import-not-found]
    import summarizedexperiment  # type: ignore[import-not-found]

_PATH_LOCKS_GUARD = threading.Lock()
_PATH_LOCKS: weakref.WeakValueDictionary[str, _WeakRefLock] = (
    weakref.WeakValueDictionary()
)


class _WeakRefLock:
    """A threading.Lock wrapper that supports weak references.

    :class:`threading.Lock` objects cannot be stored in a
    :class:`weakref.WeakValueDictionary` because they do not support weak
    references. This class wraps a lock instance, making the wrapper itself
    weakly referenceable while proxying the context-manager protocol
    (``__enter__`` / ``__exit__``) to the inner lock.
    """

    def __init__(self):
        self._lock = threading.Lock()

    def __enter__(self):
        return self._lock.__enter__()

    def __exit__(self, *args):
        return self._lock.__exit__(*args)


def _strip_extended_prefix(text: str) -> str:
    """Return a resolved path without a Windows extended-length prefix.

    :meth:`pathlib.Path.resolve` keeps the ``\\\\?\\`` prefix whenever it
    cannot confirm that the plain spelling names the same file, which is
    the case while another thread holds that file open. One destination
    then has two spellings, and callers reaching it by different spellings
    would key separate locks. Paths without the prefix, including every
    POSIX path, are returned unchanged.

    Args:
        text: A resolved filesystem path.

    Returns:
        The path with any extended-length prefix removed.
    """
    prefix = "\\\\?\\"
    unc_prefix = prefix + "UNC\\"
    if text[: len(unc_prefix)].upper() == unc_prefix:
        return "\\\\" + text[len(unc_prefix) :]
    if text.startswith(prefix):
        return text[len(prefix) :]
    return text


def _path_lock_for_path(path: Path) -> _WeakRefLock:
    """Return a weakly retained lock for a canonical filesystem destination.

    The key is normalized so that spellings Windows treats as one file share
    a single lock during registry initialization, downloads, and ZIP writes.

    Args:
        path: Target cache payload, registry directory, or archive path.

    Returns:
        A lock shared by all writers to this canonical destination.
    """
    resolved = _strip_extended_prefix(str(path.expanduser().resolve()))
    key = os.path.normcase(resolved)
    with _PATH_LOCKS_GUARD:
        lock = _PATH_LOCKS.get(key)
        if lock is None:
            lock = _WeakRefLock()
            _PATH_LOCKS[key] = lock
        return lock


def _zip_lock_for_path(zip_path: Path) -> _WeakRefLock:
    """Return the process-local lock protecting writes to this archive."""
    return _path_lock_for_path(zip_path)


@contextlib.contextmanager
def _biocfilecache(root: Path):
    """Open the optional shared R/Python registry under a process-local lock.

    Network transfers must take place outside this context. Separate processes
    and R sessions require external coordination when mutating the same cache.
    """
    module = import_optional_module("pybiocfilecache")
    root = root.expanduser().resolve()
    with _path_lock_for_path(root / "BiocFileCache.sqlite"):
        with module.BiocFileCache(root) as cache:
            yield cache


def _biocfilecache_rows(cache, url: str | None = None) -> list[dict]:
    """Read records without the upstream get() wait for missing payloads."""
    model = import_optional_module("pybiocfilecache.models").Resource
    with cache.get_session() as session:
        query = session.query(model)
        if url is not None:
            query = query.filter(model.rname == url)
        return [row.to_dict() for row in query.all()]


def _biocfilecache_path(url: str, root: Path) -> Path:
    """Resolve an R- or Python-registered URL, or its native cache destination.

    Relative database paths are resolved against the cache root, including R's
    web records. Missing payloads retain their registered destination for
    atomic repair. Duplicate URL names are ambiguous and raise ValueError.
    """
    with _biocfilecache(root) as cache:
        rows = _biocfilecache_rows(cache, url)
        if len(rows) > 1:
            raise ValueError(f"Multiple BiocFileCache entries for URL: {url}")
        if rows:
            return (cache.config.cache_dir / rows[0]["rpath"]).resolve()
        return _cache_path(url, cache.config.cache_dir)


def _register_biocfilecache(
    root: Path, url: str, path: Path, *, transferred: bool = False
) -> None:
    """Register a complete payload while preserving R resource paths/types.

    Callers hold the payload lock. Local records use the optional package's
    checksum implementation. A refreshed web record loses old HTTP validators,
    because recount3's unconditional transfer does not collect new validators.
    The registry and payload are not a single cross-process transaction.
    """
    model = import_optional_module("pybiocfilecache.models").Resource
    utilities = import_optional_module("pybiocfilecache.utils")
    with _biocfilecache(root) as cache:
        rows = _biocfilecache_rows(cache, url)
        if not rows:
            record = cache.add(
                rname=url, fpath=path.resolve(), rtype="local", action="asis"
            )
            expected_rid = f"BFC{record['id']}"
            if record["rid"] != expected_rid:
                with cache.get_session() as session:
                    session.query(model).filter(
                        model.id == record["id"]
                    ).update({"rid": expected_rid}, synchronize_session=False)
            return
        record = rows[0]
        values = {
            "access_time": datetime.datetime.now(),
            "last_modified_time": record["last_modified_time"],
        }
        if record["rtype"] == "web":
            if transferred:
                values.update(etag=None, last_modified_time=None, expires=None)
        else:
            values["etag"] = utilities.calculate_file_hash(
                path, cache.config.hash_algorithm
            )
        with cache.get_session() as session:
            session.query(model).filter(model.id == record["id"]).update(
                values, synchronize_session=False
            )


def _biocfilecache_paths(root: Path) -> list[Path]:
    """List registered paths, including payloads removed outside the cache."""
    if not (root / "BiocFileCache.sqlite").exists():
        return []
    with _biocfilecache(root) as cache:
        return [
            (cache.config.cache_dir / row["rpath"]).resolve()
            for row in _biocfilecache_rows(cache)
        ]


def _cache_owned(path: Path, root: Path) -> bool:
    """Report whether the cache owns this payload and may delete it.

    R registers external files with ``action="asis"`` and leaves them in place
    on ``bfcremove``; only payloads inside the cache directory are owned.
    """
    path = path.resolve()
    root = root.resolve()
    return path == root or root in path.parents


def _remove_biocfilecache_files(root: Path, paths: list[Path]) -> None:
    """Remove selected payloads and registry rows while the cache is idle.

    Relative R paths are resolved against the cache root before unlinking.
    Native, unregistered payloads are also removed. Database files are retained.
    Registered files outside the cache directory keep R's ``asis`` semantics:
    their rows are dropped, but the external payloads are never deleted.
    """
    model = import_optional_module("pybiocfilecache.models").Resource
    with _biocfilecache(root) as cache:
        cache_dir = cache.config.cache_dir
        selected = {path.resolve() for path in paths}
        with cache.get_session() as session:
            for row in session.query(model).all():
                path = (cache_dir / row.rpath).resolve()
                if path in selected:
                    if _cache_owned(path, cache_dir):
                        path.unlink(missing_ok=True)
                    session.delete(row)
            for path in paths:
                if _cache_owned(path, cache_dir):
                    path.unlink(missing_ok=True)


def _cache_internal(path: Path, root: Path) -> bool:
    """Identify SQLite databases, journals, and R's database lock file."""
    return path.parent.resolve() == root.resolve() and (
        path.name == "BiocFileCache.sqlite"
        or path.name.startswith("BiocFileCache.sqlite-")
        or path.name == "BiocFileCache.sqlite.LOCK"
    )


def _sha256(text: str) -> str:
    """Return the hex SHA256 digest of input text.

    Args:
        text: Input string to hash.

    Returns:
        64-character hexadecimal SHA256 digest.
    """
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _cache_key_for_url(url: str) -> str:
    """Generate stable cache key derived from full URL.

    Creates a deterministic cache filename using a short SHA256 prefix
    combined with the URL path's basename. This ensures unique but
    recognizable cache entries.

    Args:
        url: Full resource URL to generate key for.

    Returns:
        Cache key in format: `{sha256_prefix}__{url_basename}`.

    Example:
        >>> _cache_key_for_url("http://example.com/path/to/file.tsv.gz")
        'a1b2c3d4e5f67890__file.tsv.gz'
    """
    parsed = urllib.parse.urlparse(url)
    key = parsed.path.lstrip("/")
    digest = _sha256(url)[:16]  # First 16 chars for reasonable uniqueness
    return f"{digest}__{Path(key).name}"


def _cache_path(url: str, cache_root: str | Path) -> Path:
    """Return full filesystem path for cached URL content.

    Args:
        url: Resource URL to locate in cache.
        cache_root: Base cache directory path.

    Returns:
        Absolute Path where URL content should be cached.

    Raises:
        NotADirectoryError: If cache_root exists but is not a directory.
    """
    _ensure_dir(cache_root)
    return Path(cache_root) / _cache_key_for_url(url)


def _ensure_dir(path: str | Path) -> None:
    """Ensure directory exists, creating parents as needed.

    Args:
        path: Directory path to create.

    Raises:
        NotADirectoryError: If path exists as a non-directory.
        OSError: If directory creation fails due to permissions or other system
          error.
    """
    p = Path(path)
    try:
        p.mkdir(parents=True, exist_ok=True)
    except FileExistsError as exc:
        raise NotADirectoryError(f"Exists but is not a directory: {p}") from exc


def _hardlink_or_copy(src: Path, dst: Path) -> None:
    """Materialize file by hardlink, falling back to copy on errors.

    Attempts to create a hardlink for efficiency. Falls back to copy operation
    on cross-device, permission, or too-many-links errors.

    Args:
        src: Source file path.
        dst: Destination file path.

    Raises:
        OSError: On unexpected filesystem errors beyond the handled cases.
        FileNotFoundError: If source file doesn't exist.
    """
    _ensure_dir(dst.parent)
    tmp_dst = dst.parent / (
        f".{dst.name}.{os.getpid()}_"
        f"{threading.get_ident()}_{time.time_ns()}.tmp"
    )

    try:
        try:
            os.link(src, tmp_dst)
        except OSError as e:
            if e.errno in (
                errno.EXDEV,
                errno.EPERM,
                errno.EACCES,
                errno.EMLINK,
            ):
                shutil.copy2(src, tmp_dst)
            else:
                raise
        _atomic_replace(tmp_dst, dst)
    finally:
        if tmp_dst.exists():
            try:
                tmp_dst.unlink()
            except OSError:
                pass


def _atomic_replace(src_tmp: Path, final_path: Path) -> None:
    """Atomically replace final path with temporary file.

    Uses filesystem atomic replace operation to ensure the destination
    file either exists completely or not at all, preventing partial
    file states.

    Args:
        src_tmp: Temporary file path containing new content.
        final_path: Final destination path for the content.

    Raises:
        OSError: If the replace operation fails.
    """
    _ensure_dir(final_path.parent)
    os.replace(src_tmp, final_path)


def _ssl_insecure_context() -> ssl.SSLContext:
    """Return SSL context with verification disabled (not recommended).

    Warning:
        This disables certificate verification and should only be used
        for testing or with trusted networks.

    Returns:
        SSL context with hostname verification and certificate checking
        disabled.
    """
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def http_open(
    url: str,
    *,
    timeout: float,
    headers: dict[str, str] | None,
    insecure_ssl: bool,
    user_agent: str,
) -> BinaryIO:
    """Open URL with urllib, ensuring all headers are valid strings.

    Args:
        url: Absolute URL to open.
        timeout: Socket timeout in seconds.
        headers: Optional extra headers to send. None values are dropped.
        insecure_ssl: If True, disable TLS certificate verification. This adds
            an unverified ``HTTPSHandler``, which urllib uses only for
            ``https://`` URLs; it has no effect on plain ``http://`` URLs
            (e.g. the default Duffel mirror).
        user_agent: User-Agent header override.

    Returns:
        Binary HTTP response object.

    Raises:
        URLError: For network-level errors.
        HTTPError: For HTTP protocol errors.
        socket.timeout: For connection timeouts.
        ValueError: If URL is empty.
    """
    if not url:
        raise ValueError("Empty URL passed to http_open().")

    merged: dict[str, str] = {"User-Agent": str(user_agent)}
    if headers:
        for k, v in headers.items():
            if v is None:
                continue
            merged[str(k)] = str(v)

    opener = urllib.request.build_opener()
    if insecure_ssl:
        opener.add_handler(
            urllib.request.HTTPSHandler(context=_ssl_insecure_context())
        )

    req = urllib.request.Request(url, headers=merged)
    return opener.open(req, timeout=timeout)


def with_retries(func, *, attempts: int, base_sleep: float = 0.5):
    """Execute function with simple exponential backoff on transient errors.

    Implements retry logic for network operations that may fail transiently.
    Sleeps between attempts follow exponential backoff pattern.

    Args:
        func: Callable to execute with retries.
        attempts: Maximum number of execution attempts.
        base_sleep: Base sleep time in seconds for backoff calculation.

    Returns:
        Return value of the successful function execution.

    Raises:
        Exception: The last exception encountered if all attempts fail.
    """
    last_exc: BaseException | None = None
    for i in range(max(1, attempts)):
        try:
            return func()
        except (
            urllib.error.URLError,
            socket.timeout,
            TimeoutError,
            ConnectionError,
            http.client.HTTPException,
            ssl.SSLError,
        ) as exc:
            last_exc = exc
            if i < attempts - 1:
                time.sleep(base_sleep * (2**i))
    if last_exc:  # pragma: no branch
        raise last_exc


def _stream_copy(src_fh, dst_fh, *, chunk_size: int) -> int:
    """Copy bytes in chunks from source to destination filehandle.

    Args:
        src_fh: Source file-like object supporting read().
        dst_fh: Destination file-like object supporting write().
        chunk_size: Number of bytes to read/write per chunk.

    Returns:
        Total number of bytes copied.
    """
    total = 0
    while True:
        chunk = src_fh.read(chunk_size)
        if not chunk:
            break
        dst_fh.write(chunk)
        total += len(chunk)
    return total


def download_to_file(
    url: str,
    out_path: Path,
    *,
    chunk_size: int,
    timeout: float,
    insecure_ssl: bool,
    user_agent: str,
    attempts: int,
) -> None:
    """Download URL atomically into output path.

    Streams content directly to temporary file then atomically replaces
    destination to avoid partial files on failure. Implements retry logic
    for transient network errors.

    Args:
        url: Resource URL to download.
        out_path: Destination file path.
        chunk_size: Bytes to read/write per chunk during streaming.
        timeout: Network timeout in seconds.
        insecure_ssl: If True, disable TLS verification.
        user_agent: HTTP User-Agent header.
        attempts: Maximum retry attempts for transient errors.

    Raises:
        DownloadError: If download fails after all retry attempts.
        OSError: If filesystem operations fail.
    """

    def _do():
        with contextlib.closing(
            http_open(
                url,
                timeout=timeout,
                headers=None,
                insecure_ssl=insecure_ssl,
                user_agent=user_agent,
            )
        ) as resp:
            _ensure_dir(out_path.parent)

            tmp = out_path.parent / (
                f".{out_path.name}.{os.getpid()}_"
                f"{threading.get_ident()}_{time.time_ns()}.downloading"
            )

            try:
                with open(tmp, "wb") as fh:
                    _stream_copy(resp, fh, chunk_size=chunk_size)
                _atomic_replace(tmp, out_path)
            finally:
                if tmp.exists():
                    try:
                        tmp.unlink()
                    except OSError:
                        pass

    try:
        with_retries(_do, attempts=attempts)
    except Exception as exc:
        raise errors.DownloadError(
            f"Failed to download {url!r} -> {str(out_path)!r}"
        ) from exc


def _write_or_replace_in_zip(
    zip_path: Path,
    source_path: Path,
    arcname: str,
    overwrite: bool,
) -> None:
    """Write ``source_path`` into ``zip_path`` as ``arcname``.

    If ``zip_path`` already contains a member named ``arcname``:

    * When ``overwrite`` is ``False``, this function returns without changes.
    * When ``overwrite`` is ``True``, the ZIP is rewritten without the old
      member and then the new file is added.

    Rewriting is necessary because ZIP archives can contain multiple members
    with the same name; appending another member does not remove the existing
    one.

    Args:
        zip_path: Path to the destination ZIP archive.
        source_path: Local path to the file to add.
        arcname: Member name to use inside the ZIP.
        overwrite: If ``True``, replace an existing member named ``arcname``.
            If ``False``, this function is a no-op when that member already
            exists.

    Raises:
        FileNotFoundError: If ``source_path`` does not exist.
        DownloadError: If ``zip_path`` exists but is not a valid ZIP archive.
        zipfile.BadZipFile: If the ZIP file is malformed.
        OSError: For filesystem errors while writing the ZIP.
    """
    _ensure_dir(zip_path.parent)

    with _zip_lock_for_path(zip_path):
        if not zip_path.exists():
            with zipfile.ZipFile(
                zip_path, "w", compression=zipfile.ZIP_DEFLATED
            ) as zf:
                zf.write(source_path, arcname)
            return

        if not zipfile.is_zipfile(zip_path):
            raise errors.DownloadError(
                f"Destination {zip_path} exists but is not a valid ZIP."
            )

        member_exists = False
        with zipfile.ZipFile(zip_path, "r") as zf:
            try:
                zf.getinfo(arcname)
                member_exists = True
            except KeyError:
                pass

        if member_exists and not overwrite:
            return

        if member_exists:
            tmp_zip = zip_path.parent / (
                f".{zip_path.name}.{os.getpid()}_"
                f"{threading.get_ident()}_{time.time_ns()}.tmpzip"
            )
            try:
                with zipfile.ZipFile(zip_path, "r") as zf_in:
                    with zipfile.ZipFile(
                        tmp_zip, "w", compression=zipfile.ZIP_DEFLATED
                    ) as zf_out:
                        zf_out.comment = zf_in.comment

                        for info in zf_in.infolist():
                            if info.filename == arcname:
                                continue

                            out_info = zipfile.ZipInfo(
                                info.filename, date_time=info.date_time
                            )
                            out_info.compress_type = info.compress_type
                            out_info.comment = info.comment
                            out_info.extra = info.extra
                            out_info.create_system = info.create_system
                            out_info.create_version = info.create_version
                            out_info.extract_version = info.extract_version
                            out_info.flag_bits = info.flag_bits
                            out_info.internal_attr = info.internal_attr
                            out_info.external_attr = info.external_attr

                            with zf_in.open(info, "r") as f_in:
                                with zf_out.open(out_info, "w") as f_out:
                                    shutil.copyfileobj(f_in, f_out)

                        zf_out.write(source_path, arcname)

                _atomic_replace(tmp_zip, zip_path)

            finally:
                if tmp_zip.exists():
                    try:
                        tmp_zip.unlink()
                    except OSError:
                        pass
        else:
            with zipfile.ZipFile(
                zip_path, "a", compression=zipfile.ZIP_DEFLATED
            ) as zf:
                zf.write(source_path, arcname)


def download_stream_to_zip(
    url: str,
    zip_path: Path,
    arcname: str,
    *,
    chunk_size: int,
    overwrite: bool,
    timeout: float,
    insecure_ssl: bool,
    user_agent: str,
    attempts: int,
) -> None:
    """Download a URL to a temporary file, then write it into a ZIP archive.

    This uses a two-phase approach:

    1) Download the URL to a temporary file (with retries).
    2) Add that file to ``zip_path`` under ``arcname``.

    The temporary file avoids leaving a partially-written ZIP member if the
    network request fails mid-stream.

    If ``overwrite`` is ``True`` and the ZIP already contains ``arcname``, the
    ZIP is rewritten to avoid duplicate entries.

    Args:
        url: URL to download.
        zip_path: Path to the destination ZIP archive.
        arcname: Member name to use inside the ZIP.
        chunk_size: Number of bytes to read/write per chunk while downloading.
        overwrite: If ``True``, replace an existing member named ``arcname``.
            If ``False``, this function is a no-op when that member already
            exists.
        timeout: Network timeout in seconds.
        insecure_ssl: If ``True``, TLS certificate verification is disabled.
        user_agent: HTTP ``User-Agent`` header to send.
        attempts: Maximum retry attempts for transient network errors.

    Raises:
        DownloadError: If the download fails after all retry attempts, or if the
            destination exists but is not a valid ZIP file.
        zipfile.BadZipFile: If the ZIP file is malformed.
        OSError: For filesystem errors while writing the temporary file or ZIP.
    """
    _ensure_dir(zip_path.parent)

    if not overwrite and zip_path.exists():
        with _zip_lock_for_path(zip_path):
            if zip_path.exists() and zipfile.is_zipfile(zip_path):
                try:
                    with zipfile.ZipFile(zip_path, "r") as zf:
                        zf.getinfo(arcname)
                        return
                except KeyError:
                    pass

    tmp_path = zip_path.parent / (
        f".r3_dl_{os.getpid()}_{threading.get_ident()}_{time.time_ns()}.tmp"
    )
    try:

        def _do() -> None:
            with contextlib.closing(
                http_open(
                    url,
                    timeout=timeout,
                    headers=None,
                    insecure_ssl=insecure_ssl,
                    user_agent=user_agent,
                )
            ) as resp:
                with open(tmp_path, "wb") as f:
                    _stream_copy(resp, f, chunk_size=chunk_size)

        with_retries(_do, attempts=attempts)

        _write_or_replace_in_zip(zip_path, tmp_path, arcname, overwrite)

    finally:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass


def write_cached_file_to_zip(
    cached_file: Path,
    zip_path: Path,
    arcname: str,
    *,
    overwrite: bool,
) -> None:
    """Write an existing on-disk file into a ZIP archive.

    If the ZIP does not exist, it is created. If the ZIP already contains a
    member named ``arcname``:

    * When ``overwrite`` is ``False``, this function returns without changes.
    * When ``overwrite`` is ``True``, the ZIP is rewritten without the old
      member and then the new file is added to avoid duplicate entries.

    Args:
        cached_file: Path to a local file that already exists on disk.
        zip_path: Path to the destination ZIP archive.
        arcname: Member name to use inside the ZIP.
        overwrite: If ``True``, replace an existing member named ``arcname``.
            If ``False``, this function is a no-op when that member already
            exists.

    Raises:
        FileNotFoundError: If ``cached_file`` does not exist.
        DownloadError: If ``zip_path`` exists but is not a valid ZIP archive.
        zipfile.BadZipFile: If the ZIP file is malformed.
        OSError: For filesystem errors while reading ``cached_file`` or writing
            the ZIP.
    """
    _write_or_replace_in_zip(zip_path, cached_file, arcname, overwrite)


def _normalize_genomic_unit(genomic_unit: str) -> str:
    """Return a normalized genomic unit string and validate it.

    Args:
      genomic_unit: Requested feature level.

    Returns:
      Lowercase genomic unit string.

    Raises:
      ValueError: If the genomic unit is not one of ``"gene"``,
        ``"exon"``, or ``"junction"``.
    """
    gu = str(genomic_unit).strip().lower()
    valid = {"gene", "exon", "junction"}
    if gu not in valid:
        raise ValueError(
            f"Invalid genomic_unit {genomic_unit!r}; expected one of "
            f"{sorted(valid)!r}."
        )
    return gu


def _resolve_counts_assay_name(
    se_like: Any,
    *,
    preferred_assay_name: str = "raw_counts",
    fallback_assay_name: str = "counts",
) -> str:
    """Resolve the assay name that carries the recount3 coverage-sum matrix.

    Prefers ``preferred_assay_name`` when present; otherwise falls back to
    ``fallback_assay_name`` with a warning for backwards compatibility.
    """
    assay_names = getattr(se_like, "assay_names", None)
    if assay_names and preferred_assay_name in assay_names:
        return preferred_assay_name
    if assay_names and fallback_assay_name in assay_names:
        logging.warning(
            "Assay %r not found; falling back to legacy assay %r. "
            "Rebuild the object with assay_name=%r to silence this warning.",
            preferred_assay_name,
            fallback_assay_name,
            preferred_assay_name,
        )
        return fallback_assay_name
    raise ValueError(
        f"Object must contain a {preferred_assay_name!r} assay"
        + (
            f" (or legacy {fallback_assay_name!r})"
            if fallback_assay_name
            else ""
        )
        + "."
    )


def _coerce_col_data_to_pandas(sample_metadata_source: Any) -> pd.DataFrame:
    """Coerce sample metadata into a :class:`~pandas.DataFrame`.

    Args:
        sample_metadata_source: Either a BiocPy
          ``(Ranged)SummarizedExperiment``-like object with a
          `.col_data.to_pandas()` method, or a
          :class:`~pandas.DataFrame` already.

    Returns:
        A :class:`~pandas.DataFrame` of sample metadata.

    Raises:
        TypeError: If `sample_metadata_source` cannot be coerced to a
          :class:`~pandas.DataFrame`.
    """
    if isinstance(sample_metadata_source, pd.DataFrame):
        return sample_metadata_source

    if hasattr(sample_metadata_source, "col_data") and hasattr(
        sample_metadata_source.col_data, "to_pandas"
    ):
        return sample_metadata_source.col_data.to_pandas()

    raise TypeError(
        "Expected a pandas.DataFrame or a SummarizedExperiment-like "
        "object with `col_data.to_pandas()`."
    )


def _coerce_numeric_column(series: pd.Series, column_name: str) -> pd.Series:
    """Coerce a Series to numeric, erroring on non-numeric non-missing values.

    Args:
        series: Input Series.
        column_name: Name used for error messages.

    Returns:
        Float Series.

    Raises:
        ValueError: If non-missing values cannot be coerced to numeric.
    """
    cleaned = series.replace(r"^\s*$", pd.NA, regex=True)

    numeric = pd.to_numeric(cleaned, errors="coerce")

    invalid = cleaned.notna() & numeric.isna()
    if invalid.any():
        examples = cleaned[invalid].head(3).tolist()
        raise ValueError(
            f"Metadata column {column_name!r} contains non-numeric values "
            f"(examples: {examples!r})."
        )
    return numeric.astype(float)


_MAX_EXACT_FLOAT_INTEGER = 2**53


def canonical_identifier_series(values: pd.Series) -> pd.Series:
    """Render identifiers as text in a way that does not depend on dtype.

    recount3 ships each metadata table for a project as its own TSV, and
    each is parsed independently. A numeric identifier such as ``rail_id``
    can therefore land as ``int64`` in one table and as ``float64`` in
    another, because a single blank cell is enough to make pandas widen
    the column. Stringifying those two columns directly yields ``"123488"``
    and ``"123488.0"``, which no longer compare equal, so an inner join on
    the identifier silently collapses to zero rows. R never sees this
    because ``merge()`` compares the parsed numbers rather than their text
    form.

    Numeric columns whose present values are all whole numbers (and, for
    floats, small enough to be represented exactly) are therefore rendered
    through an integer, so the text depends only on the value. Everything
    else -- including non-integral numbers, text and booleans -- is
    stringified unchanged.

    Args:
        values: An identifier column parsed from a recount3 table.

    Returns:
        A ``string``-dtype Series with the same index as ``values``, with
        missing entries preserved as :data:`pandas.NA`.

    Examples:
        Two tables that parsed the same identifier differently still agree
        on the canonical text::

            >>> import pandas as pd
            >>> list(canonical_identifier_series(pd.Series([123488])))
            ['123488']
            >>> list(canonical_identifier_series(pd.Series([123488.0])))
            ['123488']
    """
    if pd.api.types.is_bool_dtype(values):
        return values.astype("string")
    if not pd.api.types.is_numeric_dtype(values):
        return values.astype("string")

    numeric = pd.to_numeric(values, errors="coerce")
    present = numeric.dropna()
    if not present.empty:
        if not bool((present == present.round()).all()):
            return values.astype("string")
        if pd.api.types.is_float_dtype(numeric) and not bool(
            (present.abs() <= _MAX_EXACT_FLOAT_INTEGER).all()
        ):
            return values.astype("string")

    canonical = pd.Series(pd.NA, index=values.index, dtype="string")
    if not present.empty:
        canonical.loc[present.index] = present.astype("int64").astype("string")
    return canonical


def _resolve_metadata_column(
    metadata_df: pd.DataFrame,
    column_name: str,
) -> pd.Series:
    """Resolve a metadata column name robustly.

    This mirrors the strictness of the recount3 R implementation (which expects
    exact column names), but also supports the Python-side convention where the
    namespace separator may be `__` instead of `.` for the first separator
    (e.g., `recount_qc.star.average_mapped_length` vs
    `recount_qc__star.average_mapped_length`).

    Args:
        metadata_df: Sample metadata.
        column_name: Column name to resolve.

    Returns:
        The resolved column as a pandas Series.

    Raises:
        ValueError: If the column cannot be found.
    """
    lower_to_actual = {str(c).lower(): c for c in metadata_df.columns}
    key = column_name.lower()
    if key in lower_to_actual:
        return metadata_df[lower_to_actual[key]]

    # Try swapping only the first namespace separator '.' -> '__'
    if "." in column_name:
        namespace, rest = column_name.split(".", 1)
        alt = f"{namespace}__{rest}".lower()
        if alt in lower_to_actual:
            return metadata_df[lower_to_actual[alt]]

    raise ValueError(
        f"Required metadata column {column_name!r} not found. "
        "If your metadata uses '__' as a namespace separator, pass the "
        "actual column name explicitly."
    )


_OPTIONAL_DEPENDENCY_INSTALL_COMMANDS = types.MappingProxyType(
    {
        "biocframe": 'pip install "recount3[biocpy]"',
        "genomicranges": 'pip install "recount3[biocpy]"',
        "summarizedexperiment": 'pip install "recount3[biocpy]"',
        "pyBigWig": (
            'pip install "recount3[bigwig]"\n'
            "  conda install -c conda-forge -c bioconda pybigwig"
        ),
        "pybiocfilecache": 'pip install "recount3[pybiocfilecache]"',
        "pybiocfilecache.models": 'pip install "recount3[pybiocfilecache]"',
        "pybiocfilecache.utils": 'pip install "recount3[pybiocfilecache]"',
        "pyarrow": 'pip install "recount3[parquet]"',
        "anndata": 'pip install "recount3[anndata]"',
        "delayedarray": 'pip install "recount3[anndata]"',
    }
)

_PARQUET_INSTALL_COMMAND = _OPTIONAL_DEPENDENCY_INSTALL_COMMANDS["pyarrow"]


def _format_optional_dependency_import_error(
    module_name: str,
    exc: BaseException | None = None,
) -> str:
    """Format a standardized ImportError message for an optional dependency.

    Args:
        module_name: Import name used by Python (for example, ``pyBigWig``).
        exc: Underlying exception raised during import, if available.

    Returns:
        A user-facing error message suitable for raising as ImportError.
    """
    command = _OPTIONAL_DEPENDENCY_INSTALL_COMMANDS.get(
        module_name,
        f"pip install {module_name}",
    )

    detail = ""
    if exc is not None:
        detail = f"\n\nOriginal import error: {exc!r}"

    return (
        f"Optional dependency {module_name!r} is required for this feature."
        f"{detail}\n\nInstall it with:\n\n  {command}\n"
    )


def _format_optional_dependency_import_failure(
    module_name: str,
    exc: BaseException,
) -> str:
    """Return an error message for an optional dependency that failed to import.

    Some optional dependencies are native extensions. In those cases, importing
    the module can fail even when it is installed (for example, due to missing
    shared libraries). This helper surfaces the original failure while still
    including installation guidance.

    Args:
        module_name: Import name used by Python.
        exc: The underlying exception raised during import.

    Returns:
        A detailed message suitable for raising as an :exc:`ImportError`.
    """
    return (
        f"Optional dependency {module_name!r} could not be imported.\n"
        f"Import error: {exc!r}\n\n"
        f"{_format_optional_dependency_import_error(module_name)}"
    )


@functools.lru_cache(maxsize=None)
def import_optional_module(module_name: str) -> types.ModuleType:
    """Import and cache an optional dependency.

    This is the single entry point for optional runtime imports.

    Args:
        module_name: Import name used by Python (for example, "biocframe").

    Returns:
        The imported module.

    Raises:
        ImportError: If the dependency is missing or fails to import.
    """
    try:
        return importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        raise ImportError(
            _format_optional_dependency_import_error(module_name),
        ) from exc
    except Exception as exc:  # pylint: disable=broad-except
        raise ImportError(
            _format_optional_dependency_import_failure(module_name, exc),
        ) from exc


def _get_module_attribute(
    module: types.ModuleType,
    attribute_name: str,
    *,
    module_name: str,
) -> Any:
    """Return an attribute from an imported module with a stable error message.

    Args:
        module: Imported module returned by import_optional_module.
        attribute_name: Attribute to retrieve from the module.
        module_name: Import name used to load the module. This is used only for
          error messaging.

    Returns:
        The attribute value.

    Raises:
        ImportError: If the optional dependency is missing or fails to import.
    """
    try:
        return getattr(module, attribute_name)
    except AttributeError as exc:
        raise ImportError(
            _format_optional_dependency_import_error(module_name, exc),
        ) from exc


def get_biocframe_class() -> type["biocframe.BiocFrame"]:
    """Return the BiocPy "biocframe.BiocFrame" class.

    Returns:
        The "biocframe.BiocFrame" class.

    Raises:
        ImportError: If the optional dependency is missing or fails to import.
    """
    module = import_optional_module("biocframe")
    return cast(
        type["biocframe.BiocFrame"],
        _get_module_attribute(
            module,
            "BiocFrame",
            module_name="biocframe",
        ),
    )


def get_genomicranges_class() -> type["genomicranges.GenomicRanges"]:
    """Return the BiocPy "genomicranges.GenomicRanges" class.

    Returns:
        The "genomicranges.GenomicRanges" class.

    Raises:
        ImportError: If the optional dependency is missing or fails to import.
    """
    module = import_optional_module("genomicranges")
    return cast(
        type["genomicranges.GenomicRanges"],
        _get_module_attribute(
            module,
            "GenomicRanges",
            module_name="genomicranges",
        ),
    )


def get_summarizedexperiment_class() -> (
    type["summarizedexperiment.SummarizedExperiment"]
):
    """Return the BiocPy "summarizedexperiment.SummarizedExperiment" class.

    Returns:
        The "summarizedexperiment.SummarizedExperiment" class.

    Raises:
        ImportError: If the optional dependency is missing or fails to import.
    """
    module = import_optional_module("summarizedexperiment")
    return cast(
        type["summarizedexperiment.SummarizedExperiment"],
        _get_module_attribute(
            module,
            "SummarizedExperiment",
            module_name="summarizedexperiment",
        ),
    )


def get_ranged_summarizedexperiment_class() -> (
    type["summarizedexperiment.RangedSummarizedExperiment"]
):
    """Return the BiocPy RangedSummarizedExperiment class.

    Returns:
        The "summarizedexperiment.RangedSummarizedExperiment" class.

    Raises:
        ImportError: If the optional dependency is missing or fails to import.
    """
    module = import_optional_module("summarizedexperiment")
    return cast(
        type["summarizedexperiment.RangedSummarizedExperiment"],
        _get_module_attribute(
            module,
            "RangedSummarizedExperiment",
            module_name="summarizedexperiment",
        ),
    )


def get_pybigwig_module() -> types.ModuleType:
    """Return the optional ``pyBigWig`` module.

    This is a small convenience wrapper around import_optional_module so that
    callers do not need to hard-code the import name.

    Returns:
        The imported ``pyBigWig`` module.

    Raises:
        ImportError: If the optional dependency is missing or fails to import.
    """
    return import_optional_module("pyBigWig")


def get_anndata_module() -> types.ModuleType:
    """Return the optional ``anndata`` module.

    Returns:
        The imported ``anndata`` module.

    Raises:
        ImportError: If the optional dependency is missing or fails to import.
    """
    return import_optional_module("anndata")


def ensure_anndata_support() -> None:
    """Verify that AnnData conversion is possible before doing expensive work.

    :meth:`summarizedexperiment.SummarizedExperiment.to_anndata` imports both
    ``anndata`` and ``delayedarray``, and ``summarizedexperiment`` declares
    neither as a required dependency. Probing both up front lets callers fail
    before downloading and assembling data that cannot then be written.

    Raises:
        ImportError: If either dependency is missing or fails to import.
    """
    for module_name in ("anndata", "delayedarray"):
        import_optional_module(module_name)


def _anndata_frames(adata: Any) -> list[tuple[str, Any]]:
    """Return an AnnData object's ``obs`` and ``var`` frames, in that order.

    The AnnData helpers below accept :data:`~typing.Any` because they also run
    against test doubles, so each pair is returned only when the attribute is
    present and frame-shaped.

    Args:
        adata: AnnData object to read ``obs`` and ``var`` from.

    Returns:
        ``(name, frame)`` pairs for whichever of the two frames are present.
    """
    frames: list[tuple[str, Any]] = []
    for frame_name in ("obs", "var"):
        frame = getattr(adata, frame_name, None)
        if frame is not None and hasattr(frame, "columns"):
            frames.append((frame_name, frame))
    return frames


def normalize_anndata_for_hdf5(adata: Any) -> list[str]:
    """Cast all-missing object columns so an AnnData object can be written.

    :mod:`h5py` cannot serialize an object-dtype column whose values are all
    :data:`None`: there is no string for ``anndata`` to infer a type from, and
    the write fails with ``TypeError: Can't implicitly convert non-string
    objects to strings``. recount3 produces such columns whenever a GTF
    attribute (for example, ``phase``) or a sample-metadata field is absent for
    every row, so real projects hit this on the default export path.

    Casting those columns to all-NaN ``float64`` preserves "missing for every
    row" and round-trips through :func:`anndata.read_h5ad`. Columns with at
    least one string are already written correctly and are left alone, as are
    columns holding genuinely non-string objects, which still fail loudly.

    Args:
        adata: AnnData object. Its ``obs`` and ``var`` frames are modified
          in place.

    Returns:
        ``"frame.column"`` labels that were cast, in ``obs`` then ``var``
        order, for logging.
    """
    converted: list[str] = []
    for frame_name, frame in _anndata_frames(adata):
        for column in frame.columns:
            series = frame[column]
            if series.dtype != object or series.empty:
                continue
            if series.isna().all():
                frame[column] = series.astype("float64")
                converted.append(f"{frame_name}.{column}")
    return converted


def hdf5_unsafe_column_names(adata: Any) -> list[str]:
    """Return ``obs``/``var`` column names that HDF5 cannot use as group keys.

    ``anndata`` stores each column of ``obs`` and ``var`` as an HDF5 dataset
    named after the column, and HDF5 treats ``"/"`` as a path separator, so a
    column name containing one fails with ``ValueError: Forward slashes are
    not allowed in keys``. recount3 metadata hits this: the STAR QC fields are
    named after splice-site motifs, for example
    ``recount_qc__star.number_of_splices:_gt/ag``.

    Args:
        adata: AnnData object to inspect.

    Returns:
        ``"frame.column"`` labels containing a forward slash, in ``obs`` then
        ``var`` order.
    """
    unsafe: list[str] = []
    for frame_name, frame in _anndata_frames(adata):
        unsafe.extend(
            f"{frame_name}.{column}"
            for column in frame.columns
            if "/" in str(column)
        )
    return unsafe


def sanitize_anndata_column_names(adata: Any) -> list[tuple[str, str]]:
    """Replace forward slashes in ``obs``/``var`` column names with ``"_"``.

    This renames the columns an analysis indexes by, so callers should apply it
    only on an explicit request and report every rename.

    Args:
        adata: AnnData object. Its ``obs`` and ``var`` frames are renamed
          in place.

    Returns:
        ``(old, new)`` name pairs that were applied, in ``obs`` then ``var``
        order.

    Raises:
        ValueError: If a sanitized name would collide with another column in
          the same frame, which would silently merge two distinct fields.
    """
    renames: list[tuple[str, str]] = []
    for frame_name, frame in _anndata_frames(adata):
        mapping: dict[Any, str] = {}
        taken = {
            str(column) for column in frame.columns if "/" not in str(column)
        }
        for column in frame.columns:
            name = str(column)
            if "/" not in name:
                continue
            new_name = name.replace("/", "_")
            if new_name in taken:
                raise ValueError(
                    f"Cannot sanitize {frame_name} column {name!r} for HDF5: "
                    f"the sanitized name {new_name!r} is already used by "
                    "another column, and renaming would merge two distinct "
                    "fields. Write a .pkl instead, or rename the column "
                    "before exporting."
                )
            taken.add(new_name)
            mapping[column] = new_name
            renames.append((name, new_name))

        if mapping:
            frame.rename(columns=mapping, inplace=True)
    return renames


def _plain_metadata_value(value: Any) -> Any:
    """Return ``value`` rebuilt from containers :mod:`h5py` can write.

    ``h5py`` has no writer for a :class:`tuple`, and ``anndata`` rejects any
    ``uns`` mapping that is not a mutable mapping, so every nested container
    is rebuilt as a plain :class:`dict` or :class:`list`. Scalars pass
    through untouched.

    Args:
        value: A metadata value of any type.

    Returns:
        An equivalent value built only from dicts, lists, and scalars.
    """
    if hasattr(value, "as_dict"):
        value = value.as_dict()
    if isinstance(value, Mapping):
        return {str(k): _plain_metadata_value(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain_metadata_value(item) for item in value]
    return value


def experiment_metadata_as_dict(experiment: Any) -> dict[str, Any]:
    """Return an experiment's provenance metadata as a plain, writable dict.

    BiocPy stores experiment metadata as a
    :class:`~biocutils.NamedList.NamedList`, and recount3's provenance nests
    tuples inside it (``metadata_columns`` maps each namespaced sample-metadata
    column to its ``(table, column)`` origin). Neither survives a round trip
    through ``anndata``, so both are converted here.

    Args:
        experiment: A ``(Ranged)SummarizedExperiment``-like object exposing
          ``get_metadata()``.

    Returns:
        The metadata as a :class:`dict` of dicts, lists, and scalars. Empty
        when the experiment carries no metadata, which BiocPy represents as
        an unnamed ``NamedList``. It has no keys to map into ``uns``, and
        its ``as_dict()`` raises rather than returning an empty dict.
    """
    metadata = experiment.get_metadata()
    if metadata is None:
        return {}
    if hasattr(metadata, "as_dict"):
        if getattr(metadata, "get_names", lambda: None)() is None:
            return {}
        metadata = metadata.as_dict()
    return {
        str(key): _plain_metadata_value(value)
        for key, value in dict(metadata).items()
    }


def experiment_to_anndata(experiment: Any) -> Any:
    """Convert an experiment to AnnData, carrying provenance into ``uns``.

    ``SummarizedExperiment.to_anndata()`` forwards its metadata straight to
    ``AnnData(uns=...)``, which rejects anything that is not a mutable
    mapping. Because BiocPy holds that metadata as a ``NamedList``, the direct
    call fails for every experiment recount3 builds. Converting the metadata
    first and assigning ``uns`` afterwards keeps the provenance without
    depending on how BiocPy chooses to store it.

    The experiment itself is not modified: ``set_metadata`` returns a copy
    that shares the assays, so nothing is duplicated in memory.

    Args:
        experiment: A ``(Ranged)SummarizedExperiment`` to convert.

    Returns:
        An ``anndata.AnnData`` with samples in rows, genes in columns, and
        the experiment's provenance in ``uns``.
    """
    metadata = experiment_metadata_as_dict(experiment)
    adata = experiment.set_metadata(None).to_anndata()
    adata.uns = metadata
    return adata


def _iter_uns_mappings(
    uns: Any, prefix: str = "uns"
) -> Iterator[tuple[str, Any]]:
    """Yield every nested mapping inside ``uns`` with its dotted path."""
    if isinstance(uns, Mapping):
        yield prefix, uns
        for key, value in uns.items():
            yield from _iter_uns_mappings(value, f"{prefix}.{key}")
    elif isinstance(uns, list):
        for index, value in enumerate(uns):
            yield from _iter_uns_mappings(value, f"{prefix}[{index}]")


def hdf5_unsafe_uns_keys(adata: Any) -> list[str]:
    """Return ``uns`` keys that HDF5 cannot use as group keys.

    ``uns`` becomes a group hierarchy in the HDF5 file, so the same forward
    slash that breaks an ``obs`` column name breaks a key here. recount3 hits
    this through ``uns["metadata_columns"]``, which is keyed by the very
    sample-metadata column names that :func:`hdf5_unsafe_column_names`
    reports.

    Args:
        adata: AnnData object to inspect.

    Returns:
        Dotted ``"uns.<path>.<key>"`` labels containing a forward slash, in
        traversal order.
    """
    unsafe: list[str] = []
    for path, mapping in _iter_uns_mappings(getattr(adata, "uns", None)):
        unsafe.extend(f"{path}.{key}" for key in mapping if "/" in str(key))
    return unsafe


def sanitize_anndata_uns_keys(adata: Any) -> list[tuple[str, str]]:
    """Replace forward slashes in nested ``uns`` keys with ``"_"``.

    This is the ``uns`` counterpart of
    :func:`sanitize_anndata_column_names` and belongs to the same explicit
    request: ``uns["metadata_columns"]`` is keyed by the ``obs`` column names,
    so renaming one without the other would leave the provenance map pointing
    at columns that no longer exist. The unsanitized name is not lost;
    each value records its ``(table, column)`` origin verbatim.

    Args:
        adata: AnnData object. Its ``uns`` mappings are rewritten in place.

    Returns:
        ``(old, new)`` key pairs that were applied, in traversal order.

    Raises:
        ValueError: If a sanitized key would collide with another key in the
          same mapping, which would silently merge two distinct entries.
    """
    renames: list[tuple[str, str]] = []
    for path, mapping in _iter_uns_mappings(getattr(adata, "uns", None)):
        slashed = [key for key in mapping if "/" in str(key)]
        if not slashed:
            continue
        taken = {str(key) for key in mapping if "/" not in str(key)}
        for key in slashed:
            name = str(key)
            new_name = name.replace("/", "_")
            if new_name in taken:
                raise ValueError(
                    f"Cannot sanitize {path} key {name!r} for HDF5: the "
                    f"sanitized key {new_name!r} is already used in the same "
                    "mapping, and renaming would merge two distinct entries. "
                    "Write a .pkl instead, which keeps the keys verbatim."
                )
            taken.add(new_name)
            mapping[new_name] = mapping.pop(key)
            renames.append((name, new_name))
    return renames


_PARQUET_IMPL_ENGINE_NAMES = types.MappingProxyType(
    {
        "PyArrowImpl": "pyarrow",
        "FastParquetImpl": "fastparquet",
    }
)


def _format_parquet_engine_error(exc: BaseException) -> str:
    """Format an actionable error for an unresolvable pandas Parquet engine.

    Args:
        exc: The :exc:`ImportError` pandas raised while resolving an engine.

    Returns:
        A user-facing error message suitable for raising as an
        :exc:`ImportError`.
    """
    configured = "auto"
    with contextlib.suppress(Exception):
        configured = str(pd.get_option("io.parquet.engine"))

    tried = (
        f"pandas is configured to use the {configured!r} engine "
        "(io.parquet.engine)"
        if configured != "auto"
        else "pandas tried the 'pyarrow' and 'fastparquet' engines"
    )

    return (
        f"Writing Parquet requires a Parquet engine, but {tried} and none "
        "is usable.\n\n"
        f"Install one with:\n\n  {_PARQUET_INSTALL_COMMAND}\n\n"
        "Or write a text format instead, by choosing an output path ending "
        "in .tsv, .tsv.gz, or .csv.\n\n"
        f"Original import error: {exc!r}"
    )


def ensure_parquet_engine() -> str:
    """Verify that pandas can resolve a Parquet engine, and name it.

    This defers engine selection to pandas rather than probing import names
    directly, so it honours the ``io.parquet.engine`` option, pandas' own
    minimum-version rules, and a ``fastparquet``-only installation.

    Returns:
        The resolved engine name (``"pyarrow"`` or ``"fastparquet"``), or
        ``"auto"`` when the engine cannot be identified but pandas accepted it.

    Raises:
        ImportError: If no usable Parquet engine is installed.
        ValueError: If ``io.parquet.engine`` is set to an unknown engine.
    """
    try:
        # pylint: disable-next=import-outside-toplevel
        from pandas.io.parquet import get_engine
    except ImportError:
        return "auto"

    try:
        impl = get_engine("auto")
    except ImportError as exc:
        raise ImportError(_format_parquet_engine_error(exc)) from exc

    return _PARQUET_IMPL_ENGINE_NAMES.get(type(impl).__name__, "auto")


def sparse_column_names(frame: pd.DataFrame) -> list[str]:
    """Return the names of columns backed by a pandas sparse dtype.

    Junction MM resources load as sparse-backed DataFrames (see
    :meth:`recount3.resource.R3Resource.load`), and no Parquet engine accepts
    :class:`pandas.SparseDtype` columns.

    Args:
        frame: DataFrame to inspect.

    Returns:
        Column names with a sparse dtype, in column order.
    """
    return [
        str(name)
        for name, dtype in frame.dtypes.items()
        if isinstance(dtype, pd.SparseDtype)
    ]


def densify_sparse_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Return ``frame`` with every sparse column converted to its dense dtype.

    Each sparse column becomes a dense column of its
    :attr:`pandas.SparseDtype.subtype`. Densifying a junction matrix
    materializes every implicit zero, so callers should treat this as an
    explicit, opt-in memory cost.

    Args:
        frame: DataFrame that may contain sparse columns.

    Returns:
        The original object when no column is sparse, otherwise a new
        DataFrame with dense columns.
    """
    conversions = {
        name: dtype.subtype
        for name, dtype in frame.dtypes.items()
        if isinstance(dtype, pd.SparseDtype)
    }
    if not conversions:
        return frame
    return frame.astype(conversions)


_JXN_SIDECAR_RE = re.compile(r"\.(MM|ID|RR)\.gz$", re.IGNORECASE)


def _derive_junction_sidecar_url(url: str, new_ext: str) -> str:
    """Return a junction sidecar URL by swapping the ``.{MM,ID,RR}.gz`` suffix.

    Junction files come in a triplet sharing a common stem: a MatrixMarket
    matrix (``.MM.gz``), a sample-ID table (``.ID.gz``), and a row-ranges
    table (``.RR.gz``). Given the URL for any one of the three, this function
    produces the URL for another member of the triplet.

    Args:
        url: A non-empty URL whose path ends with ``.MM.gz``, ``.ID.gz``, or
            ``.RR.gz`` (case-insensitive).
        new_ext: The target extension token. Must be one of ``"MM"``, ``"ID"``,
            or ``"RR"`` (case-insensitive).

    Returns:
        The URL with the trailing ``.<ext>.gz`` replaced by
        ``.<new_ext>.gz``.

    Raises:
        ValueError: If ``url`` is empty, if ``new_ext`` is not one of
            ``MM``/``ID``/``RR``, or if ``url`` does not end with a
            recognised junction suffix.
    """
    if not url:
        raise ValueError("url must be non-empty")
    new_ext_u = new_ext.upper()
    if new_ext_u not in {"MM", "ID", "RR"}:
        raise ValueError(f"new_ext must be one of MM/ID/RR, got {new_ext!r}")
    if _JXN_SIDECAR_RE.search(url):
        return _JXN_SIDECAR_RE.sub(f".{new_ext_u}.gz", url)
    raise ValueError(f"Cannot derive junction sidecar URL from {url!r}")
