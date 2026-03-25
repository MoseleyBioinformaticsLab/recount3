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
#
"""Internal utilities for HTTP, filesystem, and caching operations.

The module is organized into four sections:

1. Caching Utilities: URL-based cache key generation and path management
2. Filesystem Utilities: Directory management and atomic file operations
3. HTTP Utilities: Network requests with retries, streaming, and ZIP handling
4. Optional Dependency Management: lazy import helpers for BiocPy packages
   (``biocframe``, ``summarizedexperiment``, ``genomicranges``) and
   ``pyBigWig``; raises a standardized :exc:`ImportError` when a required
   optional package is missing.

Note:
    This module is considered internal implementation detail and may change
    without notice. Users should prefer the public API in the main package.

Examples:
    Typical usage is through the public API, but direct usage might look like::

        from pathlib import Path
        from recount3._utils import _cache_path

        path = _cache_path("https://example.com/data.tsv", Path("/tmp/cache"))

Attributes:
    _ZIP_LOCKS_GUARD: Threading lock to safely interact with the weakref dictionary.
    _ZIP_LOCKS: A weak reference dictionary mapping ZIP file paths to specific
        threading locks. This synchronizes ZIP mutations per-file, preventing
        race conditions without leaking memory for inactive cache paths.
"""

from __future__ import annotations

import contextlib
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
from pathlib import Path
from typing import BinaryIO, Any, cast, TYPE_CHECKING

import pandas as pd

from recount3 import errors

if TYPE_CHECKING:  # pragma: no cover
    import biocframe  # type: ignore[import-not-found]
    import genomicranges  # type: ignore[import-not-found]
    import summarizedexperiment  # type: ignore[import-not-found]

_ZIP_LOCKS_GUARD = threading.Lock()
_ZIP_LOCKS: weakref.WeakValueDictionary[str, _WeakRefLock] = (
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


def _zip_lock_for_path(zip_path: Path) -> _WeakRefLock:
    """Return the shared lock for a specific ZIP path.

    Args:
        zip_path: Target ZIP file path.

    Returns:
        A lock object shared by all operations writing to this ZIP path.
    """
    key = str(zip_path.expanduser().resolve())
    with _ZIP_LOCKS_GUARD:
        lock = _ZIP_LOCKS.get(key)
        if lock is None:
            lock = _WeakRefLock()
            _ZIP_LOCKS[key] = lock
        return lock





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
        raise NotADirectoryError(
            f"Exists but is not a directory: {p}"
        ) from exc


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
    tmp_dst = (
        dst.parent
        / f".{dst.name}.{os.getpid()}_{threading.get_ident()}_{time.time_ns()}.tmp"
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
        insecure_ssl: If True, disable TLS verification.
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

            tmp = (
                out_path.parent
                / f".{out_path.name}.{os.getpid()}_{threading.get_ident()}_{time.time_ns()}.downloading"
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
                f".{zip_path.name}.{os.getpid()}_{threading.get_ident()}_{time.time_ns()}.tmpzip"
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
    """Coerce sample metadata into a pandas DataFrame.

    Args:
        sample_metadata_source: Either a BiocPy (Ranged)SummarizedExperiment-like object with a
          `.col_data.to_pandas()` method, or a pandas DataFrame already.

    Returns:
        A pandas DataFrame of sample metadata.

    Raises:
        TypeError: If `sample_metadata_source` cannot be coerced to a pandas DataFrame.
    """
    if isinstance(sample_metadata_source, pd.DataFrame):
        return sample_metadata_source

    if hasattr(sample_metadata_source, "col_data") and hasattr(
        sample_metadata_source.col_data, "to_pandas"
    ):
        return sample_metadata_source.col_data.to_pandas()

    raise TypeError(
        "Expected a pandas.DataFrame or a SummarizedExperiment-like object with "
        "`col_data.to_pandas()`."
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
        "biocframe": "pip install biocframe genomicranges summarizedexperiment",
        "genomicranges": "pip install biocframe genomicranges summarizedexperiment",
        "summarizedexperiment": "pip install biocframe genomicranges summarizedexperiment",
        "pyBigWig": (
            "pip install pyBigWig\n"
            "  conda install -c conda-forge -c bioconda pybigwig"
        ),
    }
)


def _format_optional_dependency_import_error(
    module_name: str,
    exc: BaseException | None = None,
) -> str:
    """Format a standardized ImportError message for an optional dependency.

    Args:
        module_name: Import name used by Python (for example, "pyBigWig").
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
        A detailed message suitable for raising as a CompatibilityError.
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
    """Return the BiocPy "summarizedexperiment.RangedSummarizedExperiment" class.

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
    """Return the optional "pyBigWig" module.

    This is a small convenience wrapper around import_optional_module so that
    callers do not need to hard-code the import name.

    Returns:
        The imported "pyBigWig" module.

    Raises:
        ImportError: If the optional dependency is missing or fails to import.
    """
    return import_optional_module("pyBigWig")




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
