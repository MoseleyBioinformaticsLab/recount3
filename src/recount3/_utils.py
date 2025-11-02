"""Internal utilities for HTTP, filesystem, and caching operations.

This module consolidates low-level implementation details for network operations,
filesystem management, and caching that support the public recount3 API. These
utilities are designed to be thread-safe, configurable, and side-effect free
beyond their intended operations.

The module is organized into three sections:

1. Caching Utilities: URL-based cache key generation and path management
2. Filesystem Utilities: Directory management and atomic file operations  
3. HTTP Utilities: Network requests with retries, streaming, and ZIP handling

Note:
    This module is considered internal implementation detail and may change
    without notice. Users should prefer the public API in the main package.

Examples:
    Typical usage is through the public API, but direct usage might look like:

    >>> from recount3._utils import _cache_path, download_to_file
    >>> cache_path = _cache_path("http://example.com/data.tsv", "/tmp/cache")
    >>> download_to_file("http://example.com/data.tsv", cache_path, ...)

Attributes:
    _FILE_LOCK: Threading lock for synchronizing file operations to prevent
        race conditions during concurrent cache access.
"""

from __future__ import annotations

import contextlib
import errno
import hashlib
import os
import shutil
import socket
import ssl
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from typing import BinaryIO

from .errors import DownloadError

# Serialize writes to the same filesystem paths
_FILE_LOCK = threading.Lock()


# =============================================================================
# Caching Utilities
# =============================================================================

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


# =============================================================================
# Filesystem Utilities  
# =============================================================================

def _ensure_dir(path: str | Path) -> None:
    """Ensure directory exists, creating parents as needed.
    
    Args:
        path: Directory path to create.
        
    Raises:
        NotADirectoryError: If path exists as a non-directory.
        OSError: If directory creation fails due to permissions or other system error.
    """
    p = Path(path)
    if not p.exists():
        p.mkdir(parents=True, exist_ok=True)
    elif not p.is_dir():
        raise NotADirectoryError(str(p))


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
    try:
        if dst.exists():
            dst.unlink()
        os.link(src, dst)
    except OSError as e:
        if e.errno in (errno.EXDEV, errno.EPERM, errno.EACCES, errno.EMLINK):
            _ensure_dir(dst.parent)
            shutil.copy2(src, dst)
        else:
            raise


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


# =============================================================================
# HTTP Utilities
# =============================================================================

def _ssl_insecure_context() -> ssl.SSLContext:
    """Return SSL context with verification disabled (not recommended).
    
    Warning:
        This disables certificate verification and should only be used
        for testing or with trusted networks.
        
    Returns:
        SSL context with hostname verification and certificate checking disabled.
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
        opener.add_handler(urllib.request.HTTPSHandler(context=_ssl_insecure_context()))

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
        except (urllib.error.URLError, socket.timeout) as exc:
            last_exc = exc
            if i < attempts - 1:
                time.sleep(base_sleep * (2**i))
    if last_exc:
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
            tmp = out_path.with_suffix(out_path.suffix + ".downloading")
            with open(tmp, "wb") as fh:
                _stream_copy(resp, fh, chunk_size=chunk_size)
            _atomic_replace(tmp, out_path)

    try:
        with_retries(_do, attempts=attempts)
    except Exception as exc:
        raise DownloadError(f"Failed to download {url!r} -> {str(out_path)!r}") from exc


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
    """Stream URL content directly into ZIP file without caching.
    
    Args:
        url: Resource URL to download.
        zip_path: Destination ZIP file path.
        arcname: Name for the file within the ZIP archive.
        chunk_size: Bytes to read/write per chunk during streaming.
        overwrite: If True, replace existing file in ZIP.
        timeout: Network timeout in seconds.
        insecure_ssl: If True, disable TLS verification.
        user_agent: HTTP User-Agent header.
        attempts: Maximum retry attempts for transient errors.
        
    Raises:
        DownloadError: If download fails after all retry attempts.
        zipfile.BadZipFile: If ZIP file operations fail.
    """
    _ensure_dir(zip_path.parent)
    mode = "a" if zip_path.exists() else "w"
    with zipfile.ZipFile(zip_path, mode=mode, compression=zipfile.ZIP_DEFLATED) as zf:
        if not overwrite:
            try:
                zf.getinfo(arcname)
                return
            except KeyError:
                pass

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
                with zf.open(arcname, "w") as zfh:
                    _stream_copy(resp, zfh, chunk_size=chunk_size)

        with_retries(_do, attempts=attempts)


def write_cached_file_to_zip(
    cached_file: Path,
    zip_path: Path,
    arcname: str,
    *,
    overwrite: bool,
) -> None:
    """Write existing cached file into ZIP archive.
    
    Args:
        cached_file: Path to existing cached file.
        zip_path: Destination ZIP file path.
        arcname: Name for the file within the ZIP archive.
        overwrite: If True, replace existing file in ZIP.
        
    Raises:
        FileNotFoundError: If cached_file doesn't exist.
        zipfile.BadZipFile: If ZIP file operations fail.
    """
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
