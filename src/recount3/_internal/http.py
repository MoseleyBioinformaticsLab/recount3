"""HTTP and streaming helpers (internal).

These helpers are parameterized and do not read environment variables. The
public API passes values from :mod:`recount3.config` to preserve behavior.
"""

from __future__ import annotations

import contextlib
import socket
import ssl
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from typing import BinaryIO

from .fs import _atomic_replace, _ensure_dir
from ..errors import DownloadError


def _ssl_insecure_context() -> ssl.SSLContext:
    """Return an SSL context with verification disabled (not recommended)."""
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
    """Open a URL with urllib, ensuring all headers are valid strings.

    Args:
      url: Absolute URL to open.
      timeout: Socket timeout in seconds.
      headers: Optional extra headers to send. ``None`` values are dropped.
      insecure_ssl: If True, disable TLS verification.
      user_agent: User-Agent header override.

    Returns:
      A binary HTTP response object.

    Raises:
      URLError, HTTPError, socket.timeout: Network/HTTP errors from urllib.
      ValueError: If ``url`` is empty.
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
    """Execute ``func`` with simple exponential backoff on transient errors."""
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
    """Copy bytes in chunks from ``src_fh`` to ``dst_fh``."""
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
    """Download ``url`` atomically into ``out_path``.

    Network reads stream directly to a temp file which is atomically moved
    into place to avoid partial files on failure.
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
    except Exception as exc:  # Re-wrap for clarity.
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
    """Stream URL content directly into a ``.zip`` file (no cache involvement)."""
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
    """Write an existing file into ``zip_path`` under ``arcname``."""
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
