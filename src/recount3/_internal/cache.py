"""Cache path helpers (internal)."""

from __future__ import annotations

import hashlib
from pathlib import Path

from .fs import _ensure_dir


def _sha256(text: str) -> str:
    """Return the hex SHA256 of ``text``."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _cache_key_for_url(url: str) -> str:
    """Return a stable cache key derived from the full URL."""
    # Keep the original naming scheme: short digest + basename.
    from urllib.parse import urlparse  # local import to keep module tiny

    parsed = urlparse(url)
    key = parsed.path.lstrip("/")
    digest = _sha256(url)[:16]
    return f"{digest}__{Path(key).name}"


def _cache_path(url: str, cache_root: str | Path) -> Path:
    """Return the full cache path for a given URL under ``cache_root``."""
    _ensure_dir(cache_root)
    return Path(cache_root) / _cache_key_for_url(url)
