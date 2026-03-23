"""Configuration helpers for recount3.

This module centralizes configuration so that environment-dependent values
are not hidden as mutable module globals. Values are read once via
:func:`default_config` and can be overridden by constructing :class:`Config`
directly. CLI flags in :mod:`recount3.cli` take the highest precedence and
override both environment variables and :class:`Config` defaults.

The module also exposes three cache utility functions:
:func:`recount3_cache`, :func:`recount3_cache_files`, and
:func:`recount3_cache_rm`.

Environment variables (all optional):
  * ``RECOUNT3_URL``: base URL of the duffel mirror
    (default: ``https://duffel.rail.bio/recount3/``).
  * ``RECOUNT3_CACHE_DIR``: directory for the on-disk file cache
    (default: ``~/.cache/recount3``).
  * ``RECOUNT3_CACHE_DISABLE``: set to ``"1"`` to disable caching entirely.
  * ``RECOUNT3_HTTP_TIMEOUT``: HTTP request timeout in seconds (default: 30).
  * ``RECOUNT3_MAX_RETRIES``: maximum retry attempts for transient errors
    (default: 3).
  * ``RECOUNT3_INSECURE_SSL``: set to ``"1"`` to skip TLS verification
    (unsafe; use only for debugging certificate issues).
  * ``RECOUNT3_USER_AGENT``: custom ``User-Agent`` header string.
  * ``RECOUNT3_CHUNK_SIZE``: streaming chunk size in bytes (default: 65536).

Typical usage example::

    from pathlib import Path
    from recount3.config import Config, recount3_cache, recount3_cache_rm

    # Read the current cache directory (creates it if absent):
    cache_dir = recount3_cache()

    # Use a custom cache location for this session:
    cfg = Config(cache_dir=Path("/scratch/recount3_cache"))

    # Remove cached files matching a pattern (dry run first):
    to_delete = recount3_cache_rm(dry_run=True)
    recount3_cache_rm(predicate=lambda p: "sra" in str(p))
"""

from __future__ import annotations

import os
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from recount3.version import __version__

from recount3._utils import _ensure_dir

_DEFAULT_CHUNK_SIZE: int = 1024 * 1024  # 1 MiB


@dataclass(frozen=True, slots=True)
class Config:
    """Immutable configuration bag.

    Attributes:
      base_url: Base URL for the duffel mirror (ends with a trailing slash).
      timeout: Network timeout in seconds.
      insecure_ssl: True to disable TLS verification (not recommended).
      max_retries: Max HTTP retry attempts for transient errors.
      user_agent: Custom HTTP User-Agent.
      cache_dir: Cache directory for downloaded files.
      cache_disabled: If True, disable cache behavior globally.
      chunk_size: Default chunk size in bytes for streaming copies.
    """

    base_url: str
    timeout: int
    insecure_ssl: bool
    max_retries: int
    user_agent: str
    cache_dir: Path
    cache_disabled: bool
    chunk_size: int = _DEFAULT_CHUNK_SIZE


def default_config() -> Config:
    """Return configuration constructed from environment variables.

    Returns:
      A :class:`Config` populated from the environment.

    Notes:
      Values are parsed to sensible types and the base URL is normalized to
      include a trailing slash (matching the original behavior).
    """
    base = (
        os.environ.get(
            "RECOUNT3_URL", "http://duffel.rail.bio/recount3/"
        ).rstrip("/")
        + "/"
    )
    cache_dir = Path(
        os.environ.get(
            "RECOUNT3_CACHE_DIR",
            os.path.join(
                os.path.expanduser("~"), ".cache", "recount3", "files"
            ),
        )
    )
    return Config(
        base_url=base,
        timeout=int(os.environ.get("RECOUNT3_HTTP_TIMEOUT", "60")),
        insecure_ssl=os.environ.get("RECOUNT3_INSECURE_SSL", "0") == "1",
        max_retries=int(os.environ.get("RECOUNT3_MAX_RETRIES", "3")),
        user_agent=(
            os.environ.get("RECOUNT3_USER_AGENT")
            or (
                f"recount3-python/{__version__} "
                "(+https://github.com/MoseleyBioinformaticsLab/recount3)"
            )
        ),
        cache_dir=cache_dir,
        cache_disabled=os.environ.get("RECOUNT3_CACHE_DISABLE", "0") == "1",
        chunk_size=int(
            os.environ.get("RECOUNT3_CHUNK_SIZE", str(_DEFAULT_CHUNK_SIZE))
        ),
    )


def recount3_cache(config: Config | None = None) -> Path:
    """Return the cache directory used for recount3 downloads.

    This helper normalizes and materializes the cache directory based on
    the provided configuration (or the default configuration when omitted).

    Args:
      config: Optional configuration. If None, :func:`default_config` is
        used.

    Returns:
      Absolute :class:`pathlib.Path` to the cache directory.
    """
    cfg = config or default_config()
    _ensure_dir(cfg.cache_dir)
    return cfg.cache_dir


def recount3_cache_files(
    config: Config | None = None,
    *,
    pattern: str | None = None,
) -> list[Path]:
    """List cached files managed by recount3.

    Args:
      config: Optional configuration. If None, :func:`default_config` is
        used.
      pattern: Optional glob-style pattern (as accepted by
        :meth:`pathlib.Path.rglob`) to filter files relative to the cache
        root, for example ``"*.tsv.gz"`` or ``"*__SRP123456*"``. If None,
        all files are returned.

    Returns:
      A sorted list of :class:`pathlib.Path` objects pointing to cached
      files. If the cache directory does not exist yet, an empty list is
      returned.
    """
    cfg = config or default_config()
    root = cfg.cache_dir

    if not root.exists() or not root.is_dir():
        return []

    glob_pattern = pattern if pattern is not None else "*"
    files: list[Path] = []

    for path in root.rglob(glob_pattern):
        if path.is_file():
            files.append(path)

    return sorted(files, key=str)


def recount3_cache_rm(
    *,
    config: Config | None = None,
    predicate: Callable[[Path], bool] | None = None,
    dry_run: bool = False,
) -> list[Path]:
    """Remove cached files that match a predicate.

    This helper is analogous to the R-side ``recount3_cache_rm()``: it
    walks the cache directory and removes any file for which ``predicate``
    returns True. Directories are left in place.

    Args:
      config: Optional configuration. If None, :func:`default_config` is
        used.
      predicate: Callable taking a :class:`pathlib.Path` and returning
        True if the file should be removed. If None, all cached files are
        selected.
      dry_run: If True, do not delete any files and only report which
        paths would be removed.

    Returns:
      A sorted list of :class:`pathlib.Path` objects that were removed (or
      would be removed when ``dry_run`` is True).

    Raises:
      OSError: If filesystem operations fail during deletion.

    Examples:
        List everything that would be removed, without deleting::

            to_delete = recount3_cache_rm(dry_run=True)

        Remove all cached files (empty the cache)::

            recount3_cache_rm()

        Remove only files related to the ``"sra"`` data source::

            recount3_cache_rm(predicate=lambda p: "sra" in str(p))
    """
    cfg = config or default_config()
    root = cfg.cache_dir

    if not root.exists() or not root.is_dir():
        return []

    def _select(path: Path) -> bool:
        if predicate is None:
            return True
        return predicate(path)

    candidates: list[Path] = []
    for path in root.rglob("*"):
        if path.is_file() and _select(path):
            candidates.append(path)

    candidates = sorted(candidates, key=str)

    if dry_run:
        return candidates

    for path in candidates:
        path.unlink()

    return candidates
