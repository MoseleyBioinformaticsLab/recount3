"""Configuration helpers for recount3.

This module centralizes configuration so that environment-dependent values
are not hidden as mutable module globals. Values are read once via
:func:`default_config` and can be overridden by constructing :class:`Config`.

The defaults match the original script's behavior.

Environment variables:
  * RECOUNT3_URL
  * RECOUNT3_CACHE_DIR
  * RECOUNT3_CACHE_DISABLE
  * RECOUNT3_HTTP_TIMEOUT
  * RECOUNT3_MAX_RETRIES
  * RECOUNT3_INSECURE_SSL
  * RECOUNT3_USER_AGENT
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path


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
    chunk_size: int = 1024 * 1024  # 1 MiB


def default_config() -> Config:
    """Return configuration constructed from environment variables.

    Returns:
      A :class:`Config` populated from the environment.

    Notes:
      Values are parsed to sensible types and the base URL is normalized to
      include a trailing slash (matching the original behavior).
    """
    base = os.environ.get("RECOUNT3_URL", "http://duffel.rail.bio/recount3/").rstrip("/") + "/"
    cache_dir = Path(
        os.environ.get(
            "RECOUNT3_CACHE_DIR",
            os.path.join(os.path.expanduser("~"), ".cache", "recount3", "files"),
        )
    )
    return Config(
        base_url=base,
        timeout=int(os.environ.get("RECOUNT3_HTTP_TIMEOUT", "60")),
        insecure_ssl=os.environ.get("RECOUNT3_INSECURE_SSL", "0") == "1",
        max_retries=int(os.environ.get("RECOUNT3_MAX_RETRIES", "3")),
        user_agent=(
            os.environ.get("RECOUNT3_USER_AGENT")
            or "recount3-python/0.2 (+https://github.com/MoseleyBioinformaticsLab/recount3)"
        ),
        cache_dir=cache_dir,
        cache_disabled=os.environ.get("RECOUNT3_CACHE_DISABLE", "0") == "1",
        chunk_size=1024 * 1024,
    )
