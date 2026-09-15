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
  * ``RECOUNT3_URL``: base URL of the recount3 mirror
    (default: ``http://duffel.rail.bio/recount3/``).
  * ``RECOUNT3_CACHE_DIR``: directory for the on-disk file cache
    (default: ``~/.cache/recount3/files`` for ``filesystem``; R's
    ``tools::R_user_dir("recount3", "cache")`` for ``pybiocfilecache``).
  * ``RECOUNT3_CACHE_BACKEND``: ``filesystem`` (default) or optional
    ``pybiocfilecache`` registry.
  * ``RECOUNT3_CACHE_DISABLE``: set to ``"1"`` to disable caching entirely.
  * ``RECOUNT3_HTTP_TIMEOUT``: HTTP request timeout in seconds (default: 60).
  * ``RECOUNT3_MAX_RETRIES``: maximum retry attempts for transient errors
    (default: 3).
  * ``RECOUNT3_INSECURE_SSL``: set to ``"1"`` to skip TLS certificate
    verification. This only affects ``https://`` base URLs; it is a no-op for
    the default ``http://`` Duffel mirror (see "Mirrors" below).
  * ``RECOUNT3_USER_AGENT``: custom ``User-Agent`` header string.
  * ``RECOUNT3_CHUNK_SIZE``: streaming chunk size in bytes
    (default: 1048576, i.e. 1 MiB).

Mirrors:
  recount3 publishes the **same relative file layout** on several
  interchangeable public mirrors; ``RECOUNT3_URL`` / ``base_url`` works with
  any of them unchanged:

  * Duffel load balancer (default): ``http://duffel.rail.bio/recount3/``
  * AWS Open Data:
    ``https://recount-opendata.s3.amazonaws.com/recount3/release/``
  * JHU IDIES (Dataverse): ``https://data.idies.jhu.edu/recount3/data/``

  The package is coupled to recount3's layout convention rather than to any one
  host, so switching mirrors is purely a ``base_url`` change. TLS settings apply
  only to ``https://`` endpoints: the default Duffel mirror is plain ``http``
  (no TLS, so ``insecure_ssl`` is irrelevant), the AWS and JHU mirrors are
  ``https`` with valid certificates (no flag needed), and
  ``RECOUNT3_INSECURE_SSL`` / ``--insecure-ssl`` is meaningful only for an
  ``https`` endpoint presenting an untrusted or self-signed certificate.

Typical usage example::

    import dataclasses
    from pathlib import Path
    import recount3 as r3

    # Read the current cache directory (creates it if absent):
    cache_dir = r3.recount3_cache()

    # Use a custom cache location for this session (override one field of
    # the environment-derived defaults; Config is immutable):
    cfg = dataclasses.replace(
        r3.default_config(), cache_dir=Path("/scratch/recount3_cache")
    )
    custom_cache_dir = r3.recount3_cache(cfg)

    # Remove cached files matching a pattern (dry run first):
    to_delete = r3.recount3_cache_rm(dry_run=True)
    r3.recount3_cache_rm(predicate=lambda p: "sra" in str(p))
"""

from __future__ import annotations

import os
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from recount3.version import __version__

from recount3._utils import (
    _ensure_dir,
    _cache_internal,
    _biocfilecache_paths,
    _remove_biocfilecache_files,
)

_DEFAULT_CHUNK_SIZE: int = 1024 * 1024  # 1 MiB


@dataclass(frozen=True, slots=True)
class Config:
    """Immutable configuration bag.

    Attributes:
      base_url: Base URL for the recount3 mirror (ends with a trailing slash).
          Defaults to the Duffel load balancer; the AWS Open Data and JHU IDIES
          mirrors serve the same layout (see the module docstring "Mirrors").
      timeout: Network timeout in seconds.
      insecure_ssl: True to disable TLS certificate verification (not
          recommended). Only affects ``https://`` base URLs; a no-op for the
          default ``http://`` Duffel mirror.
      max_retries: Max HTTP retry attempts for transient errors.
      user_agent: Custom HTTP User-Agent.
      cache_dir: Cache directory for downloaded files.
      cache_disabled: If True, disable cache behavior globally.
      cache_backend: "filesystem" or opt-in "pybiocfilecache" registry.
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
    cache_backend: str = "filesystem"

    def __post_init__(self) -> None:
        """Reject unsupported cache backends before performing any I/O."""
        if self.cache_backend not in ("filesystem", "pybiocfilecache"):
            raise ValueError(f"Unknown cache backend: {self.cache_backend!r}")


def _default_cache_dir(backend: str) -> Path:
    """Choose the native directory or mirror R's R_user_dir cache rules."""
    if backend != "pybiocfilecache":
        return Path.home() / ".cache" / "recount3" / "files"
    root = os.environ.get("R_USER_CACHE_DIR") or os.environ.get(
        "XDG_CACHE_HOME"
    )
    if root:
        base = Path(root)
    elif sys.platform == "win32":
        base = Path(os.environ.get("LOCALAPPDATA", "")) / "R" / "cache"
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches" / "org.R-project.R"
    else:
        base = Path.home() / ".cache"
    return base / "R" / "recount3"


def default_config(*, cache_backend: str | None = None) -> Config:
    """Return configuration constructed from environment variables.

    Args:
      cache_backend: Override the environment's backend before selecting its
        default directory. An explicit ``RECOUNT3_CACHE_DIR`` still wins.

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
    if cache_backend is None:
        cache_backend = os.environ.get("RECOUNT3_CACHE_BACKEND", "filesystem")
    cache_dir = Path(
        os.environ.get("RECOUNT3_CACHE_DIR", _default_cache_dir(cache_backend))
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
        cache_backend=cache_backend,
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
      A sorted list of payload paths, excluding database internals. The
      optional backend includes registered external or missing paths as well
      as native files. An absent cache directory yields an empty list.
    """
    cfg = config or default_config()
    root = cfg.cache_dir

    if not root.exists() or not root.is_dir():
        return []

    glob_pattern = pattern if pattern is not None else "*"
    files: list[Path] = []

    for path in root.rglob(glob_pattern):
        if path.is_file() and not _cache_internal(path, root):
            files.append(path)
    if cfg.cache_backend == "pybiocfilecache":
        for path in _biocfilecache_paths(root):
            if not _cache_internal(path, root) and path.match(glob_pattern):
                files.append(path)
        files = list({path.resolve() for path in files})
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
      A sorted list of :class:`pathlib.Path` objects that were removed from
      the cache (or would be removed when ``dry_run`` is True). With the
      optional backend, a registered payload outside the cache directory has
      its registry row dropped while the external file is left in place,
      matching R's ``action="asis"`` semantics.

    Raises:
      OSError: If filesystem operations fail during deletion.

    Examples:
        List everything that would be removed, without deleting::

            to_delete = r3.recount3_cache_rm(dry_run=True)

        Remove all cached files (empty the cache)::

            r3.recount3_cache_rm()

        Remove only files related to the ``"sra"`` data source::

            r3.recount3_cache_rm(predicate=lambda p: "sra" in str(p))
    """
    cfg = config or default_config()
    root = cfg.cache_dir

    if not root.exists() or not root.is_dir():
        return []

    def _select(path: Path) -> bool:
        if predicate is None:
            return True
        return predicate(path)

    candidates = [p for p in recount3_cache_files(cfg) if _select(p)]

    if dry_run:
        return candidates

    if cfg.cache_backend == "pybiocfilecache":
        _remove_biocfilecache_files(root, candidates)
    else:
        for path in candidates:
            path.unlink()

    return candidates
