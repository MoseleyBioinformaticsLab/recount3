"""Filesystem helpers (internal).

These helpers are intentionally small and side-effect free beyond interacting
with the filesystem. They preserve the behavior of the original script.
"""

from __future__ import annotations

import errno
import os
import shutil
from pathlib import Path


def _ensure_dir(path: str | Path) -> None:
    """Ensure a directory exists, creating parents as needed.

    Args:
      path: Directory path to create.

    Raises:
      NotADirectoryError: If the path exists as a non-directory.
    """
    p = Path(path)
    if not p.exists():
        p.mkdir(parents=True, exist_ok=True)
    elif not p.is_dir():
        raise NotADirectoryError(str(p))


def _hardlink_or_copy(src: Path, dst: Path) -> None:
    """Materialize by hardlink, falling back to copy on cross-device/perm errors.

    Mirrors the original logic: attempt to hardlink, and on known errno values
    (EXDEV, EPERM, EACCES, EMLINK) copy with metadata preservation.
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
        else:  # unexpected
            raise


def _atomic_replace(src_tmp: Path, final_path: Path) -> None:
    """Atomically replace the final path with the temporary file."""
    _ensure_dir(final_path.parent)
    os.replace(src_tmp, final_path)
