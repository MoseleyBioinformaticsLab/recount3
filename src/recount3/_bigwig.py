"""Thin, safe wrapper around a pyBigWig handle.

This wrapper mirrors the original behavior and surfaces a small typed API.

Dependency:
  * pyBigWig (optional) - load lazily and raise ImportError with guidance.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Any, cast


try:  # Optional dependency; import lazily in _ensure_open for clarity.
    import pyBigWig  # type: ignore[import-not-found]
except Exception:  # pragma: no cover - optional dep
    pyBigWig = None  # type: ignore[assignment]


@dataclasses.dataclass(slots=True)
class BigWigFile:
    """Manage opening/closing a BigWig file with a small, typed API.

    Attributes:
      path: Filesystem path to a `.bw` file (must exist).
      mode: File mode for `pyBigWig.open`. Reading is the default ("r").
    """

    path: Path
    mode: str = "r"

    _bw: Any | None = dataclasses.field(default=None, init=False, repr=False)

    # ---- lifecycle -----------------------------------------------------

    def _ensure_open(self):  # -> pyBigWig.pyBigWig (annotated dynamically)
        """Open the file if needed and return the live pyBigWig handle.

        Returns:
          Live `pyBigWig` handle.

        Raises:
          FileNotFoundError: If ``path`` does not exist.
          ImportError: If ``pyBigWig`` is not installed.
          RuntimeError: If opening the file fails.
        """
        if self._bw is not None:
            return cast("pyBigWig.pyBigWig", self._bw)  # type: ignore

        if pyBigWig is None:
            raise ImportError(
                "pyBigWig is required to read BigWig files. "
                "Install with `pip install pyBigWig` or "
                "`conda install -c conda-forge -c bioconda pybigwig`."
            )
        if not self.path.exists():
            raise FileNotFoundError(str(self.path))

        bw = pyBigWig.open(str(self.path), self.mode)  # type: ignore[attr-defined]
        if bw is None:  # pyBigWig returns None on failure
            raise RuntimeError(f"Failed to open BigWig file: {self.path}")
        self._bw = bw
        return bw

    def close(self) -> None:
        """Close the underlying pyBigWig handle if open."""
        if self._bw is not None:
            try:
                self._bw.close()  # type: ignore[call-arg]
            finally:
                self._bw = None

    def __enter__(self):  # -> pyBigWig.pyBigWig
        return self._ensure_open()

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    # ---- convenience API -----------------------------------------------

    def is_open(self) -> bool:
        """Return True if the underlying handle is open."""
        return self._bw is not None

    def chroms(self, chrom: str | None = None):
        """Return chromosome lengths, or a single length."""
        bw = self._ensure_open()
        return bw.chroms(chrom) if chrom else bw.chroms()

    def header(self) -> dict[str, Any]:
        """Return the BigWig header metadata."""
        bw = self._ensure_open()
        return cast(dict[str, Any], bw.header())

    def values(self, chrom: str, start: int, end: int, *, numpy: bool | None = None):
        """Return per-base values over a half-open interval [start, end)."""
        bw = self._ensure_open()
        return bw.values(chrom, int(start), int(end), numpy=numpy)

    def stats(
        self,
        chrom: str,
        start: int | None = None,
        end: int | None = None,
        *,
        type: str = "mean",
        n_bins: int | None = None,
        exact: bool | None = None,
    ) -> list[float | None]:
        """Return summary statistic(s) over an interval or whole chromosome."""
        bw = self._ensure_open()
        kwargs: dict[str, Any] = {"type": type}
        if n_bins is not None:
            kwargs["nBins"] = int(n_bins)  # pyright: ignore[reportArgumentType]
        if exact is not None:
            kwargs["exact"] = bool(exact)  # pyright: ignore[reportArgumentType]
        if start is None or end is None:
            return bw.stats(chrom, **kwargs)
        return bw.stats(chrom, int(start), int(end), **kwargs)

    def intervals(self, chrom: str, start: int | None = None, end: int | None = None):
        """Return interval triples (start, end, value) overlapping a region."""
        bw = self._ensure_open()
        if start is None or end is None:
            return bw.intervals(chrom)
        return bw.intervals(chrom, int(start), int(end))
