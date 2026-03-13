"""BigWig file access via a small, safe wrapper around ``pyBigWig``.

This module provides :class:`BigWigFile`, a lightweight convenience wrapper that:

* Defers importing the optional ``pyBigWig`` dependency until the first read.
* Opens the BigWig file lazily and caches the live handle.
* Ensures the handle is closed via an explicit :meth:`close` method and
  context-manager support.

The optional dependency is imported through :func:`recount3._utils.get_pybigwig_module`,
which standardizes error handling for missing/failed optional imports.

Typical usage:

  >>> from pathlib import Path
  >>> from recount3._bigwig import BigWigFile
  >>> with BigWigFile(Path("example.bw")) as bw:
  ...     lengths = bw.chroms()
  ...     mean = bw.stats("chr1", 0, 1000)[0]

Dependency:
  * pyBigWig (optional)
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Any, Mapping, TYPE_CHECKING, cast
from types import TracebackType

from recount3 import _utils

if TYPE_CHECKING:  # pragma: no cover
    import pyBigWig  # type: ignore[import-not-found]


@dataclasses.dataclass(slots=True)
class BigWigFile:
    """A lazily-opened BigWig reader with a small, typed API.

    Instances are cheap to construct and do not open the file until the first
    method call that requires a live ``pyBigWig`` handle. The handle is cached
    for subsequent calls and can be explicitly released with :meth:`close`.

    Attributes:
        path: Filesystem path to a BigWig file (typically ``.bw``). The file must
          exist when the handle is opened.
        mode: File mode passed to ``pyBigWig.open``. Reading is the default
          (``"r"``).
    """

    path: Path
    mode: str = "r"

    _bw: Any | None = dataclasses.field(default=None, init=False, repr=False)

    def _ensure_open(self) -> pyBigWig.pyBigWig:
        """Open the file if needed and return the live ``pyBigWig`` handle.

        This method is the single place that opens the underlying file. It is
        called internally by all public read methods.

        Returns:
            A live ``pyBigWig`` handle.

        Raises:
            FileNotFoundError: If ``path`` does not exist at open time.
            ImportError: If the optional ``pyBigWig`` dependency is missing or
              fails to import.
            RuntimeError: If ``pyBigWig.open`` fails (it returns ``None``).
        """
        if self._bw is not None:
            return cast("pyBigWig.pyBigWig", self._bw)

        if not self.path.exists():
            raise FileNotFoundError(str(self.path))

        pybigwig = _utils.get_pybigwig_module()
        bw = pybigwig.open(str(self.path), self.mode)
        if bw is None:  # pyBigWig returns None on failure
            raise RuntimeError(f"Failed to open BigWig file: {self.path}")

        self._bw = bw
        return bw

    def close(self) -> None:
        """Close the underlying ``pyBigWig`` handle if it is open.

        This method is idempotent.
        """
        if self._bw is None:
            return
        try:
            self._bw.close()
        finally:
            self._bw = None

    def __enter__(self) -> pyBigWig.pyBigWig:
        """Enter a context manager and return the live ``pyBigWig`` handle."""
        return self._ensure_open()

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        """Exit a context manager and close the file handle."""
        del exc_type, exc, tb
        self.close()

    def is_open(self) -> bool:
        """Return ``True`` if the underlying handle is currently open."""
        return self._bw is not None

    def chroms(self, chrom: str | None = None) -> Mapping[str, int] | int | None:
        """Return chromosome lengths, or a single chromosome length.

        Args:
            chrom: If provided, return only the length for this chromosome.

        Returns:
            If ``chrom`` is None, a mapping from chromosome name to length.
            Otherwise, the length for the requested chromosome. ``None`` may be
            returned if the chromosome is not present.
        """
        bw = self._ensure_open()
        return bw.chroms(chrom) if chrom else bw.chroms()

    def header(self) -> dict[str, Any]:
        """Return the BigWig header metadata."""
        bw = self._ensure_open()
        return cast(dict[str, Any], bw.header())

    def values(
        self,
        chrom: str,
        start: int,
        end: int,
        *,
        numpy: bool | None = None,
    ) -> list[float] | Any:
        """Return per-base values over a half-open interval ``[start, end)``.

        Args:
            chrom: Chromosome name.
            start: 0-based start coordinate (inclusive).
            end: 0-based end coordinate (exclusive).
            numpy: Forwarded to ``pyBigWig.values``. When True, ``pyBigWig`` may
              return a NumPy array depending on its configuration.

        Returns:
            Values returned by ``pyBigWig.values``. When ``numpy`` is not True, this is
            typically a ``list[float]``. When ``numpy`` is True, ``pyBigWig`` may
            return a NumPy array depending on its configuration.
        """
        bw = self._ensure_open()
        return bw.values(chrom, int(start), int(end), numpy=numpy)

    def stats(
        self,
        chrom: str,
        start: int | None = None,
        end: int | None = None,
        *,
        type: str = "mean",  # pylint: disable=redefined-builtin
        n_bins: int | None = None,
        exact: bool | None = None,
    ) -> list[float | None]:
        """Return summary statistic(s) over an interval or whole chromosome.

        Args:
            chrom: Chromosome name.
            start: 0-based start coordinate (inclusive). If omitted, stats are
              computed over the whole chromosome.
            end: 0-based end coordinate (exclusive). If omitted, stats are
              computed over the whole chromosome.
            type: Statistic name understood by ``pyBigWig.stats`` (for example,
              ``"mean"``, ``"min"``, ``"max"``, ``"coverage"``).
            n_bins: If provided, request binned stats via the ``nBins`` argument.
            exact: If provided, forward to the ``exact`` argument.

        Returns:
            A list of statistic values as returned by ``pyBigWig.stats``. Values
            may be ``None`` for missing data regions.
        """
        bw = self._ensure_open()
        kwargs: dict[str, Any] = {"type": type}
        if n_bins is not None:
            kwargs["nBins"] = int(n_bins)
        if exact is not None:
            kwargs["exact"] = bool(exact)

        if start is None or end is None:
            return bw.stats(chrom, **kwargs)

        return bw.stats(chrom, int(start), int(end), **kwargs)

    def intervals(
        self,
        chrom: str,
        start: int | None = None,
        end: int | None = None,
    ) -> list[tuple[int, int, float]] | None | Any:
        """Return (start, end, value) intervals overlapping a region.

        Args:
            chrom: Chromosome name.
            start: 0-based start coordinate (inclusive). If omitted, intervals
              for the entire chromosome may be returned.
            end: 0-based end coordinate (exclusive). If omitted, intervals for
              the entire chromosome may be returned.

        Returns:
            Intervals returned by ``pyBigWig.intervals``. This is often a list of
            ``(start, end, value)`` tuples. ``None`` may be returned when no intervals
            overlap the requested region.
        """
        bw = self._ensure_open()
        if start is None or end is None:
            return bw.intervals(chrom)
        return bw.intervals(chrom, int(start), int(end))
