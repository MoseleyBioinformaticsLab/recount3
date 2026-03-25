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
"""BigWig file access via a small wrapper around ``pyBigWig``.

This is an internal module. For BigWig access via the public API, use
:meth:`~recount3.resource.R3Resource.load` on a BigWig
:class:`~recount3.resource.R3Resource`.

The optional dependency is imported through
:func:`~recount3._utils.get_pybigwig_module`.

Typical usage::

  >>> from pathlib import Path
  >>> from recount3._bigwig import BigWigFile
  >>> with BigWigFile(Path("example.bw")) as bw:
  ...     lengths = bw.chroms()
  ...     mean = bw.stats("chr1", 0, 1000)[0]

Note:
    Requires the optional ``pyBigWig`` package. Install with
    ``pip install pyBigWig``. An :exc:`ImportError` is raised on first use
    if the package is not available. ``pyBigWig`` can be difficult to
    install on non-Linux systems.
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

    def chroms(
        self, chrom: str | None = None
    ) -> Mapping[str, int] | int | None:
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
