# Copyright (c) 2026, Alexander Alsalihi, Robert M. Flight, Hunter N.B. Moseley.
# All rights reserved.
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
"""Read local BigWig coverage with :class:`recount3.BigWigFile`.

BigWig resources loaded through :meth:`~recount3.resource.R3Resource.load`
also return this wrapper. The optional ``pyBigWig`` dependency is imported
when opening a file; install it with ``pip install "recount3[bigwig]"``.

Example::

    >>> import recount3 as r3
    >>> reader = r3.BigWigFile("example.bw")
    >>> with reader:
    ...     lengths = reader.chroms()
    ...     mean = reader.stats("chr1", 0, 1000, exact=True)[0]
"""

from __future__ import annotations

from collections.abc import Mapping
import dataclasses
import operator
import pathlib
import types
from typing import Any, cast, TYPE_CHECKING

import numpy as np

from recount3 import _utils
from recount3 import errors

if TYPE_CHECKING:  # pragma: no cover
    import pyBigWig  # type: ignore[import-not-found]


def _native_string(value: str, name: str) -> str:
    """Return a plain string encodable as UTF-8 and free of NUL characters."""
    if not isinstance(value, str):
        raise TypeError(f"{name} must be a string")
    # Match the contents C reads; str(value) would call an overridden __str__.
    value = str.__str__(value)
    if "\0" in value:
        raise ValueError(f"{name} must not contain NUL characters")
    value.encode("utf-8")
    return value


def _coordinate(value: int, name: str) -> int:
    """Return an integer index in the unsigned 32-bit coordinate range."""
    result = operator.index(value)
    if not 0 <= result <= 2**32 - 1:
        raise ValueError(f"{name} must be between 0 and 2**32 - 1")
    return result


def _bounds(start: int | None, end: int | None) -> dict[str, int]:
    """Return validated coordinate keywords, omitting bounds set to None."""
    kwargs = {}
    if start is not None:
        kwargs["start"] = _coordinate(start, "start")
    if end is not None:
        kwargs["end"] = _coordinate(end, "end")
    return kwargs


@dataclasses.dataclass(slots=True)
class BigWigFile:
    """A lazily opened reader for local BigWig files.

    The first read or context entry opens a native ``pyBigWig`` handle, which
    is cached until :meth:`close`. Changes to ``path`` and ``mode`` take effect
    on the next open, after the cached handle closes.

    A context yields the native handle and closes it on exit. Call methods on
    the wrapper inside the context to use its validation and unsigned interval
    coordinates. Close through the wrapper: closing the native handle directly
    leaves a closed handle cached in the wrapper.

    Read methods raise :exc:`RuntimeError` in write mode before opening or
    querying. Use the native context-manager handle for writing.

    String arguments must encode as UTF-8 and contain no embedded NUL
    characters. Coordinates must be integers in ``0..2**32 - 1``; NumPy
    integer scalars are accepted.

    Attributes:
        path: Local path to a BigWig file. Must exist when opened for reading.
        mode: Mode passed to ``pyBigWig.open``; defaults to ``"r"``. A mode
          containing ``"w"`` creates or overwrites the file on each open.

    Note:
        Opening errors propagate from read methods and context entry:
        :exc:`FileNotFoundError` for a missing read path, :exc:`ImportError`
        for an unavailable ``pyBigWig`` dependency, :exc:`RuntimeError` for a
        native open failure, and :exc:`~recount3.errors.LoadError` for BigBed
        input.
    """

    path: pathlib.Path | str
    mode: str = "r"

    _bw: Any | None = dataclasses.field(default=None, init=False, repr=False)
    _bw_writable: bool = dataclasses.field(
        default=False, init=False, repr=False
    )

    def _ensure_open(self) -> pyBigWig.pyBigWig:
        """Return the cached handle, or open and cache a validated handle.

        Cache a reader only after format validation succeeds; close rejected
        handles. Write mode bypasses existence and format checks.

        Raises:
            TypeError: If ``path`` cannot be converted to a Path or ``mode``
              is not a string.
            ValueError: If the path or mode contains an embedded NUL.
            UnicodeEncodeError: If the path or mode cannot be encoded as UTF-8.
            FileNotFoundError: If the path does not exist when reading.
            ImportError: If ``pyBigWig`` is unavailable or fails to import.
            RuntimeError: If ``pyBigWig.open`` fails.
            recount3.errors.LoadError: If a reader opens a BigBed file.
        """
        if self._bw is not None:
            return cast("pyBigWig.pyBigWig", self._bw)

        mode = _native_string(self.mode, "mode")
        path = pathlib.Path(self.path)
        filename = _native_string(str(path), "path")
        writable = "w" in mode
        if not writable and not path.exists():
            raise FileNotFoundError(str(self.path))

        pybigwig = _utils.get_pybigwig_module()
        bw = pybigwig.open(filename, mode)
        if bw is None:
            raise RuntimeError(f"Failed to open BigWig file: {self.path}")

        try:
            if not writable and not bw.isBigWig():
                raise errors.LoadError(f"Not a BigWig file: {self.path}")
        except BaseException:
            bw.close()
            raise

        self._bw = bw
        self._bw_writable = writable
        return bw

    def _ensure_readable(self) -> pyBigWig.pyBigWig:
        """Return a read handle, rejecting writers before opening or querying.

        Opening a writer would truncate the file, and native ``values()`` can
        crash on a headerless writer. Check the cached handle's mode when
        present, since ``mode`` configures the next open.
        """
        writable = (
            self._bw_writable
            if self._bw is not None
            else "w" in _native_string(self.mode, "mode")
        )
        if writable:
            raise RuntimeError(
                "BigWig read methods cannot be used in write mode"
            )
        return self._ensure_open()

    def close(self) -> None:
        """Close and clear the cached handle, if any; repeated calls are safe."""
        if self._bw is None:
            return
        try:
            self._bw.close()
        finally:
            self._bw = None
            self._bw_writable = False

    def __enter__(self) -> pyBigWig.pyBigWig:
        """Return the cached or newly opened native ``pyBigWig`` handle."""
        return self._ensure_open()

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: types.TracebackType | None,
    ) -> None:
        """Close the cached handle on normal or exceptional context exit."""
        del exc_type, exc, tb
        self.close()

    def is_open(self) -> bool:
        """Return whether this wrapper holds a cached native handle.

        Closing the native handle directly does not update this result.
        """
        return self._bw is not None

    def chroms(
        self, chrom: str | None = None
    ) -> Mapping[str, int] | int | None:
        """Return chromosome lengths in bases.

        Args:
            chrom: Chromosome name, or None to request all chromosomes.

        Returns:
            A dictionary of names to lengths when ``chrom`` is None. Otherwise,
            the requested chromosome's length, or None if the name is unknown.

        Raises:
            TypeError: If ``chrom`` is neither a string nor None.
            ValueError: If ``chrom`` contains an embedded NUL.
            UnicodeEncodeError: If ``chrom`` cannot be encoded as UTF-8.
            RuntimeError: If the handle is in write mode.
        """
        if chrom is not None:
            chrom = _native_string(chrom, "chrom")
        bw = self._ensure_readable()
        return bw.chroms(chrom) if chrom is not None else bw.chroms()

    def header(self) -> dict[str, Any]:
        """Return BigWig header metadata as supplied by pyBigWig.

        Floating-point summaries may be truncated to integers, and extrema
        from the native writer may be inaccurate. Use :meth:`stats` with
        ``exact=True`` to compute statistics from the stored signal.

        Raises:
            RuntimeError: If the handle is in write mode.
        """
        bw = self._ensure_readable()
        return cast(dict[str, Any], bw.header())

    def values(
        self,
        chrom: str,
        start: int,
        end: int,
        *,
        numpy: bool | None = None,
    ) -> list[float] | np.ndarray:
        """Return per-base values over the half-open interval ``[start, end)``.

        Args:
            chrom: Chromosome name.
            start: 0-based integer start coordinate (inclusive).
            end: 0-based integer end coordinate (exclusive).
            numpy: True requests a NumPy array and requires ``pyBigWig`` built
              with NumPy support. False or None returns a list with either
              build.

        Returns:
            A list of floats, or an array with dtype ``float32``, in genomic
            order. Uncovered bases are NaN. A file with no entries returns an
            empty result, subject to the older NumPy-output bug noted below.

        Raises:
            TypeError: If ``chrom`` is not a string or a coordinate is not an
              integer, or NumPy output is requested from a build without
              NumPy support.
            ValueError: If ``chrom`` contains an embedded NUL, a coordinate is
              outside ``0..2**32 - 1``, or a list query exceeds ``2**31 - 1``
              bases.
            UnicodeEncodeError: If ``chrom`` cannot be encoded as UTF-8.
            RuntimeError: If the handle is in write mode or pyBigWig rejects
              the chromosome or interval.

        Notes:
            List output is limited to ``2**31 - 1`` bases by pyBigWig's signed
            loop counter. Both output forms allocate memory proportional to
            the query length; use smaller intervals to limit memory use.

            For NumPy reads of files with no entries, pyBigWig versions before
            0.3.25 return a scalar array containing an uninitialized value.
            Use 0.3.25 or newer for a correctly empty array.
        """
        chrom = _native_string(chrom, "chrom")
        start = _coordinate(start, "start")
        end = _coordinate(end, "end")
        if not numpy and end - start > 2**31 - 1:
            raise ValueError("List queries cannot exceed 2**31 - 1 bases")
        bw = self._ensure_readable()
        if numpy:
            return bw.values(chrom, start, end, numpy=True)
        return bw.values(chrom, start, end)

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
        """Return binned statistics over a half-open interval or chromosome.

        Args:
            chrom: Chromosome name.
            start: 0-based integer start coordinate (inclusive). Defaults to 0.
            end: 0-based integer end coordinate (exclusive). Defaults to the
              chromosome length independently of ``start``.
            type: One of ``"mean"``, ``"std"``, ``"min"``, ``"max"``,
              ``"coverage"``, or ``"sum"``. Native aliases ``"dev"`` for
              ``"std"`` and ``"cov"`` for ``"coverage"`` are also accepted.
            n_bins: Integer bin count in ``1..2**31 - 1``, forwarded as native
              ``nBins``. None requests one bin.
            exact: True computes from the stored intervals. False or None
              allows approximate zoom summaries when available.

        Returns:
            A list containing one float or None per bin, in genomic order.
            Missing-data bins produce None; native undefined results also
            appear as None.

        Raises:
            TypeError: If ``chrom`` or ``type`` is not a string, or a coordinate
              or bin count is not an integer.
            ValueError: If ``chrom`` or ``type`` contains an embedded NUL, a
              coordinate is outside ``0..2**32 - 1``, or ``n_bins`` is outside
              ``1..2**31 - 1``.
            UnicodeEncodeError: If ``chrom`` or ``type`` cannot be encoded as
              UTF-8.
            RuntimeError: If the handle is in write mode or pyBigWig rejects
              the chromosome, interval, or statistic type.

        Notes:
            Statistics use covered bases; ``coverage`` is the covered fraction
            of each bin. ``std`` uses sample standard deviation, returning zero
            for a single covered base when ``exact=True``.

            Use ``exact=True`` for sums: native zoom-summary bugs can
            substantially undercount them or return None for covered
            zero-signal regions. Keep ``n_bins`` at or below the number of
            queried bases; finer bins can give inconsistent native results
            even with ``exact=True``.
        """
        chrom = _native_string(chrom, "chrom")
        kwargs: dict[str, Any] = {
            "type": _native_string(type, "type"),
            **_bounds(start, end),
        }
        if n_bins is not None:
            n_bins = operator.index(n_bins)
            if not 1 <= n_bins <= 2**31 - 1:
                raise ValueError("n_bins must be between 1 and 2**31 - 1")
            kwargs["nBins"] = n_bins
        if exact is not None:
            kwargs["exact"] = bool(exact)

        bw = self._ensure_readable()
        return bw.stats(chrom, **kwargs)

    def intervals(
        self,
        chrom: str,
        start: int | None = None,
        end: int | None = None,
    ) -> tuple[tuple[int, int, float], ...] | None:
        """Return stored intervals overlapping a half-open query region.

        Args:
            chrom: Chromosome name.
            start: 0-based integer start coordinate (inclusive). Defaults to 0.
            end: 0-based integer end coordinate (exclusive). Defaults to the
              chromosome length independently of ``start``.

        Returns:
            A tuple of ``(start, end, value)`` tuples in coordinate order, or
            None when no intervals overlap. Endpoints are unsigned coordinates
            from the file and are not clipped to the query bounds.

        Raises:
            TypeError: If ``chrom`` is not a string or a coordinate is not an
              integer.
            ValueError: If ``chrom`` contains an embedded NUL or a coordinate
              is outside ``0..2**32 - 1``.
            UnicodeEncodeError: If ``chrom`` cannot be encoded as UTF-8.
            RuntimeError: If the handle is in write mode or pyBigWig rejects
              the chromosome or interval.
        """
        chrom = _native_string(chrom, "chrom")
        kwargs = _bounds(start, end)
        bw = self._ensure_readable()
        result = bw.intervals(chrom, **kwargs)
        if result and result[-1][1] < 0:
            return tuple(
                (start % 2**32, end % 2**32, value)
                for start, end, value in result
            )
        return result
