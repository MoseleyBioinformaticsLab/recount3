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

from __future__ import annotations

import math
import pathlib
from typing import Any
from unittest import mock

import numpy as np
import pytest

import recount3
from recount3 import _utils
from recount3 import errors
from recount3._bigwig import BigWigFile

_DATA_DIR = pathlib.Path(__file__).parent / "data"
_BW_PATH = (
    _DATA_DIR
    / "recount3_mirror"
    / "recount3"
    / "human"
    / "data_sources"
    / "sra"
    / "base_sums"
    / "15"
    / "SRP009615"
    / "77"
    / "sra.base_sums.SRP009615_SRR387777.ALL.bw"
)

_BW_CHROM = "chr1"
_BW_CHROM_LEN = 2_000_000

_SENTINEL: Any = object()


class _StringSubclass(str):
    """Python overrides must not disguise the contents seen by C parsers."""

    def __contains__(self, item: object) -> bool:
        return False

    def encode(self, *args: Any, **kwargs: Any) -> bytes:
        return b"disguised"

    def __str__(self) -> str:
        return "disguised"


def _mock_pybigwig(
    open_return: Any = _SENTINEL,
) -> tuple[mock.MagicMock, mock.MagicMock]:
    """Return ``(mock_module, mock_handle)`` for patching get_pybigwig_module.

    Args:
        open_return: Value returned by ``mock_module.open()``. Defaults to a
            fresh :class:`unittest.mock.MagicMock` acting as the file handle.

    Returns:
        A 2-tuple ``(mock_module, mock_handle)``.
    """
    handle = mock.MagicMock()
    handle.isBigWig.return_value = True
    handle.intervals.return_value = None
    module = mock.MagicMock()
    module.open.return_value = (
        handle if open_return is _SENTINEL else open_return
    )
    return module, handle


def test_construction_defaults() -> None:
    """BigWigFile stores path and defaults mode to 'r' without opening."""
    bw = BigWigFile(path=_BW_PATH)
    assert bw.path == _BW_PATH
    assert bw.mode == "r"
    assert not bw.is_open()


def test_public_export_is_the_wrapper() -> None:
    assert recount3.BigWigFile is BigWigFile


def test_construction_custom_mode() -> None:
    """BigWigFile stores a caller-supplied mode string."""
    bw = BigWigFile(path=_BW_PATH, mode="w")
    assert bw.mode == "w"


def test_ensure_open_returns_cached_handle() -> None:
    """Second _ensure_open call returns the same handle without re-opening."""
    mock_module, _ = _mock_pybigwig()
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        h1 = bw._ensure_open()
        h2 = bw._ensure_open()
    bw.close()
    assert h1 is h2
    mock_module.open.assert_called_once_with(str(_BW_PATH), "r")


def test_ensure_open_missing_path_raises_file_not_found(
    tmp_path: pathlib.Path,
) -> None:
    """FileNotFoundError is raised when the path does not exist on disk."""
    bw = BigWigFile(path=tmp_path / "nonexistent.bw")
    with pytest.raises(FileNotFoundError):
        bw._ensure_open()


def test_ensure_open_pybigwig_returns_none_raises_runtime_error() -> None:
    """RuntimeError is raised when pyBigWig.open() returns None."""
    mock_module, _ = _mock_pybigwig(open_return=None)
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        with pytest.raises(RuntimeError, match="Failed to open BigWig file"):
            bw._ensure_open()


def test_ensure_open_import_error_propagates() -> None:
    """ImportError from get_pybigwig_module propagates unchanged."""
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils,
        "get_pybigwig_module",
        side_effect=ImportError("pyBigWig not found"),
    ):
        with pytest.raises(ImportError, match="pyBigWig not found"):
            bw._ensure_open()


def test_close_idempotent_on_unopened_instance() -> None:
    """close() on a never-opened BigWigFile is a safe no-op."""
    bw = BigWigFile(path=_BW_PATH)
    bw.close()
    assert not bw.is_open()


def test_close_calls_bw_close_and_clears_reference() -> None:
    """close() delegates to _bw.close() and sets _bw to None."""
    mock_module, mock_handle = _mock_pybigwig()
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        bw._ensure_open()
    bw.close()
    mock_handle.close.assert_called_once()
    assert not bw.is_open()


def test_close_clears_reference_when_bw_close_raises() -> None:
    """_bw is set to None in the finally block even if _bw.close() raises."""
    mock_module, mock_handle = _mock_pybigwig()
    mock_handle.close.side_effect = OSError("handle gone")
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        bw._ensure_open()
    with pytest.raises(OSError, match="handle gone"):
        bw.close()
    assert not bw.is_open()


def test_is_open_false_before_first_read() -> None:
    """is_open() returns False on a freshly constructed instance."""
    assert BigWigFile(path=_BW_PATH).is_open() is False


def test_is_open_true_after_ensure_open() -> None:
    """is_open() returns True once _ensure_open() has succeeded."""
    mock_module, _ = _mock_pybigwig()
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        bw._ensure_open()
        assert bw.is_open() is True
    bw.close()


def test_is_open_false_after_close() -> None:
    """is_open() returns False after close()."""
    mock_module, _ = _mock_pybigwig()
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        bw._ensure_open()
    bw.close()
    assert bw.is_open() is False


def test_context_manager_returns_handle_and_closes() -> None:
    """Context manager yields the pyBigWig handle and closes on exit."""
    mock_module, mock_handle = _mock_pybigwig()
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        with bw as handle:
            assert handle is mock_handle
            assert bw.is_open()
    assert not bw.is_open()
    mock_handle.close.assert_called_once()


def test_context_manager_closes_on_exception() -> None:
    """Context manager closes the handle even when the body raises."""
    mock_module, mock_handle = _mock_pybigwig()
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        with pytest.raises(ValueError):
            with bw:
                raise ValueError("body error")
    assert not bw.is_open()
    mock_handle.close.assert_called_once()


def test_exit_passes_exception_info_and_closes() -> None:
    """__exit__ receives exception arguments, closes the handle, and returns."""
    mock_module, _ = _mock_pybigwig()
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        bw._ensure_open()
        bw.__exit__(TypeError, TypeError("e"), None)
    assert not bw.is_open()


@pytest.mark.requires_pybigwig
def test_chroms_all_returns_full_mapping() -> None:
    """chroms() returns a dict mapping every chromosome to its length."""
    bw = BigWigFile(path=_BW_PATH)
    result = bw.chroms()
    assert isinstance(result, dict)
    assert result.get(_BW_CHROM) == _BW_CHROM_LEN
    bw.close()


@pytest.mark.requires_pybigwig
def test_chroms_single_returns_integer_length() -> None:
    """chroms(chrom) returns the integer length for the named chromosome."""
    bw = BigWigFile(path=_BW_PATH)
    result = bw.chroms(_BW_CHROM)
    assert result == _BW_CHROM_LEN
    bw.close()


@pytest.mark.requires_pybigwig
def test_header_returns_dict_with_standard_fields() -> None:
    """header() returns a dict that includes standard BigWig header keys."""
    bw = BigWigFile(path=_BW_PATH)
    h = bw.header()
    assert isinstance(h, dict)
    assert "version" in h
    assert "nBasesCovered" in h
    bw.close()


@pytest.mark.requires_pybigwig
def test_values_returns_list_of_per_base_values() -> None:
    """values() returns exactly one float per base in the requested interval."""
    bw = BigWigFile(path=_BW_PATH)
    result = bw.values(_BW_CHROM, 0, 10)
    assert isinstance(result, list)
    assert len(result) == 10
    bw.close()


def test_values_numpy_request_forwarded() -> None:
    """True requests native NumPy output."""
    mock_module, mock_handle = _mock_pybigwig()
    mock_handle.values.return_value = [1.0, 2.0]
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        bw.values(_BW_CHROM, 0, 2, numpy=True)
    mock_handle.values.assert_called_once_with(_BW_CHROM, 0, 2, numpy=True)


@pytest.mark.requires_pybigwig
def test_stats_whole_chrom_no_coords() -> None:
    """stats(chrom) with no start/end summarises the whole chromosome."""
    bw = BigWigFile(path=_BW_PATH)
    result = bw.stats(_BW_CHROM)
    assert isinstance(result, list)
    assert len(result) == 1
    bw.close()


@pytest.mark.requires_pybigwig
def test_stats_with_start_and_end_coords() -> None:
    """stats(chrom, start, end) returns stats for the given interval."""
    bw = BigWigFile(path=_BW_PATH)
    result = bw.stats(_BW_CHROM, 0, 10_000)
    assert isinstance(result, list)
    assert len(result) == 1
    bw.close()


def test_stats_n_bins_forwarded_as_nbins() -> None:
    """Providing n_bins adds the 'nBins' keyword to the pyBigWig.stats call."""
    mock_module, mock_handle = _mock_pybigwig()
    mock_handle.stats.return_value = [1.0, 2.0, 3.0]
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        bw.stats(_BW_CHROM, 0, 3_000, n_bins=3)
    assert mock_handle.stats.call_args.kwargs.get("nBins") == 3


def test_stats_exact_forwarded_to_pybigwig() -> None:
    """Providing exact=True adds 'exact' to the pyBigWig.stats call."""
    mock_module, mock_handle = _mock_pybigwig()
    mock_handle.stats.return_value = [0.5]
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        bw.stats(_BW_CHROM, 0, 1_000, exact=True)
    assert mock_handle.stats.call_args.kwargs.get("exact") is True


def test_stats_n_bins_and_exact_combined() -> None:
    """n_bins and exact can be combined; both appear in the underlying call."""
    mock_module, mock_handle = _mock_pybigwig()
    mock_handle.stats.return_value = [1.0, 2.0]
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        result = bw.stats(
            _BW_CHROM, 0, 2_000, type="sum", n_bins=2, exact=False
        )
    assert result is mock_handle.stats.return_value
    mock_handle.stats.assert_called_once_with(
        _BW_CHROM, start=0, end=2_000, type="sum", nBins=2, exact=False
    )


def test_intervals_whole_chrom_omits_coords_from_call() -> None:
    """intervals(chrom) calls pyBigWig.intervals without start/end coords."""
    mock_module, mock_handle = _mock_pybigwig()
    mock_handle.intervals.return_value = ((0, 1000, 0.0),)
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        result = bw.intervals(_BW_CHROM)
    mock_handle.intervals.assert_called_once_with(_BW_CHROM)
    assert result == ((0, 1000, 0.0),)


@pytest.mark.requires_pybigwig
def test_intervals_with_coords_returns_region_tuples() -> None:
    """intervals(chrom, start, end) returns (start, end, value) tuples."""
    bw = BigWigFile(path=_BW_PATH)
    result = bw.intervals(_BW_CHROM, 0, 1_000)
    assert result is not None
    assert isinstance(result, tuple)
    assert result[0] == (0, 1000, 0.0)
    bw.close()


@pytest.fixture
def tiny_bw_path(tmp_path: pathlib.Path) -> pathlib.Path:
    """A 16-base chromosome with gaps, zeros, and signed coverage."""
    pybigwig = pytest.importorskip("pyBigWig")
    path = tmp_path / "tiny.bw"
    with pybigwig.open(str(path), "w") as handle:
        handle.addHeader([("chr1", 16), ("empty", 8)])
        handle.addEntries(
            ["chr1"] * 4,
            [0, 4, 8, 12],
            ends=[2, 6, 10, 14],
            values=[1.0, 0.0, 3.0, -2.0],
        )
    return path


@pytest.mark.requires_pybigwig
@pytest.mark.parametrize("method", ["stats", "intervals"])
@pytest.mark.parametrize(
    "bounds, expected_bounds",
    [
        ({}, (0, 16)),
        ({"start": 8}, (8, 16)),
        ({"end": 6}, (0, 6)),
        ({"start": 4, "end": 10}, (4, 10)),
    ],
)
def test_optional_bounds_match_native(
    tiny_bw_path: pathlib.Path,
    method: str,
    bounds: dict[str, int],
    expected_bounds: tuple[int, int],
) -> None:
    """Omitting one coordinate must not discard the coordinate supplied."""
    reader = BigWigFile(tiny_bw_path)
    with reader as native:
        expected = getattr(native, method)("chr1", *expected_bounds)
        assert getattr(reader, method)("chr1", **bounds) == expected


@pytest.mark.parametrize("numpy_flag", [None, False])
def test_values_without_numpy_support(numpy_flag: bool | None) -> None:
    """Non-NumPy builds expose a positional-only values() method."""
    module, handle = _mock_pybigwig()
    handle.values.side_effect = lambda chrom, start, end, /: [1.0]
    reader = BigWigFile(_BW_PATH)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        try:
            assert reader.values("chr1", 0, 1, numpy=numpy_flag) == [1.0]
        finally:
            reader.close()
    handle.values.assert_called_once_with("chr1", 0, 1)


def test_chroms_empty_string_is_a_chromosome_query() -> None:
    """Only None requests all chromosomes; an empty name is still a name."""
    module, handle = _mock_pybigwig()
    reader = BigWigFile(_BW_PATH)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        try:
            assert reader.chroms("") is handle.chroms.return_value
        finally:
            reader.close()
    handle.chroms.assert_called_once_with("")


@pytest.mark.parametrize("method", ["values", "stats", "intervals"])
@pytest.mark.parametrize("coordinate", [0.9, "1", float("nan")])
@pytest.mark.parametrize("bound", ["start", "end"])
def test_coordinates_must_be_integers(
    method: str,
    coordinate: Any,
    bound: str,
) -> None:
    """Reject fractional/string coordinates instead of silently converting."""
    module, handle = _mock_pybigwig()
    reader = BigWigFile(_BW_PATH)
    bounds = {"start": 0, "end": 10, bound: coordinate}
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        try:
            with pytest.raises(TypeError):
                getattr(reader, method)("chr1", **bounds)
        finally:
            reader.close()
    getattr(handle, method).assert_not_called()


@pytest.mark.parametrize("method", ["values", "stats", "intervals"])
@pytest.mark.parametrize("coordinate", [-1, 2**32, 2**32 + 1])
@pytest.mark.parametrize("bound", ["start", "end"])
def test_coordinates_cannot_wrap_in_native_code(
    method: str,
    coordinate: int,
    bound: str,
) -> None:
    """Validate the unsigned 32-bit range before passing coordinates to C."""
    module, handle = _mock_pybigwig()
    reader = BigWigFile(_BW_PATH)
    bounds = {"start": 0, "end": 10, bound: coordinate}
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        try:
            with pytest.raises(ValueError, match=bound):
                getattr(reader, method)("chr1", **bounds)
        finally:
            reader.close()
    getattr(handle, method).assert_not_called()


@pytest.mark.parametrize("n_bins", [0, -1, 2**31])
def test_n_bins_must_be_positive_c_int(n_bins: int) -> None:
    """Invalid bin counts must not reach native allocation code."""
    module, handle = _mock_pybigwig()
    reader = BigWigFile(_BW_PATH)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        try:
            with pytest.raises(ValueError, match="n_bins"):
                reader.stats("chr1", 0, 10, n_bins=n_bins)
        finally:
            reader.close()
    handle.stats.assert_not_called()


@pytest.mark.parametrize("n_bins", [1.9, "2"])
def test_n_bins_must_be_an_integer(n_bins: Any) -> None:
    reader = BigWigFile(_BW_PATH)
    module, handle = _mock_pybigwig()
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        try:
            with pytest.raises(TypeError):
                reader.stats("chr1", 0, 10, n_bins=n_bins)
        finally:
            reader.close()
    handle.stats.assert_not_called()


@pytest.mark.parametrize("method", ["values", "stats", "intervals"])
def test_largest_coordinate_is_forwarded_without_allocation(
    method: str,
) -> None:
    """Exercise valid uint32 endpoints using mocks, avoiding a huge query."""
    module, handle = _mock_pybigwig()
    handle.intervals.return_value = ((2**32 - 2, 2**32 - 1, 1.0),)
    reader = BigWigFile(_BW_PATH)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        try:
            result = getattr(reader, method)("chr1", 2**32 - 2, 2**32 - 1)
        finally:
            reader.close()
    assert result is getattr(handle, method).return_value
    if method == "values":
        handle.values.assert_called_once_with("chr1", 2**32 - 2, 2**32 - 1)
    else:
        expected: dict[str, Any] = {"start": 2**32 - 2, "end": 2**32 - 1}
        if method == "stats":
            expected["type"] = "mean"
        getattr(handle, method).assert_called_once_with("chr1", **expected)


def test_largest_bin_count_is_forwarded_without_allocation() -> None:
    module, handle = _mock_pybigwig()
    reader = BigWigFile(_BW_PATH)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        try:
            reader.stats("chr1", n_bins=2**31 - 1)
        finally:
            reader.close()
    handle.stats.assert_called_once_with("chr1", type="mean", nBins=2**31 - 1)


@pytest.mark.parametrize("span", [2**31, 2**32 - 1])
@pytest.mark.parametrize("numpy_flag", [None, False])
def test_values_list_span_cannot_overflow_native_loop(
    span: int, numpy_flag: bool | None
) -> None:
    """Exercise dangerous sizes with mocks, never allocating their output."""
    module, handle = _mock_pybigwig()
    reader = BigWigFile(_BW_PATH)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        try:
            with pytest.raises(ValueError, match=r"2\*\*31 - 1"):
                reader.values("chr1", 0, span, numpy=numpy_flag)
        finally:
            reader.close()
    handle.values.assert_not_called()
    module.open.assert_not_called()


@pytest.mark.parametrize(
    "span, numpy_flag", [(2**31 - 1, False), (2**31, True)]
)
def test_large_value_spans_forwarded_without_allocation(
    span: int, numpy_flag: bool
) -> None:
    module, handle = _mock_pybigwig()
    reader = BigWigFile(_BW_PATH)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        try:
            assert (
                reader.values("chr1", 0, span, numpy=numpy_flag)
                is handle.values.return_value
            )
        finally:
            reader.close()
    expected = {"numpy": True} if numpy_flag else {}
    handle.values.assert_called_once_with("chr1", 0, span, **expected)


@pytest.mark.requires_pybigwig
def test_numpy_integer_coordinates_and_bins(tiny_bw_path: pathlib.Path) -> None:
    reader = BigWigFile(tiny_bw_path)
    with reader:
        assert reader.values("chr1", np.int64(0), np.uint32(2)) == [1.0, 1.0]
        assert reader.stats(
            "chr1", np.int64(0), np.int64(2), n_bins=np.int64(2), exact=True
        ) == [1.0, 1.0]
        assert reader.intervals("chr1", np.int64(0), np.int64(2)) == (
            (0, 2, 1.0),
        )


@pytest.mark.requires_pybigwig
def test_string_path_and_reopen(tiny_bw_path: pathlib.Path) -> None:
    reader = BigWigFile(str(tiny_bw_path))
    with reader as first:
        assert reader.chroms("chr1") == 16
    assert not reader.is_open()
    with reader as second:
        assert second is not first
        assert reader.values("chr1", 0, 2) == [1.0, 1.0]


@pytest.mark.requires_pybigwig
@pytest.mark.parametrize("mode", ["w", "wb"])
def test_write_mode_can_create_a_new_file(
    tmp_path: pathlib.Path, mode: str
) -> None:
    path = tmp_path / "created.bw"
    writer = BigWigFile(path, mode=mode)
    assert not path.exists()
    with writer as native:
        native.addHeader([("chr1", 2)])
        native.addEntries(["chr1"], [0], ends=[2], values=[2.0])
    assert not writer.is_open()
    reader = BigWigFile(path)
    with reader:
        assert reader.values("chr1", 0, 2) == [2.0, 2.0]


@pytest.mark.parametrize(
    "method", ["chroms", "header", "values", "stats", "intervals"]
)
@pytest.mark.parametrize("mode", ["w", "wb"])
def test_read_methods_reject_writers_without_opening(
    tmp_path: pathlib.Path, method: str, mode: str
) -> None:
    """A mistaken read must not truncate a file or reach unsafe native code."""
    path = tmp_path / "existing.bw"
    path.write_bytes(b"do not overwrite")
    reader = BigWigFile(path, mode=mode)
    module, handle = _mock_pybigwig()
    args = () if method in {"chroms", "header"} else ("chr1", 0, 1)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        with pytest.raises(RuntimeError, match="writ"):
            getattr(reader, method)(*args)
    module.open.assert_not_called()
    getattr(handle, method).assert_not_called()
    assert path.read_bytes() == b"do not overwrite"
    assert not reader.is_open()


@pytest.mark.parametrize(
    "method", ["chroms", "header", "values", "stats", "intervals"]
)
def test_read_methods_reject_an_open_writer(
    method: str, tmp_path: pathlib.Path
) -> None:
    """Cached handles retain their mode even if the next-open mode changes."""
    module, handle = _mock_pybigwig()
    writer = BigWigFile(tmp_path / "writer.bw", mode="w")
    args = () if method in {"chroms", "header"} else ("chr1", 0, 1)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        with writer:
            writer.mode = "r"
            with pytest.raises(RuntimeError, match="writ"):
                getattr(writer, method)(*args)
            getattr(handle, method).assert_not_called()
            assert writer.is_open()
    handle.close.assert_called_once_with()


@pytest.mark.requires_pybigwig
def test_reading_cached_handle_uses_its_original_mode(
    tiny_bw_path: pathlib.Path,
) -> None:
    reader = BigWigFile(tiny_bw_path)
    with reader:
        reader.mode = "w"
        assert reader.values("chr1", 0, 2) == [1.0, 1.0]


@pytest.mark.requires_pybigwig
def test_read_guard_resets_on_reopen(tmp_path: pathlib.Path) -> None:
    reader = BigWigFile(tmp_path / "reopen.bw", mode="w")
    with reader as native:
        native.addHeader([("chr1", 2)])
        native.addEntries(["chr1"], [0], ends=[2], values=[2.0])
    reader.mode = "r"
    with reader:
        assert reader.values("chr1", 0, 2) == [2.0, 2.0]


def test_open_error_is_not_cached_and_can_be_retried() -> None:
    module, handle = _mock_pybigwig()
    module.open.side_effect = [RuntimeError("bad file"), handle]
    reader = BigWigFile(_BW_PATH)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        with pytest.raises(RuntimeError, match="bad file"):
            reader._ensure_open()
        assert not reader.is_open()
        assert reader._ensure_open() is handle
        reader.close()
    assert module.open.call_count == 2


def test_bigbed_is_rejected_and_handle_is_closed() -> None:
    module, handle = _mock_pybigwig()
    handle.isBigWig.return_value = False
    reader = BigWigFile(_BW_PATH)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        with pytest.raises(errors.LoadError, match="Not a BigWig file"):
            reader._ensure_open()
    assert not reader.is_open()
    handle.close.assert_called_once_with()


def test_type_check_failure_closes_unowned_handle() -> None:
    module, handle = _mock_pybigwig()
    handle.isBigWig.side_effect = RuntimeError("unreadable type")
    reader = BigWigFile(_BW_PATH)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        with pytest.raises(RuntimeError, match="unreadable type"):
            reader._ensure_open()
    assert not reader.is_open()
    handle.close.assert_called_once_with()


@pytest.mark.requires_pybigwig
def test_tiny_file_values_gaps_and_exact_statistics(
    tiny_bw_path: pathlib.Path,
) -> None:
    """Check every nonempty subregion against an independent dense oracle."""
    signal = [
        1.0,
        1.0,
        math.nan,
        math.nan,
        0.0,
        0.0,
        math.nan,
        math.nan,
        3.0,
        3.0,
        math.nan,
        math.nan,
        -2.0,
        -2.0,
        math.nan,
        math.nan,
    ]
    reader = BigWigFile(tiny_bw_path)
    with reader:
        for start in range(16):
            for end in range(start + 1, 17):
                expected = signal[start:end]
                np.testing.assert_allclose(
                    reader.values("chr1", start, end), expected, equal_nan=True
                )
                covered = [value for value in expected if not math.isnan(value)]
                for kind, statistic in [
                    ("mean", sum(covered) / len(covered) if covered else None),
                    ("min", min(covered) if covered else None),
                    ("max", max(covered) if covered else None),
                    ("sum", sum(covered) if covered else None),
                    (
                        "coverage",
                        len(covered) / len(expected) if covered else None,
                    ),
                ]:
                    assert reader.stats(
                        "chr1", start, end, type=kind, exact=True
                    ) == pytest.approx([statistic])


@pytest.mark.requires_pybigwig
def test_binned_statistics_preserve_zero_and_missing(
    tiny_bw_path: pathlib.Path,
) -> None:
    reader = BigWigFile(tiny_bw_path)
    with reader:
        assert reader.stats("chr1", n_bins=8, exact=True) == [
            1.0,
            None,
            0.0,
            None,
            3.0,
            None,
            -2.0,
            None,
        ]
        assert reader.stats("chr1", type="std", exact=True) == pytest.approx(
            [
                math.sqrt(26 / 7),
            ]
        )
        header = reader.header()
        assert header["nBasesCovered"] == 8
        assert header["sumData"] == 4
        assert header["sumSquared"] == 28


@pytest.mark.requires_pybigwig
def test_exact_statistics_bypass_zoom_summaries(tmp_path: pathlib.Path) -> None:
    """Use enough entries for zooms; tiny files cannot exercise exact=True."""
    pybigwig = pytest.importorskip("pyBigWig")
    path = tmp_path / "zoom.bw"
    signal = [float((i * 7) % 13) for i in range(4096)]
    with pybigwig.open(str(path), "w") as native:
        native.addHeader([("chr1", len(signal))])
        native.addEntries(
            ["chr1"] * len(signal),
            list(range(len(signal))),
            ends=list(range(1, len(signal) + 1)),
            values=signal,
        )
    reader = BigWigFile(path)
    with reader:
        assert reader.header()["nLevels"] > 0
        expected = [
            sum(signal[start : start + 200]) / 200
            for start in range(3, 1003, 200)
        ]
        assert reader.stats(
            "chr1", 3, 1003, n_bins=5, exact=True
        ) == pytest.approx(expected)
        approximate = reader.stats("chr1", 3, 1003, n_bins=5)
        assert approximate != pytest.approx(expected)
        assert (
            reader.stats("chr1", 3, 1003, n_bins=5, exact=False) == approximate
        )


@pytest.mark.requires_pybigwig
def test_intervals_keep_original_endpoints(tiny_bw_path: pathlib.Path) -> None:
    reader = BigWigFile(tiny_bw_path)
    with reader:
        assert reader.intervals("chr1", 1, 9) == (
            (0, 2, 1.0),
            (4, 6, 0.0),
            (8, 10, 3.0),
        )


@pytest.mark.requires_pybigwig
def test_large_coordinates_are_unsigned_without_large_arrays(
    tmp_path: pathlib.Path,
) -> None:
    """A sparse 32-bit chromosome needs only a few hundred bytes on disk."""
    pybigwig = pytest.importorskip("pyBigWig")
    path = tmp_path / "large-coordinates.bw"
    intervals = (
        (2**31 - 2, 2**31 + 2, 1.0),
        (2**31 + 10, 2**31 + 12, 2.0),
        (2**32 - 3, 2**32 - 1, 3.0),
    )
    with pybigwig.open(str(path), "w") as native:
        # Disable zoom generation; only three sparse intervals are written.
        native.addHeader([("chr1", 2**32 - 1)], maxZooms=0)
        native.addEntries(
            ["chr1"] * 3,
            [start for start, _, _ in intervals],
            ends=[end for _, end, _ in intervals],
            values=[value for _, _, value in intervals],
        )
    reader = BigWigFile(path)
    with reader:
        assert reader.chroms("chr1") == 2**32 - 1
        assert reader.intervals("chr1") == intervals
        # The first interval ends beyond the query and crosses the sign bit.
        assert reader.intervals("chr1", 2**31 - 2, 2**31) == intervals[:1]
        for start, end, value in intervals:
            assert reader.intervals("chr1", start, end) == (
                (start, end, value),
            )
            assert reader.values("chr1", start, end) == [value] * (end - start)
            assert reader.stats("chr1", start, end, exact=True) == [value]


@pytest.mark.requires_pybigwig
def test_missing_regions_and_empty_file(
    tmp_path: pathlib.Path, tiny_bw_path: pathlib.Path
) -> None:
    reader = BigWigFile(tiny_bw_path)
    with reader:
        assert reader.chroms("missing") is None
        assert reader.intervals("chr1", 2, 4) is None
        assert reader.intervals("empty") is None
        assert reader.stats("empty", n_bins=2, exact=True) == [None, None]
        assert all(math.isnan(value) for value in reader.values("empty", 0, 8))
    pybigwig = pytest.importorskip("pyBigWig")
    path = tmp_path / "empty.bw"
    with pybigwig.open(str(path), "w") as native:
        native.addHeader([("chr1", 8)])
    reader = BigWigFile(path)
    with reader:
        assert reader.header()["nBasesCovered"] == 0
        assert reader.values("chr1", 0, 8) == []
        assert reader.intervals("chr1") is None
        assert reader.stats("chr1", n_bins=2) == [None, None]


@pytest.mark.requires_pybigwig
def test_numpy_values_match_list(tiny_bw_path: pathlib.Path) -> None:
    pybigwig = pytest.importorskip("pyBigWig")
    if not pybigwig.numpy:
        pytest.skip("pyBigWig was compiled without NumPy support")
    reader = BigWigFile(tiny_bw_path)
    with reader:
        array = reader.values("chr1", 0, 16, numpy=True)
        assert isinstance(array, np.ndarray)
        assert array.dtype == np.float32
        np.testing.assert_allclose(
            array, reader.values("chr1", 0, 16), equal_nan=True
        )


@pytest.mark.requires_pybigwig
@pytest.mark.parametrize("method", ["values", "stats", "intervals"])
@pytest.mark.parametrize("bounds", [(2, 2), (4, 2), (0, 17)])
def test_invalid_native_bounds_propagate(
    tiny_bw_path: pathlib.Path,
    method: str,
    bounds: tuple[int, int],
) -> None:
    reader = BigWigFile(tiny_bw_path)
    with reader:
        with pytest.raises(RuntimeError, match="Invalid interval bounds"):
            getattr(reader, method)("chr1", *bounds)
        with pytest.raises(RuntimeError, match="Invalid interval bounds"):
            getattr(reader, method)("missing", 0, 2)


@pytest.mark.parametrize("method", ["chroms", "values", "stats", "intervals"])
@pytest.mark.parametrize(
    "chrom, error",
    [
        (123, TypeError),
        (b"chr1", TypeError),
        (["chr1"], TypeError),
        ("chr1\0suffix", ValueError),
        ("chr\ud800", UnicodeEncodeError),
    ],
)
def test_invalid_chromosome_rejected_before_open(
    method: str, chrom: Any, error: type[Exception]
) -> None:
    """Do not let native string-parser errors escape as SystemError."""
    module, handle = _mock_pybigwig()
    reader = BigWigFile(_BW_PATH)
    args = () if method == "chroms" else (0, 1)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        with pytest.raises(error):
            getattr(reader, method)(chrom, *args)
    module.open.assert_not_called()
    getattr(handle, method).assert_not_called()
    assert not reader.is_open()


@pytest.mark.parametrize("argument", ["path", "mode"])
@pytest.mark.parametrize(
    "suffix, error", [("\0suffix", ValueError), ("\ud800", UnicodeEncodeError)]
)
def test_invalid_open_strings_leave_existing_file_untouched(
    tmp_path: pathlib.Path, argument: str, suffix: str, error: type[Exception]
) -> None:
    """Validate strings before handing a potentially destructive open to C."""
    path = tmp_path / "existing.bw"
    path.write_bytes(b"original contents")
    args = {"path": str(path), "mode": "w"}
    args[argument] += suffix
    reader = BigWigFile(**args)
    module, _ = _mock_pybigwig()
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        with pytest.raises(error):
            with reader:
                pytest.fail("Invalid strings must never reach an open handle")
    module.open.assert_not_called()
    assert not reader.is_open()
    assert path.read_bytes() == b"original contents"


@pytest.mark.requires_pybigwig
def test_unicode_names_and_paths(tmp_path: pathlib.Path) -> None:
    """Valid UTF-8 and str subclasses remain supported at the C boundary."""
    pybigwig = pytest.importorskip("pyBigWig")
    path = tmp_path / "covérage.bw"
    chrom = np.str_("chrα")
    writer = BigWigFile(path, mode="w")
    with writer as native:
        native.addHeader([(str(chrom), 4)], maxZooms=0)
        native.addEntries([str(chrom)], [0], ends=[4], values=[1.25])
    reader = BigWigFile(path)
    with reader:
        assert reader.chroms() == {str(chrom): 4}
        assert reader.chroms(chrom) == 4
        assert reader.values(chrom, 0, 4) == [1.25] * 4
        assert reader.stats(chrom, type="sum", exact=True) == [5.0]
        assert reader.intervals(chrom) == ((0, 4, 1.25),)


@pytest.mark.parametrize(
    "kind, error",
    [
        (123, TypeError),
        (b"mean", TypeError),
        (["mean"], TypeError),
        ("mean\0suffix", ValueError),
        ("\ud800", UnicodeEncodeError),
    ],
)
def test_invalid_statistic_string_rejected_before_open(
    kind: Any, error: type[Exception]
) -> None:
    module, handle = _mock_pybigwig()
    reader = BigWigFile(_BW_PATH)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        with pytest.raises(error):
            reader.stats("chr1", type=kind)
    module.open.assert_not_called()
    handle.stats.assert_not_called()


@pytest.mark.requires_pybigwig
def test_unknown_statistic_propagates(tiny_bw_path: pathlib.Path) -> None:
    reader = BigWigFile(tiny_bw_path)
    with reader:
        with pytest.raises(RuntimeError, match="Invalid type"):
            reader.stats("chr1", type="unknown")


@pytest.mark.parametrize("method", ["chroms", "values", "stats", "intervals"])
@pytest.mark.parametrize("suffix", ["\0suffix", "\ud800"])
def test_string_subclasses_cannot_bypass_chromosome_validation(
    method: str, suffix: str
) -> None:
    module, handle = _mock_pybigwig()
    reader = BigWigFile(_BW_PATH)
    args = () if method == "chroms" else (0, 1)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        with pytest.raises(ValueError):
            getattr(reader, method)(_StringSubclass("chr1" + suffix), *args)
    module.open.assert_not_called()
    getattr(handle, method).assert_not_called()


@pytest.mark.parametrize("suffix", ["\0suffix", "\ud800"])
def test_string_subclasses_cannot_bypass_statistic_validation(
    suffix: str,
) -> None:
    module, handle = _mock_pybigwig()
    reader = BigWigFile(_BW_PATH)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        with pytest.raises(ValueError):
            reader.stats("chr1", type=_StringSubclass("mean" + suffix))
    module.open.assert_not_called()
    handle.stats.assert_not_called()


@pytest.mark.parametrize("suffix", ["\0suffix", "\ud800"])
def test_string_subclasses_cannot_bypass_mode_validation(
    tmp_path: pathlib.Path, suffix: str
) -> None:
    path = tmp_path / "existing.bw"
    path.write_bytes(b"original contents")
    reader = BigWigFile(path, mode=_StringSubclass("w" + suffix))
    module, _ = _mock_pybigwig()
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        with pytest.raises(ValueError):
            with reader:
                pytest.fail("Invalid modes must never open a native handle")
    module.open.assert_not_called()
    assert path.read_bytes() == b"original contents"


@pytest.mark.parametrize(
    "method", ["chroms", "header", "values", "stats", "intervals"]
)
@pytest.mark.parametrize("already_open", [False, True])
def test_string_subclass_write_mode_cannot_bypass_read_guard(
    tmp_path: pathlib.Path, method: str, already_open: bool
) -> None:
    path = tmp_path / "existing.bw"
    path.write_bytes(b"original contents")
    reader = BigWigFile(path, mode=_StringSubclass("w"))
    module, handle = _mock_pybigwig()
    args = () if method in {"chroms", "header"} else ("chr1", 0, 1)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        try:
            if already_open:
                reader.__enter__()
                reader.mode = "r"
            with pytest.raises(RuntimeError, match="write mode"):
                getattr(reader, method)(*args)
            if not already_open:
                module.open.assert_not_called()
            getattr(handle, method).assert_not_called()
        finally:
            reader.close()
    assert path.read_bytes() == b"original contents"


def test_string_subclasses_preserve_native_contents() -> None:
    """Normalization must not call overridden __str__ or change valid names."""
    module, handle = _mock_pybigwig()
    reader = BigWigFile(_BW_PATH, mode=_StringSubclass("r"))
    chrom = _StringSubclass("chr1")
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        with reader:
            reader.chroms(chrom)
            reader.values(chrom, 0, 1)
            reader.stats(chrom, type=_StringSubclass("sum"))
            reader.intervals(chrom)
    assert type(module.open.call_args.args[1]) is str
    assert module.open.call_args.args[1] == "r"
    for method in (
        handle.chroms,
        handle.values,
        handle.stats,
        handle.intervals,
    ):
        assert type(method.call_args.args[0]) is str
        assert method.call_args.args[0] == "chr1"
    assert type(handle.stats.call_args.kwargs["type"]) is str
    assert handle.stats.call_args.kwargs["type"] == "sum"


@pytest.mark.parametrize("mode", [None, b"w", ["w"], {"w"}])
@pytest.mark.parametrize("method", ["__enter__", "header"])
def test_invalid_mode_types_rejected_before_open(
    mode: Any, method: str
) -> None:
    module, _ = _mock_pybigwig()
    reader = BigWigFile(_BW_PATH, mode=mode)
    with mock.patch.object(_utils, "get_pybigwig_module", return_value=module):
        with pytest.raises(TypeError, match="mode must be a string"):
            getattr(reader, method)()
    module.open.assert_not_called()


@pytest.mark.requires_pybigwig
def test_native_string_subclasses(tiny_bw_path: pathlib.Path) -> None:
    """Preserve valid native strings and reject disguised writes intact."""
    reader = BigWigFile(tiny_bw_path, mode=_StringSubclass("r"))
    with reader:
        chrom = _StringSubclass("chr1")
        assert reader.chroms(chrom) == 16
        assert reader.values(chrom, 0, 2) == [1.0, 1.0]
        assert reader.stats(chrom, type=_StringSubclass("sum"), exact=True) == [
            4.0
        ]
        assert reader.intervals(chrom, 0, 2) == ((0, 2, 1.0),)

    original = tiny_bw_path.read_bytes()
    reader.mode = _StringSubclass("w")
    try:
        with pytest.raises(RuntimeError, match="write mode"):
            # stats() safely rejects native writers, even without a header.
            reader.stats("chr1", 0, 1)
    finally:
        reader.close()
    assert not reader.is_open()
    assert tiny_bw_path.read_bytes() == original


@pytest.mark.requires_pybigwig
@pytest.mark.parametrize("encoding", ["bedGraph", "variableStep", "fixedStep"])
@pytest.mark.parametrize("offset", [0, 2**31 - 2])
def test_queries_across_storage_encodings(
    tmp_path: pathlib.Path, encoding: str, offset: int
) -> None:
    """Exercise all native block decoders with tiny, sparse signed data."""
    pybigwig = pytest.importorskip("pyBigWig")
    path = tmp_path / "encoding.bw"
    starts = [offset + i for i in (2, 6, 10, 14)]
    values = [-2.0, 0.0, 1.25, 3.5]
    stored = tuple((a, a + 2, v) for a, v in zip(starts, values))
    with pybigwig.open(str(path), "w") as native:
        native.addHeader([("chr1", offset + 20)], maxZooms=0)
        if encoding == "bedGraph":
            native.addEntries(
                ["chr1"] * 4,
                starts,
                ends=[a + 2 for a in starts],
                values=values,
            )
        elif encoding == "variableStep":
            native.addEntries("chr1", starts, values=values, span=2)
        else:
            native.addEntries("chr1", starts[0], values=values, span=2, step=4)

    dense = [math.nan] * 20
    for start, end, value in stored:
        dense[start - offset : end - offset] = [value] * (end - start)
    reader = BigWigFile(path)
    with reader:
        assert reader.intervals("chr1") == stored
        np.testing.assert_allclose(
            reader.values("chr1", offset, offset + 20), dense, equal_nan=True
        )
        # Four unevenly covered bins, each five bases wide.
        expected: dict[str, list[float]] = {
            kind: [] for kind in ("mean", "sum", "min", "max", "coverage")
        }
        for left in range(0, 20, 5):
            data = [v for v in dense[left : left + 5] if not math.isnan(v)]
            for kind, value in (
                ("mean", sum(data) / len(data)),
                ("sum", sum(data)),
                ("min", min(data)),
                ("max", max(data)),
                ("coverage", len(data) / 5),
            ):
                expected[kind].append(value)
        for kind, result in expected.items():
            assert reader.stats(
                "chr1", offset, offset + 20, type=kind, n_bins=4, exact=True
            ) == pytest.approx(result)
