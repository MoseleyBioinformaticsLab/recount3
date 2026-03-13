from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest import mock

import pytest

from recount3 import _utils
from recount3._bigwig import BigWigFile

_DATA_DIR = Path(__file__).parent / "data"
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
        h1 = bw._ensure_open()  # pylint: disable=protected-access
        h2 = bw._ensure_open()  # pylint: disable=protected-access
    bw.close()
    assert h1 is h2
    assert mock_module.open.call_count == 1


def test_ensure_open_missing_path_raises_file_not_found(
    tmp_path: Path,
) -> None:
    """FileNotFoundError is raised when the path does not exist on disk."""
    bw = BigWigFile(path=tmp_path / "nonexistent.bw")
    with pytest.raises(FileNotFoundError):
        bw._ensure_open()  # pylint: disable=protected-access


def test_ensure_open_pybigwig_returns_none_raises_runtime_error() -> None:
    """RuntimeError is raised when pyBigWig.open() returns None."""
    mock_module, _ = _mock_pybigwig(open_return=None)
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        with pytest.raises(RuntimeError, match="Failed to open BigWig file"):
            bw._ensure_open()  # pylint: disable=protected-access


def test_ensure_open_import_error_propagates() -> None:
    """ImportError from get_pybigwig_module propagates unchanged."""
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils,
        "get_pybigwig_module",
        side_effect=ImportError("pyBigWig not found"),
    ):
        with pytest.raises(ImportError, match="pyBigWig not found"):
            bw._ensure_open()  # pylint: disable=protected-access


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
        bw._ensure_open()  # pylint: disable=protected-access
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
        bw._ensure_open()  # pylint: disable=protected-access
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
        bw._ensure_open()  # pylint: disable=protected-access
        assert bw.is_open() is True
    bw.close()


def test_is_open_false_after_close() -> None:
    """is_open() returns False after close()."""
    mock_module, _ = _mock_pybigwig()
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        bw._ensure_open()  # pylint: disable=protected-access
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
        bw._ensure_open()  # pylint: disable=protected-access
        bw.__exit__(TypeError, TypeError("e"), None)
    assert not bw.is_open()


def test_chroms_all_returns_full_mapping() -> None:
    """chroms() returns a dict mapping every chromosome to its length."""
    pytest.importorskip("pyBigWig")
    bw = BigWigFile(path=_BW_PATH)
    result = bw.chroms()
    assert isinstance(result, dict)
    assert result.get(_BW_CHROM) == _BW_CHROM_LEN
    bw.close()


def test_chroms_single_returns_integer_length() -> None:
    """chroms(chrom) returns the integer length for the named chromosome."""
    pytest.importorskip("pyBigWig")
    bw = BigWigFile(path=_BW_PATH)
    result = bw.chroms(_BW_CHROM)
    assert result == _BW_CHROM_LEN
    bw.close()


def test_header_returns_dict_with_standard_fields() -> None:
    """header() returns a dict that includes standard BigWig header keys."""
    pytest.importorskip("pyBigWig")
    bw = BigWigFile(path=_BW_PATH)
    h = bw.header()
    assert isinstance(h, dict)
    assert "version" in h
    assert "nBasesCovered" in h
    bw.close()


def test_values_returns_list_of_per_base_values() -> None:
    """values() returns exactly one float per base in the requested interval."""
    pytest.importorskip("pyBigWig")
    bw = BigWigFile(path=_BW_PATH)
    result = bw.values(_BW_CHROM, 0, 10)
    assert isinstance(result, list)
    assert len(result) == 10
    bw.close()


def test_values_numpy_flag_forwarded_verbatim() -> None:
    """The numpy keyword is forwarded exactly as-is to pyBigWig.values."""
    mock_module, mock_handle = _mock_pybigwig()
    mock_handle.values.return_value = [1.0, 2.0]
    bw = BigWigFile(path=_BW_PATH)
    with mock.patch.object(
        _utils, "get_pybigwig_module", return_value=mock_module
    ):
        bw.values(_BW_CHROM, 0, 2, numpy=True)
    mock_handle.values.assert_called_once_with(_BW_CHROM, 0, 2, numpy=True)


def test_stats_whole_chrom_no_coords() -> None:
    """stats(chrom) with no start/end summarises the whole chromosome."""
    pytest.importorskip("pyBigWig")
    bw = BigWigFile(path=_BW_PATH)
    result = bw.stats(_BW_CHROM)
    assert isinstance(result, list)
    assert len(result) == 1
    bw.close()


def test_stats_with_start_and_end_coords() -> None:
    """stats(chrom, start, end) returns stats for the given interval."""
    pytest.importorskip("pyBigWig")
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
        bw.stats(_BW_CHROM, 0, 2_000, n_bins=2, exact=False)
    kw = mock_handle.stats.call_args.kwargs
    assert "nBins" in kw
    assert "exact" in kw


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


def test_intervals_with_coords_returns_region_tuples() -> None:
    """intervals(chrom, start, end) returns (start, end, value) tuples."""
    pytest.importorskip("pyBigWig")
    bw = BigWigFile(path=_BW_PATH)
    result = bw.intervals(_BW_CHROM, 0, 1_000)
    assert result is not None
    assert isinstance(result, tuple)
    assert result[0] == (0, 1000, 0.0)
    bw.close()
