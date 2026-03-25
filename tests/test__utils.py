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

from __future__ import annotations

import errno
import io
import os
import ssl
import types
import urllib.error
import zipfile
from pathlib import Path
from unittest import mock

import pandas as pd
import pytest

import recount3._utils as _utils  # pylint: disable=protected-access
from recount3 import errors

_DATA_DIR = Path(__file__).parent / "data"
_MIRROR = _DATA_DIR / "recount3_mirror" / "recount3"

_GTF_GZ = (
    _MIRROR
    / "human"
    / "annotations"
    / "gene_sums"
    / "human.gene_sums.G026.gtf.gz"
)

_EXON_GZ = (
    _MIRROR
    / "human"
    / "data_sources"
    / "sra"
    / "exon_sums"
    / "65"
    / "SRP014565"
    / "sra.exon_sums.SRP014565.G026.gz"
)

_SAMPLE_URL = (
    "http://duffel.rail.bio/recount3"
    "/human/data_sources/sra/gene_sums/65/SRP014565"
    "/sra.gene_sums.SRP014565.G026.gz"
)

_BASE_JXN_URL = (
    "http://duffel.rail.bio/recount3"
    "/human/data_sources/sra/junctions/65/SRP014565"
    "/sra.junctions.SRP014565.ALL"
)


def test_weakreflock_context_manager_normal() -> None:
    """_WeakRefLock enters and exits without error in normal flow."""
    lock = _utils._WeakRefLock()  # pylint: disable=protected-access
    entered = []
    with lock:
        entered.append(True)
    assert entered == [True]


def test_weakreflock_context_manager_on_exception() -> None:
    """_WeakRefLock releases its inner lock even when body raises."""
    lock = _utils._WeakRefLock()  # pylint: disable=protected-access
    with pytest.raises(RuntimeError, match="test error"):
        with lock:
            raise RuntimeError("test error")
    acquired = lock._lock.acquire(  # pylint: disable=protected-access
        blocking=False
    )
    assert acquired
    lock._lock.release()  # pylint: disable=protected-access


def test_zip_lock_for_path_creates_lock(tmp_path: Path) -> None:
    """A new ZIP path gets its own _WeakRefLock."""
    zip_path = tmp_path / "archive.zip"
    lock = _utils._zip_lock_for_path(  # pylint: disable=protected-access
        zip_path
    )
    assert isinstance(
        lock,
        _utils._WeakRefLock,  # pylint: disable=protected-access
    )


def test_zip_lock_for_path_reuses_lock(tmp_path: Path) -> None:
    """The same ZIP path always returns the identical lock object."""
    zip_path = tmp_path / "archive.zip"
    lock_a = _utils._zip_lock_for_path(  # pylint: disable=protected-access
        zip_path
    )
    lock_b = _utils._zip_lock_for_path(  # pylint: disable=protected-access
        zip_path
    )
    assert lock_a is lock_b


def test_zip_lock_for_path_different_paths(tmp_path: Path) -> None:
    """Different ZIP paths produce independent lock objects."""
    lock_a = _utils._zip_lock_for_path(  # pylint: disable=protected-access
        tmp_path / "a.zip"
    )
    lock_b = _utils._zip_lock_for_path(  # pylint: disable=protected-access
        tmp_path / "b.zip"
    )
    assert lock_a is not lock_b


def test_sha256_returns_64_char_hex_digest() -> None:
    """_sha256 returns a 64-character lowercase hex string."""
    digest = _utils._sha256("hello")  # pylint: disable=protected-access
    assert len(digest) == 64
    assert all(c in "0123456789abcdef" for c in digest)


def test_sha256_deterministic() -> None:
    """Same input always produces the same digest."""
    a = _utils._sha256("recount3")  # pylint: disable=protected-access
    b = _utils._sha256("recount3")  # pylint: disable=protected-access
    assert a == b


def test_sha256_differs_for_different_inputs() -> None:
    """Different inputs yield different digests."""
    assert (
        _utils._sha256("foo")  # pylint: disable=protected-access
        != _utils._sha256("bar")  # pylint: disable=protected-access
    )


def test_sha256_empty_string() -> None:
    """_sha256 handles an empty string without error."""
    digest = _utils._sha256("")  # pylint: disable=protected-access
    assert len(digest) == 64


def test_cache_key_for_url_format() -> None:
    """Key is '<16-hex>__<basename>' and basename matches the URL."""
    key = _utils._cache_key_for_url(  # pylint: disable=protected-access
        _SAMPLE_URL
    )
    prefix, basename = key.split("__", 1)
    assert len(prefix) == 16
    assert all(c in "0123456789abcdef" for c in prefix)
    assert basename == "sra.gene_sums.SRP014565.G026.gz"


def test_cache_key_for_url_deterministic() -> None:
    """Same URL always produces the same key."""
    assert _utils._cache_key_for_url(  # pylint: disable=protected-access
        _SAMPLE_URL
    ) == _utils._cache_key_for_url(  # pylint: disable=protected-access
        _SAMPLE_URL
    )


def test_cache_key_for_url_different_urls_same_basename() -> None:
    """Two URLs sharing a basename but different paths produce different keys."""
    url_a = "https://host-a.example.com/path/A/file.gz"
    url_b = "https://host-b.example.com/path/B/file.gz"
    key_a = _utils._cache_key_for_url(  # pylint: disable=protected-access
        url_a
    )
    key_b = _utils._cache_key_for_url(  # pylint: disable=protected-access
        url_b
    )
    assert key_a != key_b
    assert key_a.endswith("__file.gz")
    assert key_b.endswith("__file.gz")


def test_cache_path_returns_path_under_root(tmp_path: Path) -> None:
    """Returned Path lives directly under the given cache root."""
    p = _utils._cache_path(  # pylint: disable=protected-access
        _SAMPLE_URL, tmp_path
    )
    assert p.parent == tmp_path


def test_cache_path_creates_cache_root(tmp_path: Path) -> None:
    """cache_root is created when it does not yet exist."""
    root = tmp_path / "new" / "cache"
    assert not root.exists()
    _utils._cache_path(  # pylint: disable=protected-access
        _SAMPLE_URL, root
    )
    assert root.is_dir()


def test_cache_path_accepts_string_root(tmp_path: Path) -> None:
    """cache_root may be passed as a str."""
    p = _utils._cache_path(  # pylint: disable=protected-access
        _SAMPLE_URL, str(tmp_path)
    )
    assert isinstance(p, Path)
    assert p.parent == tmp_path


def test_ensure_dir_creates_nested_directories(tmp_path: Path) -> None:
    """_ensure_dir creates a deeply nested directory tree."""
    target = tmp_path / "a" / "b" / "c"
    _utils._ensure_dir(target)  # pylint: disable=protected-access
    assert target.is_dir()


def test_ensure_dir_idempotent_on_existing_dir(tmp_path: Path) -> None:
    """_ensure_dir does not raise when the directory already exists."""
    _utils._ensure_dir(tmp_path)  # pylint: disable=protected-access
    _utils._ensure_dir(tmp_path)  # Second call – must not raise.


def test_ensure_dir_raises_when_file_occupies_path(
    tmp_path: Path,
) -> None:
    """_ensure_dir raises NotADirectoryError when a file is at path."""
    file_path = tmp_path / "regular_file"
    file_path.write_bytes(b"content")
    with pytest.raises(NotADirectoryError):
        _utils._ensure_dir(  # pylint: disable=protected-access
            file_path
        )


def test_atomic_replace_replaces_content(tmp_path: Path) -> None:
    """_atomic_replace atomically moves src_tmp over final_path."""
    src = tmp_path / "tmp_file"
    dst = tmp_path / "final_file"
    src.write_bytes(b"new content")
    dst.write_bytes(b"old content")
    _utils._atomic_replace(src, dst)  # pylint: disable=protected-access
    assert dst.read_bytes() == b"new content"
    assert not src.exists()


def test_atomic_replace_creates_parent_directory(tmp_path: Path) -> None:
    """_atomic_replace creates the parent directory if it is absent."""
    src = tmp_path / "src"
    src.write_bytes(b"data")
    dst = tmp_path / "subdir" / "final"
    _utils._atomic_replace(src, dst)  # pylint: disable=protected-access
    assert dst.read_bytes() == b"data"


def test_hardlink_or_copy_hardlink_success(tmp_path: Path) -> None:
    """_hardlink_or_copy creates a hardlink on the same filesystem."""
    dst = tmp_path / "dst.gtf.gz"
    _utils._hardlink_or_copy(  # pylint: disable=protected-access
        _GTF_GZ, dst
    )
    assert dst.exists()
    assert dst.read_bytes() == _GTF_GZ.read_bytes()


def test_hardlink_or_copy_cross_device_fallback(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """EXDEV from os.link triggers shutil.copy2 fallback."""
    monkeypatch.setattr(
        _utils.os,
        "link",
        mock.Mock(side_effect=OSError(errno.EXDEV, "cross-device")),
    )
    dst = tmp_path / "dst.gtf.gz"
    _utils._hardlink_or_copy(  # pylint: disable=protected-access
        _GTF_GZ, dst
    )
    assert dst.read_bytes() == _GTF_GZ.read_bytes()


@pytest.mark.parametrize(
    "err_no",
    [errno.EPERM, errno.EACCES, errno.EMLINK],
)
def test_hardlink_or_copy_permission_errors_fall_back_to_copy(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    err_no: int,
) -> None:
    """EPERM / EACCES / EMLINK all trigger the copy fallback."""
    monkeypatch.setattr(
        _utils.os,
        "link",
        mock.Mock(
            side_effect=OSError(err_no, os.strerror(err_no))
        ),
    )
    dst = tmp_path / "dst.gtf.gz"
    _utils._hardlink_or_copy(  # pylint: disable=protected-access
        _GTF_GZ, dst
    )
    assert dst.read_bytes() == _GTF_GZ.read_bytes()


def test_hardlink_or_copy_unexpected_oserror_propagates(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An unhandled OSError from os.link propagates to the caller."""
    monkeypatch.setattr(
        _utils.os,
        "link",
        mock.Mock(side_effect=OSError(errno.EIO, "I/O error")),
    )
    with pytest.raises(OSError, match="I/O error"):
        _utils._hardlink_or_copy(  # pylint: disable=protected-access
            _GTF_GZ, tmp_path / "dst"
        )


def test_hardlink_or_copy_tmp_cleaned_up_on_replace_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Temp file is removed when _atomic_replace fails."""
    monkeypatch.setattr(
        _utils.os,
        "link",
        mock.Mock(side_effect=OSError(errno.EXDEV, "cross-device")),
    )
    monkeypatch.setattr(
        _utils.os,
        "replace",
        mock.Mock(side_effect=OSError("replace failed")),
    )
    with pytest.raises(OSError, match="replace failed"):
        _utils._hardlink_or_copy(  # pylint: disable=protected-access
            _GTF_GZ, tmp_path / "dst.gtf.gz"
        )
    assert list(tmp_path.glob(".*.tmp")) == []


def test_hardlink_or_copy_unlink_oserror_in_finally_swallowed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """OSError during tmp cleanup in finally is silently swallowed."""
    monkeypatch.setattr(
        _utils.os,
        "link",
        mock.Mock(side_effect=OSError(errno.EXDEV, "cross-device")),
    )
    monkeypatch.setattr(
        _utils.os,
        "replace",
        mock.Mock(side_effect=OSError("replace failed")),
    )

    original_unlink = Path.unlink

    def _raise_for_tmp(
        self: Path, missing_ok: bool = False
    ) -> None:
        if ".tmp" in self.name:
            raise OSError("cannot unlink tmp")
        original_unlink(self, missing_ok=missing_ok)

    monkeypatch.setattr(Path, "unlink", _raise_for_tmp)

    with pytest.raises(OSError, match="replace failed"):
        _utils._hardlink_or_copy(  # pylint: disable=protected-access
            _GTF_GZ, tmp_path / "dst"
        )


def test_ssl_insecure_context_disables_verification() -> None:
    """_ssl_insecure_context returns a context with verification off."""
    ctx = _utils._ssl_insecure_context()  # pylint: disable=protected-access
    assert isinstance(ctx, ssl.SSLContext)
    assert not ctx.check_hostname
    assert ctx.verify_mode == ssl.CERT_NONE


def test_http_open_raises_on_empty_url() -> None:
    """http_open raises ValueError for an empty URL."""
    with pytest.raises(ValueError, match="Empty URL"):
        _utils.http_open(
            "",
            timeout=10,
            headers=None,
            insecure_ssl=False,
            user_agent="test-agent",
        )


def test_http_open_sends_user_agent(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """http_open includes the User-Agent header in the request."""
    fake_resp = io.BytesIO(b"body")
    mock_opener = mock.MagicMock()
    mock_opener.open.return_value = fake_resp

    monkeypatch.setattr(
        _utils.urllib.request,
        "build_opener",
        mock.Mock(return_value=mock_opener),
    )

    result = _utils.http_open(
        "http://example.com/file",
        timeout=5,
        headers=None,
        insecure_ssl=False,
        user_agent="recount3-test/1.0",
    )

    assert result is fake_resp
    req = mock_opener.open.call_args[0][0]
    assert req.get_header("User-agent") == "recount3-test/1.0"


def test_http_open_merges_extra_headers_skips_none(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Extra headers are merged; None-valued entries are dropped."""
    fake_resp = io.BytesIO(b"")
    mock_opener = mock.MagicMock()
    mock_opener.open.return_value = fake_resp
    monkeypatch.setattr(
        _utils.urllib.request,
        "build_opener",
        mock.Mock(return_value=mock_opener),
    )

    _utils.http_open(
        "http://example.com/x",
        timeout=5,
        headers={
            "X-Custom": "yes",
            "X-Skip": None,  # type: ignore[dict-item]
        },
        insecure_ssl=False,
        user_agent="ua",
    )

    req = mock_opener.open.call_args[0][0]
    assert req.get_header("X-custom") == "yes"
    assert req.get_header("X-skip") is None


def test_http_open_adds_https_handler_for_insecure_ssl(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """insecure_ssl=True calls add_handler with an HTTPSHandler."""
    fake_resp = io.BytesIO(b"")
    mock_opener = mock.MagicMock()
    mock_opener.open.return_value = fake_resp
    monkeypatch.setattr(
        _utils.urllib.request,
        "build_opener",
        mock.Mock(return_value=mock_opener),
    )

    _utils.http_open(
        "https://example.com/x",
        timeout=5,
        headers=None,
        insecure_ssl=True,
        user_agent="ua",
    )

    # The handler is added via opener.add_handler(), not build_opener().
    mock_opener.add_handler.assert_called_once()
    handler = mock_opener.add_handler.call_args[0][0]
    assert isinstance(
        handler, _utils.urllib.request.HTTPSHandler
    )


def test_with_retries_success_first_attempt() -> None:
    """with_retries returns immediately when func succeeds."""
    result = _utils.with_retries(lambda: 42, attempts=3)
    assert result == 42


def test_with_retries_retries_after_url_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """with_retries retries on URLError and returns on success."""
    monkeypatch.setattr(_utils.time, "sleep", mock.Mock())
    call_count = 0

    def _flaky() -> str:
        nonlocal call_count
        call_count += 1
        if call_count == 1:
            raise urllib.error.URLError("temporary failure")
        return "ok"

    result = _utils.with_retries(_flaky, attempts=2)
    assert result == "ok"
    assert call_count == 2
    _utils.time.sleep.assert_called_once()  # type: ignore[attr-defined]


def test_with_retries_all_attempts_exhausted_reraises(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """with_retries re-raises the last exception after all attempts."""
    monkeypatch.setattr(_utils.time, "sleep", mock.Mock())
    with pytest.raises(urllib.error.URLError, match="always fails"):
        _utils.with_retries(
            mock.Mock(
                side_effect=urllib.error.URLError("always fails")
            ),
            attempts=2,
        )


def test_with_retries_zero_attempts_treated_as_one(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """attempts<=0 is coerced to 1 via max(1, attempts)."""
    monkeypatch.setattr(_utils.time, "sleep", mock.Mock())
    called: list[int] = []

    def _func() -> None:
        called.append(1)
        raise urllib.error.URLError("fail")

    with pytest.raises(urllib.error.URLError):
        _utils.with_retries(_func, attempts=0)

    assert len(called) == 1


def test_with_retries_non_retryable_propagates_immediately() -> None:
    """A ValueError (not retryable) propagates on the first attempt."""
    calls: list[int] = []

    def _func() -> None:
        calls.append(1)
        raise ValueError("not retryable")

    with pytest.raises(ValueError, match="not retryable"):
        _utils.with_retries(_func, attempts=5)

    assert len(calls) == 1


def test_with_retries_sleeps_between_but_not_after_last(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Sleep fires between attempts but NOT after the final failure."""
    sleep_mock = mock.Mock()
    monkeypatch.setattr(_utils.time, "sleep", sleep_mock)

    with pytest.raises(ConnectionError):
        _utils.with_retries(
            mock.Mock(side_effect=ConnectionError("down")),
            attempts=3,
            base_sleep=0.1,
        )

    # 3 attempts, 2 sleeps.
    assert sleep_mock.call_count == 2


def test_stream_copy_empty_source_returns_zero() -> None:
    """_stream_copy returns 0 bytes for an empty source."""
    src = io.BytesIO(b"")
    dst = io.BytesIO()
    n = _utils._stream_copy(  # pylint: disable=protected-access
        src, dst, chunk_size=1024
    )
    assert n == 0
    assert dst.getvalue() == b""


def test_stream_copy_single_chunk(tmp_path: Path) -> None:
    """_stream_copy copies data smaller than one chunk correctly."""
    data = _GTF_GZ.read_bytes()
    src = io.BytesIO(data)
    dst = io.BytesIO()
    n = _utils._stream_copy(  # pylint: disable=protected-access
        src, dst, chunk_size=len(data) + 1
    )
    assert n == len(data)
    assert dst.getvalue() == data


def test_stream_copy_multiple_chunks(tmp_path: Path) -> None:
    """_stream_copy reassembles data spread over many small chunks."""
    data = _EXON_GZ.read_bytes()
    src = io.BytesIO(data)
    dst = io.BytesIO()
    n = _utils._stream_copy(  # pylint: disable=protected-access
        src, dst, chunk_size=64
    )
    assert n == len(data)
    assert dst.getvalue() == data


def test_stream_copy_returns_total_byte_count() -> None:
    """Return value exactly equals the number of bytes copied."""
    payload = b"x" * 500
    assert (
        _utils._stream_copy(  # pylint: disable=protected-access
            io.BytesIO(payload), io.BytesIO(), chunk_size=128
        )
        == 500
    )


def test_download_to_file_success(tmp_path: Path) -> None:
    """download_to_file writes URL content to the output path."""
    content = _GTF_GZ.read_bytes()
    out = tmp_path / "output.gtf.gz"

    with mock.patch(
        "recount3._utils.http_open",
        return_value=io.BytesIO(content),
    ):
        _utils.download_to_file(
            "http://example.com/file.gtf.gz",
            out,
            chunk_size=512,
            timeout=30,
            insecure_ssl=False,
            user_agent="test",
            attempts=1,
        )

    assert out.read_bytes() == content


def test_download_to_file_raises_download_error_on_network_failure(
    tmp_path: Path,
) -> None:
    """A network error is wrapped and re-raised as DownloadError."""
    with mock.patch(
        "recount3._utils.http_open",
        side_effect=urllib.error.URLError("connection refused"),
    ):
        with pytest.raises(errors.DownloadError):
            _utils.download_to_file(
                "http://example.com/file.gtf.gz",
                tmp_path / "out",
                chunk_size=512,
                timeout=30,
                insecure_ssl=False,
                user_agent="test",
                attempts=1,
            )


def test_download_to_file_cleans_up_tmp_on_replace_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Temp file is deleted when _atomic_replace fails mid-download."""
    content = _GTF_GZ.read_bytes()
    out = tmp_path / "output.gtf.gz"

    monkeypatch.setattr(
        _utils.os,
        "replace",
        mock.Mock(side_effect=OSError("disk full")),
    )

    with mock.patch(
        "recount3._utils.http_open",
        return_value=io.BytesIO(content),
    ):
        with pytest.raises(errors.DownloadError):
            _utils.download_to_file(
                "http://example.com/file.gtf.gz",
                out,
                chunk_size=512,
                timeout=30,
                insecure_ssl=False,
                user_agent="test",
                attempts=1,
            )

    assert list(tmp_path.glob(".*.downloading")) == []


def test_download_to_file_tmp_not_created_when_http_open_fails(
    tmp_path: Path,
) -> None:
    """If http_open fails, no temp file is created."""
    out = tmp_path / "output.gtf.gz"

    with mock.patch(
        "recount3._utils.http_open",
        side_effect=urllib.error.URLError("refused"),
    ):
        with pytest.raises(errors.DownloadError):
            _utils.download_to_file(
                "http://example.com/file.gtf.gz",
                out,
                chunk_size=512,
                timeout=30,
                insecure_ssl=False,
                user_agent="test",
                attempts=1,
            )

    assert not out.exists()


def test_download_to_file_tmp_unlink_error_swallowed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """OSError during tmp unlink in finally is silently ignored."""
    content = b"data"
    out = tmp_path / "out"

    monkeypatch.setattr(
        _utils.os,
        "replace",
        mock.Mock(side_effect=OSError("replace fail")),
    )

    original_unlink = Path.unlink

    def _raise_for_downloading(
        self: Path, missing_ok: bool = False
    ) -> None:
        if ".downloading" in self.name:
            raise OSError("unlink refused")
        original_unlink(self, missing_ok=missing_ok)

    monkeypatch.setattr(Path, "unlink", _raise_for_downloading)

    with mock.patch(
        "recount3._utils.http_open",
        return_value=io.BytesIO(content),
    ):
        with pytest.raises(errors.DownloadError):
            _utils.download_to_file(
                "http://example.com/f",
                out,
                chunk_size=512,
                timeout=30,
                insecure_ssl=False,
                user_agent="ua",
                attempts=1,
            )


def test_write_or_replace_in_zip_creates_new_zip(
    tmp_path: Path,
) -> None:
    """Writing to a nonexistent ZIP creates the archive."""
    zip_path = tmp_path / "archive.zip"
    _utils._write_or_replace_in_zip(  # pylint: disable=protected-access
        zip_path, _GTF_GZ, "data/gene.gtf.gz", overwrite=False
    )
    assert zip_path.exists()
    with zipfile.ZipFile(zip_path, "r") as zf:
        assert "data/gene.gtf.gz" in zf.namelist()


def test_write_or_replace_in_zip_appends_new_member(
    tmp_path: Path,
) -> None:
    """A new member is appended to an existing ZIP."""
    zip_path = tmp_path / "archive.zip"
    _utils._write_or_replace_in_zip(  # pylint: disable=protected-access
        zip_path, _GTF_GZ, "gene.gtf.gz", overwrite=False
    )
    _utils._write_or_replace_in_zip(  # pylint: disable=protected-access
        zip_path, _EXON_GZ, "exon.gz", overwrite=False
    )
    with zipfile.ZipFile(zip_path, "r") as zf:
        names = zf.namelist()
    assert "gene.gtf.gz" in names
    assert "exon.gz" in names


def test_write_or_replace_in_zip_noop_when_exists_no_overwrite(
    tmp_path: Path,
) -> None:
    """overwrite=False is a no-op when the member already exists."""
    zip_path = tmp_path / "archive.zip"
    original = _GTF_GZ.read_bytes()
    _utils._write_or_replace_in_zip(  # pylint: disable=protected-access
        zip_path, _GTF_GZ, "member.gz", overwrite=False
    )
    _utils._write_or_replace_in_zip(  # pylint: disable=protected-access
        zip_path, _EXON_GZ, "member.gz", overwrite=False
    )
    with zipfile.ZipFile(zip_path, "r") as zf:
        content = zf.read("member.gz")
    assert content == original


def test_write_or_replace_in_zip_overwrites_existing_member(
    tmp_path: Path,
) -> None:
    """overwrite=True replaces the existing member without duplicates."""
    zip_path = tmp_path / "archive.zip"
    _utils._write_or_replace_in_zip(  # pylint: disable=protected-access
        zip_path, _GTF_GZ, "member.gz", overwrite=False
    )
    _utils._write_or_replace_in_zip(  # pylint: disable=protected-access
        zip_path, _EXON_GZ, "member.gz", overwrite=True
    )
    with zipfile.ZipFile(zip_path, "r") as zf:
        names = zf.namelist()
        content = zf.read("member.gz")
    assert names.count("member.gz") == 1
    assert content == _EXON_GZ.read_bytes()


def test_write_or_replace_in_zip_invalid_zip_raises_download_error(
    tmp_path: Path,
) -> None:
    """A non-ZIP file at zip_path raises DownloadError."""
    zip_path = tmp_path / "bad.zip"
    zip_path.write_bytes(b"not a zip file at all")
    with pytest.raises(errors.DownloadError, match="not a valid ZIP"):
        _utils._write_or_replace_in_zip(  # pylint: disable=protected-access
            zip_path, _GTF_GZ, "member.gz", overwrite=False
        )


def test_write_or_replace_in_zip_preserves_other_members_on_overwrite(
    tmp_path: Path,
) -> None:
    """Overwriting one member leaves all other members intact."""
    zip_path = tmp_path / "archive.zip"
    _utils._write_or_replace_in_zip(  # pylint: disable=protected-access
        zip_path, _GTF_GZ, "gene.gz", overwrite=False
    )
    _utils._write_or_replace_in_zip(  # pylint: disable=protected-access
        zip_path, _GTF_GZ, "exon.gz", overwrite=False
    )
    _utils._write_or_replace_in_zip(  # pylint: disable=protected-access
        zip_path, _EXON_GZ, "gene.gz", overwrite=True
    )
    with zipfile.ZipFile(zip_path, "r") as zf:
        assert "exon.gz" in zf.namelist()
        assert zf.read("gene.gz") == _EXON_GZ.read_bytes()


def test_write_or_replace_in_zip_tmpzip_cleaned_up_on_replace_fail(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Temp ZIP is removed even when _atomic_replace fails."""
    zip_path = tmp_path / "archive.zip"
    _utils._write_or_replace_in_zip(  # pylint: disable=protected-access
        zip_path, _GTF_GZ, "member.gz", overwrite=False
    )

    monkeypatch.setattr(
        _utils.os,
        "replace",
        mock.Mock(side_effect=OSError("disk full")),
    )

    with pytest.raises(OSError, match="disk full"):
        _utils._write_or_replace_in_zip(  # pylint: disable=protected-access
            zip_path, _EXON_GZ, "member.gz", overwrite=True
        )

    assert list(tmp_path.glob(".*.tmpzip")) == []


def test_write_or_replace_in_zip_tmpzip_unlink_error_swallowed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """OSError from tmp-ZIP unlink in finally is silently ignored."""
    zip_path = tmp_path / "archive.zip"
    _utils._write_or_replace_in_zip(  # pylint: disable=protected-access
        zip_path, _GTF_GZ, "member.gz", overwrite=False
    )

    monkeypatch.setattr(
        _utils.os,
        "replace",
        mock.Mock(side_effect=OSError("replace fail")),
    )

    original_unlink = Path.unlink

    def _raise_for_tmpzip(
        self: Path, missing_ok: bool = False
    ) -> None:
        if ".tmpzip" in self.name:
            raise OSError("unlink refused")
        original_unlink(self, missing_ok=missing_ok)

    monkeypatch.setattr(Path, "unlink", _raise_for_tmpzip)

    with pytest.raises(OSError, match="replace fail"):
        _utils._write_or_replace_in_zip(  # pylint: disable=protected-access
            zip_path, _EXON_GZ, "member.gz", overwrite=True
        )


def test_download_stream_to_zip_creates_zip(tmp_path: Path) -> None:
    """download_stream_to_zip creates a new ZIP from a download."""
    zip_path = tmp_path / "out.zip"
    content = _GTF_GZ.read_bytes()

    with mock.patch(
        "recount3._utils.http_open",
        return_value=io.BytesIO(content),
    ):
        _utils.download_stream_to_zip(
            "http://example.com/f.gtf.gz",
            zip_path,
            "gene.gtf.gz",
            chunk_size=512,
            overwrite=False,
            timeout=30,
            insecure_ssl=False,
            user_agent="test",
            attempts=1,
        )

    with zipfile.ZipFile(zip_path, "r") as zf:
        assert zf.read("gene.gtf.gz") == content


def test_download_stream_to_zip_early_return_when_member_present(
    tmp_path: Path,
) -> None:
    """overwrite=False skips download when member already in ZIP."""
    zip_path = tmp_path / "out.zip"
    content = _GTF_GZ.read_bytes()

    with mock.patch(
        "recount3._utils.http_open",
        return_value=io.BytesIO(content),
    ):
        _utils.download_stream_to_zip(
            "http://example.com/f.gz",
            zip_path,
            "gene.gz",
            chunk_size=512,
            overwrite=False,
            timeout=30,
            insecure_ssl=False,
            user_agent="test",
            attempts=1,
        )

    with mock.patch("recount3._utils.http_open") as mock_open:
        _utils.download_stream_to_zip(
            "http://example.com/f.gz",
            zip_path,
            "gene.gz",
            chunk_size=512,
            overwrite=False,
            timeout=30,
            insecure_ssl=False,
            user_agent="test",
            attempts=1,
        )
    mock_open.assert_not_called()


def test_download_stream_to_zip_overwrite_replaces_member(
    tmp_path: Path,
) -> None:
    """overwrite=True replaces an existing member."""
    zip_path = tmp_path / "out.zip"

    with mock.patch(
        "recount3._utils.http_open",
        return_value=io.BytesIO(_GTF_GZ.read_bytes()),
    ):
        _utils.download_stream_to_zip(
            "http://example.com/f.gz",
            zip_path,
            "member.gz",
            chunk_size=512,
            overwrite=False,
            timeout=30,
            insecure_ssl=False,
            user_agent="test",
            attempts=1,
        )

    replacement = _EXON_GZ.read_bytes()
    with mock.patch(
        "recount3._utils.http_open",
        return_value=io.BytesIO(replacement),
    ):
        _utils.download_stream_to_zip(
            "http://example.com/f.gz",
            zip_path,
            "member.gz",
            chunk_size=512,
            overwrite=True,
            timeout=30,
            insecure_ssl=False,
            user_agent="test",
            attempts=1,
        )

    with zipfile.ZipFile(zip_path, "r") as zf:
        assert zf.read("member.gz") == replacement


def test_download_stream_to_zip_no_overwrite_absent_member(
    tmp_path: Path,
) -> None:
    """overwrite=False downloads a member that is not yet in the ZIP."""
    zip_path = tmp_path / "out.zip"

    with mock.patch(
        "recount3._utils.http_open",
        return_value=io.BytesIO(_GTF_GZ.read_bytes()),
    ):
        _utils.download_stream_to_zip(
            "http://example.com/a.gz",
            zip_path,
            "a.gz",
            chunk_size=512,
            overwrite=False,
            timeout=30,
            insecure_ssl=False,
            user_agent="test",
            attempts=1,
        )

    content_b = _EXON_GZ.read_bytes()
    with mock.patch(
        "recount3._utils.http_open",
        return_value=io.BytesIO(content_b),
    ):
        _utils.download_stream_to_zip(
            "http://example.com/b.gz",
            zip_path,
            "b.gz",
            chunk_size=512,
            overwrite=False,
            timeout=30,
            insecure_ssl=False,
            user_agent="test",
            attempts=1,
        )

    with zipfile.ZipFile(zip_path, "r") as zf:
        assert "b.gz" in zf.namelist()
        assert zf.read("b.gz") == content_b


def test_download_stream_to_zip_tmp_cleaned_up_on_failure(
    tmp_path: Path,
) -> None:
    """Temp download file is removed when the download fails."""
    with mock.patch(
        "recount3._utils.http_open",
        side_effect=urllib.error.URLError("refused"),
    ):
        with pytest.raises(urllib.error.URLError):
            _utils.download_stream_to_zip(
                "http://example.com/f.gz",
                tmp_path / "out.zip",
                "member.gz",
                chunk_size=512,
                overwrite=False,
                timeout=30,
                insecure_ssl=False,
                user_agent="test",
                attempts=1,
            )

    assert list(tmp_path.glob(".r3_dl_*.tmp")) == []


def test_download_stream_to_zip_tmp_unlink_error_swallowed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """OSError when deleting the tmp download file is swallowed."""
    zip_path = tmp_path / "out.zip"
    content = _GTF_GZ.read_bytes()

    original_unlink = Path.unlink

    def _raise_for_r3dl(
        self: Path, missing_ok: bool = False
    ) -> None:
        if ".r3_dl_" in self.name:
            raise OSError("cannot unlink r3 tmp")
        original_unlink(self, missing_ok=missing_ok)

    monkeypatch.setattr(Path, "unlink", _raise_for_r3dl)

    with mock.patch(
        "recount3._utils.http_open",
        return_value=io.BytesIO(content),
    ):
        _utils.download_stream_to_zip(
            "http://example.com/f.gz",
            zip_path,
            "gene.gz",
            chunk_size=512,
            overwrite=False,
            timeout=30,
            insecure_ssl=False,
            user_agent="test",
            attempts=1,
        )

    with zipfile.ZipFile(zip_path, "r") as zf:
        assert zf.read("gene.gz") == content


def test_download_stream_to_zip_invalid_zip_on_disk_raises(
    tmp_path: Path,
) -> None:
    """An invalid ZIP on disk surfaces a DownloadError."""
    zip_path = tmp_path / "out.zip"
    zip_path.write_bytes(b"not a zip")
    content = _GTF_GZ.read_bytes()

    with mock.patch(
        "recount3._utils.http_open",
        return_value=io.BytesIO(content),
    ):
        with pytest.raises(errors.DownloadError):
            _utils.download_stream_to_zip(
                "http://example.com/f.gz",
                zip_path,
                "member.gz",
                chunk_size=512,
                overwrite=False,
                timeout=30,
                insecure_ssl=False,
                user_agent="test",
                attempts=1,
            )


def test_write_cached_file_to_zip_creates_zip(tmp_path: Path) -> None:
    """write_cached_file_to_zip creates a new ZIP from a real file."""
    zip_path = tmp_path / "out.zip"
    _utils.write_cached_file_to_zip(
        _GTF_GZ, zip_path, "gene.gtf.gz", overwrite=False
    )
    with zipfile.ZipFile(zip_path, "r") as zf:
        assert "gene.gtf.gz" in zf.namelist()
        assert zf.read("gene.gtf.gz") == _GTF_GZ.read_bytes()


def test_write_cached_file_to_zip_delegates_to_write_or_replace(
    tmp_path: Path,
) -> None:
    """write_cached_file_to_zip is a thin wrapper over _write_or_replace."""
    zip_path = tmp_path / "out.zip"
    with mock.patch(
        "recount3._utils._write_or_replace_in_zip"
    ) as mock_inner:
        _utils.write_cached_file_to_zip(
            _GTF_GZ, zip_path, "gene.gz", overwrite=True
        )
    mock_inner.assert_called_once_with(
        zip_path, _GTF_GZ, "gene.gz", True
    )


@pytest.mark.parametrize(
    "value,expected",
    [
        ("gene", "gene"),
        ("exon", "exon"),
        ("junction", "junction"),
        ("GENE", "gene"),
        ("  Exon  ", "exon"),
        ("JUNCTION", "junction"),
    ],
)
def test_normalize_genomic_unit_valid_values(
    value: str, expected: str
) -> None:
    """Valid genomic units are returned in lowercase / stripped."""
    result = _utils._normalize_genomic_unit(  # pylint: disable=protected-access
        value
    )
    assert result == expected


@pytest.mark.parametrize(
    "value",
    ["transcript", "intron", "", "  ", "gene_sums"],
)
def test_normalize_genomic_unit_invalid_raises(value: str) -> None:
    """Invalid genomic units raise ValueError."""
    with pytest.raises(ValueError, match="Invalid genomic_unit"):
        _utils._normalize_genomic_unit(  # pylint: disable=protected-access
            value
        )


def test_resolve_counts_assay_name_preferred_present() -> None:
    """Preferred assay is returned when it is present."""
    se = mock.MagicMock()
    se.assay_names = ["raw_counts", "counts"]
    assert (
        _utils._resolve_counts_assay_name(se)  # pylint: disable=protected-access
        == "raw_counts"
    )


def test_resolve_counts_assay_name_falls_back_to_legacy(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Fallback assay is returned with a warning."""
    import logging

    se = mock.MagicMock()
    se.assay_names = ["counts"]
    with caplog.at_level(logging.WARNING):
        result = _utils._resolve_counts_assay_name(  # pylint: disable=protected-access
            se
        )
    assert result == "counts"
    assert "falling back" in caplog.text


def test_resolve_counts_assay_name_raises_when_neither_present() -> None:
    """ValueError is raised when neither assay name is available."""
    se = mock.MagicMock()
    se.assay_names = ["other_assay"]
    with pytest.raises(ValueError, match="raw_counts"):
        _utils._resolve_counts_assay_name(se)  # pylint: disable=protected-access


def test_resolve_counts_assay_name_raises_when_assay_names_absent() -> None:
    """ValueError when the object exposes no assay_names attribute."""
    se = mock.MagicMock(spec=[])  # No attributes at all.
    with pytest.raises(ValueError):
        _utils._resolve_counts_assay_name(se)  # pylint: disable=protected-access


def test_resolve_counts_assay_name_custom_names() -> None:
    """Custom preferred and fallback names are respected."""
    se = mock.MagicMock()
    se.assay_names = ["my_counts"]
    result = _utils._resolve_counts_assay_name(  # pylint: disable=protected-access
        se,
        preferred_assay_name="my_counts",
        fallback_assay_name="other",
    )
    assert result == "my_counts"


def test_coerce_col_data_to_pandas_passthrough_dataframe() -> None:
    """A pandas DataFrame is returned unchanged."""
    df = pd.DataFrame({"a": [1, 2]})
    result = _utils._coerce_col_data_to_pandas(  # pylint: disable=protected-access
        df
    )
    assert result is df


def test_coerce_col_data_to_pandas_se_like_object() -> None:
    """SE-like object with col_data.to_pandas() is coerced correctly."""
    df = pd.DataFrame({"sample": ["s1", "s2"]})
    se_like = mock.MagicMock()
    se_like.col_data.to_pandas.return_value = df
    result = _utils._coerce_col_data_to_pandas(  # pylint: disable=protected-access
        se_like
    )
    assert result is df


def test_coerce_col_data_to_pandas_raises_for_unsupported_type() -> None:
    """TypeError is raised for types that cannot be coerced."""
    with pytest.raises(TypeError, match="pandas.DataFrame"):
        _utils._coerce_col_data_to_pandas(  # pylint: disable=protected-access
            [1, 2, 3]
        )


def test_coerce_numeric_column_all_numeric() -> None:
    """A purely numeric Series is returned as float64."""
    s = pd.Series(["1.0", "2.5", "3"])
    result = _utils._coerce_numeric_column(  # pylint: disable=protected-access
        s, "score"
    )
    assert list(result) == pytest.approx([1.0, 2.5, 3.0])


def test_coerce_numeric_column_empty_strings_become_na() -> None:
    """Whitespace-only values are coerced to NA."""
    s = pd.Series(["1.0", "  ", "", "4.0"])
    result = _utils._coerce_numeric_column(  # pylint: disable=protected-access
        s, "score"
    )
    assert result.isna().sum() == 2
    assert result.dropna().tolist() == pytest.approx([1.0, 4.0])


def test_coerce_numeric_column_preserves_na_values() -> None:
    """Explicit NA / NaN values propagate through."""
    s = pd.Series([1.0, None, 3.0])
    result = _utils._coerce_numeric_column(  # pylint: disable=protected-access
        s, "score"
    )
    assert result.isna().sum() == 1


def test_coerce_numeric_column_raises_for_non_numeric_values() -> None:
    """Non-numeric, non-NA values raise ValueError."""
    s = pd.Series(["1.0", "abc", "3.0"])
    with pytest.raises(ValueError, match="non-numeric"):
        _utils._coerce_numeric_column(  # pylint: disable=protected-access
            s, "score"
        )


def test_resolve_metadata_column_exact_match() -> None:
    """An exact column name match returns the correct Series."""
    df = pd.DataFrame({"recount_qc.star.reads": [100, 200]})
    result = _utils._resolve_metadata_column(  # pylint: disable=protected-access
        df, "recount_qc.star.reads"
    )
    assert list(result) == [100, 200]


def test_resolve_metadata_column_case_insensitive_match() -> None:
    """Column lookup is case-insensitive."""
    df = pd.DataFrame({"RECOUNT_QC.Star.Reads": [1, 2]})
    result = _utils._resolve_metadata_column(  # pylint: disable=protected-access
        df, "recount_qc.star.reads"
    )
    assert list(result) == [1, 2]


def test_resolve_metadata_column_dot_to_dunder_fallback() -> None:
    """First '.' separator falls back to '__' when exact match fails."""
    df = pd.DataFrame({"recount_qc__star.reads": [10, 20]})
    result = _utils._resolve_metadata_column(  # pylint: disable=protected-access
        df, "recount_qc.star.reads"
    )
    assert list(result) == [10, 20]


def test_resolve_metadata_column_raises_when_not_found() -> None:
    """ValueError is raised when the column cannot be resolved."""
    df = pd.DataFrame({"other_col": [1]})
    with pytest.raises(ValueError, match="not found"):
        _utils._resolve_metadata_column(  # pylint: disable=protected-access
            df, "missing.column"
        )


def test_resolve_metadata_column_no_dot_raises_directly() -> None:
    """A column with no '.' and no match raises ValueError immediately."""
    df = pd.DataFrame({"a": [1]})
    with pytest.raises(ValueError, match="not found"):
        _utils._resolve_metadata_column(  # pylint: disable=protected-access
            df, "b"
        )


def test_format_optional_dep_error_known_module() -> None:
    """Known modules use a pre-configured install command."""
    msg = _utils._format_optional_dependency_import_error(  # pylint: disable=protected-access
        "biocframe"
    )
    assert "biocframe" in msg
    assert "pip install biocframe" in msg
    assert "Original import error" not in msg


def test_format_optional_dep_error_unknown_module() -> None:
    """Unknown modules fall back to 'pip install <name>'."""
    msg = _utils._format_optional_dependency_import_error(  # pylint: disable=protected-access
        "some_unknown_pkg"
    )
    assert "pip install some_unknown_pkg" in msg


def test_format_optional_dep_error_with_exception() -> None:
    """When exc is provided, the original error appears in the message."""
    exc = ImportError("no module named x")
    msg = _utils._format_optional_dependency_import_error(  # pylint: disable=protected-access
        "mymod", exc
    )
    assert "Original import error" in msg
    assert "no module named x" in msg


def test_format_optional_dep_error_without_exception() -> None:
    """When exc is None the error detail section is absent."""
    msg = _utils._format_optional_dependency_import_error(  # pylint: disable=protected-access
        "mymod"
    )
    assert "Original import error" not in msg


def test_format_optional_dep_failure_shows_exc_and_guidance() -> None:
    """Failure message includes the exception repr and install guidance."""
    exc = RuntimeError("shared library missing")
    msg = _utils._format_optional_dependency_import_failure(  # pylint: disable=protected-access
        "pyBigWig", exc
    )
    assert "could not be imported" in msg
    assert repr(exc) in msg
    assert "pip install" in msg


def test_import_optional_module_returns_real_module() -> None:
    """A genuinely installed module is returned."""
    import pandas

    mod = _utils.import_optional_module("pandas")
    assert mod is pandas


def test_import_optional_module_not_found_raises_import_error() -> None:
    """A missing module raises ImportError with install guidance."""
    with pytest.raises(ImportError, match="pip install"):
        _utils.import_optional_module(
            "__r3test_definitely_absent_xyz__"
        )


def test_import_optional_module_load_failure_raises_import_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A non-ModuleNotFoundError during import raises ImportError."""
    _utils.import_optional_module.cache_clear()
    monkeypatch.setattr(
        _utils.importlib,
        "import_module",
        mock.Mock(side_effect=RuntimeError("broken native ext")),
    )
    with pytest.raises(ImportError, match="could not be imported"):
        _utils.import_optional_module("__r3test_broken_ext__")
    _utils.import_optional_module.cache_clear()


def test_get_module_attribute_returns_existing_attribute() -> None:
    """An existing attribute is returned directly."""
    import pandas

    result = _utils._get_module_attribute(  # pylint: disable=protected-access
        pandas, "DataFrame", module_name="pandas"
    )
    assert result is pandas.DataFrame


def test_get_module_attribute_missing_raises_import_error() -> None:
    """A missing attribute raises ImportError with guidance."""
    fake_mod = types.ModuleType("fake_module")
    with pytest.raises(ImportError, match="pip install fake_module"):
        _utils._get_module_attribute(  # pylint: disable=protected-access
            fake_mod,
            "NonExistentClass",
            module_name="fake_module",
        )


def test_get_biocframe_class_returns_biocframe() -> None:
    """get_biocframe_class returns the BiocFrame class."""
    biocframe = pytest.importorskip("biocframe")
    cls = _utils.get_biocframe_class()
    assert cls is biocframe.BiocFrame


def test_get_biocframe_class_propagates_import_error() -> None:
    """ImportError from _get_module_attribute propagates unchanged."""
    with mock.patch(
        "recount3._utils._get_module_attribute",
        side_effect=ImportError("no BiocFrame"),
    ):
        with pytest.raises(ImportError, match="no BiocFrame"):
            _utils.get_biocframe_class()


def test_get_genomicranges_class_returns_genomicranges() -> None:
    """get_genomicranges_class returns the GenomicRanges class."""
    genomicranges = pytest.importorskip("genomicranges")
    cls = _utils.get_genomicranges_class()
    assert cls is genomicranges.GenomicRanges


def test_get_genomicranges_class_propagates_import_error() -> None:
    """ImportError from _get_module_attribute propagates."""
    with mock.patch(
        "recount3._utils._get_module_attribute",
        side_effect=ImportError("no GenomicRanges"),
    ):
        with pytest.raises(ImportError, match="no GenomicRanges"):
            _utils.get_genomicranges_class()


def test_get_summarizedexperiment_class_returns_class() -> None:
    """get_summarizedexperiment_class returns SummarizedExperiment."""
    se = pytest.importorskip("summarizedexperiment")
    cls = _utils.get_summarizedexperiment_class()
    assert cls is se.SummarizedExperiment


def test_get_summarizedexperiment_class_propagates_import_error() -> None:
    """ImportError from _get_module_attribute propagates."""
    with mock.patch(
        "recount3._utils._get_module_attribute",
        side_effect=ImportError("no SE"),
    ):
        with pytest.raises(ImportError, match="no SE"):
            _utils.get_summarizedexperiment_class()


def test_get_ranged_summarizedexperiment_class_returns_class() -> None:
    """Returns the RangedSummarizedExperiment class."""
    se = pytest.importorskip("summarizedexperiment")
    cls = _utils.get_ranged_summarizedexperiment_class()
    assert cls is se.RangedSummarizedExperiment


def test_get_ranged_summarizedexperiment_class_propagates_error() -> None:
    """ImportError from _get_module_attribute propagates."""
    with mock.patch(
        "recount3._utils._get_module_attribute",
        side_effect=ImportError("no RSE"),
    ):
        with pytest.raises(ImportError, match="no RSE"):
            _utils.get_ranged_summarizedexperiment_class()


def test_get_pybigwig_module_delegates_to_import_optional() -> None:
    """get_pybigwig_module calls import_optional_module('pyBigWig')."""
    fake_mod = types.ModuleType("pyBigWig")
    with mock.patch(
        "recount3._utils.import_optional_module",
        return_value=fake_mod,
    ) as mock_import:
        result = _utils.get_pybigwig_module()
    mock_import.assert_called_once_with("pyBigWig")
    assert result is fake_mod


def test_get_pybigwig_module_propagates_import_error() -> None:
    """ImportError is propagated when pyBigWig is absent."""
    with mock.patch(
        "recount3._utils.import_optional_module",
        side_effect=ImportError("pyBigWig not installed"),
    ):
        with pytest.raises(
            ImportError, match="pyBigWig not installed"
        ):
            _utils.get_pybigwig_module()


@pytest.mark.parametrize(
    "orig_ext,new_ext,expected_suffix",
    [
        ("MM", "ID", ".ID.gz"),
        ("MM", "RR", ".RR.gz"),
        ("ID", "MM", ".MM.gz"),
        ("ID", "RR", ".RR.gz"),
        ("RR", "MM", ".MM.gz"),
        ("RR", "ID", ".ID.gz"),
    ],
)
def test_derive_junction_sidecar_url_swaps_extension(
    orig_ext: str,
    new_ext: str,
    expected_suffix: str,
) -> None:
    """Extension suffix is swapped correctly for all valid combinations."""
    url = f"{_BASE_JXN_URL}.{orig_ext}.gz"
    result = _utils._derive_junction_sidecar_url(  # pylint: disable=protected-access
        url, new_ext
    )
    assert result.endswith(expected_suffix)
    assert result.startswith(_BASE_JXN_URL)


def test_derive_junction_sidecar_url_case_insensitive_source() -> None:
    """Lowercase extension in the source URL is matched."""
    url = f"{_BASE_JXN_URL}.mm.gz"
    result = _utils._derive_junction_sidecar_url(  # pylint: disable=protected-access
        url, "ID"
    )
    assert result.endswith(".ID.gz")


def test_derive_junction_sidecar_url_new_ext_uppercased() -> None:
    """new_ext is uppercased regardless of input case."""
    url = f"{_BASE_JXN_URL}.MM.gz"
    result = _utils._derive_junction_sidecar_url(  # pylint: disable=protected-access
        url, "id"
    )
    assert result.endswith(".ID.gz")


def test_derive_junction_sidecar_url_empty_url_raises() -> None:
    """Empty URL raises ValueError."""
    with pytest.raises(ValueError, match="non-empty"):
        _utils._derive_junction_sidecar_url(  # pylint: disable=protected-access
            "", "MM"
        )


def test_derive_junction_sidecar_url_invalid_ext_raises() -> None:
    """new_ext not in {MM, ID, RR} raises ValueError."""
    with pytest.raises(ValueError, match="MM/ID/RR"):
        _utils._derive_junction_sidecar_url(  # pylint: disable=protected-access
            f"{_BASE_JXN_URL}.MM.gz", "XY"
        )


def test_derive_junction_sidecar_url_no_pattern_raises() -> None:
    """URL without a recognised .{MM,ID,RR}.gz suffix raises ValueError."""
    with pytest.raises(ValueError, match="Cannot derive"):
        _utils._derive_junction_sidecar_url(  # pylint: disable=protected-access
            "https://example.com/file.tsv.gz", "MM"
        )
