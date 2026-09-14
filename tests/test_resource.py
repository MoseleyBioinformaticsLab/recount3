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

# pylint: disable=redefined-outer-name

from __future__ import annotations

import contextlib
import dataclasses
import gzip
import hashlib
import inspect
import io
import shutil
import threading
import warnings
from concurrent import futures
from pathlib import Path
from unittest import mock

import pandas as pd
import pytest
import scipy.io
import scipy.sparse

import recount3.resource as res_module
from recount3._bigwig import BigWigFile
from recount3._descriptions import R3ResourceDescription
from recount3._utils import (
    _cache_path,
    _derive_junction_sidecar_url,
    _path_lock_for_path,
    import_optional_module,
)
from recount3.config import Config, recount3_cache_files, recount3_cache_rm
from recount3.errors import DownloadError, LoadError
from recount3.resource import (
    R3Resource,
    _detect_mmread_kwargs,
    _ensure_cached_url,
    _read_id_rail_ids,
    _read_mm_matrix,
)

_DATA_DIR = Path(__file__).parent / "data"
_MIRROR = _DATA_DIR / "recount3_mirror" / "recount3"

_GENE_GZ = (
    _MIRROR
    / "human"
    / "data_sources"
    / "sra"
    / "gene_sums"
    / "65"
    / "SRP014565"
    / "sra.gene_sums.SRP014565.G026.gz"
)
_JXN_MM_GZ = (
    _MIRROR
    / "human"
    / "data_sources"
    / "sra"
    / "junctions"
    / "65"
    / "SRP014565"
    / "sra.junctions.SRP014565.ALL.MM.gz"
)
_JXN_ID_GZ = (
    _MIRROR
    / "human"
    / "data_sources"
    / "sra"
    / "junctions"
    / "65"
    / "SRP014565"
    / "sra.junctions.SRP014565.ALL.ID.gz"
)
_JXN_RR_GZ = (
    _MIRROR
    / "human"
    / "data_sources"
    / "sra"
    / "junctions"
    / "65"
    / "SRP014565"
    / "sra.junctions.SRP014565.ALL.RR.gz"
)
_META_MD_GZ = (
    _MIRROR
    / "human"
    / "data_sources"
    / "sra"
    / "metadata"
    / "65"
    / "SRP014565"
    / "sra.recount_qc.SRP014565.MD.gz"
)
_BW_PATH = (
    _MIRROR
    / "human"
    / "data_sources"
    / "sra"
    / "base_sums"
    / "15"
    / "SRP009615"
    / "77"
    / "sra.base_sums.SRP009615_SRR387777.ALL.bw"
)

# Stand-in remote payload for the cache tests, which never reach the network.
_CACHE_URL = "https://example.org/a.gz"


@pytest.fixture()
def cfg(tmp_path: Path) -> Config:
    """Minimal Config pointing at a temporary cache directory."""
    return Config(
        base_url="https://example.org/recount3/",
        timeout=30,
        insecure_ssl=False,
        max_retries=1,
        user_agent="test/1.0",
        cache_dir=tmp_path / "cache",
        cache_disabled=False,
        chunk_size=1024,
    )


def _make(resource_type: str, cfg: Config, **kw: object) -> R3Resource:
    """Construct an R3Resource via a registered description."""
    desc = R3ResourceDescription(resource_type=resource_type, **kw)
    return R3Resource(description=desc, config=cfg)


def _seed(resource: R3Resource, src: Path) -> Path:
    """Copy src to the cache location expected for resource."""
    dest = resource._cached_path()
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    return dest


def _fake_cfg(tmp_path: Path) -> Config:
    """Minimal Config for default_config() mock targets."""
    return Config(
        base_url="https://fallback.org/recount3/",
        timeout=10,
        insecure_ssl=False,
        max_retries=1,
        user_agent="u/1",
        cache_dir=tmp_path / "dc",
        cache_disabled=False,
        chunk_size=512,
    )


# ===========================================================================
# _ensure_cached_url
# ===========================================================================


def test_ensure_cached_url_outer_cache_hit(cfg: Config) -> None:
    """Returns immediately when the cached file already exists."""
    url = "https://example.org/recount3/human/thing.tsv.gz"
    cache_p = _cache_path(url, cfg.cache_dir)
    cache_p.parent.mkdir(parents=True, exist_ok=True)
    cache_p.touch()
    with mock.patch("recount3.resource.download_to_file") as mock_dl:
        result = _ensure_cached_url(
            url=url,
            cache_root=cfg.cache_dir,
            cfg=cfg,
            chunk_size=1024,
        )
    assert result == cache_p
    mock_dl.assert_not_called()


def test_ensure_cached_url_inner_lock_check_hit(cfg: Config) -> None:
    """Inner lock check returns early when another thread pre-created file."""
    url = "https://example.org/recount3/human/thing2.tsv.gz"
    cache_p = _cache_path(url, cfg.cache_dir)
    cache_p.parent.mkdir(parents=True, exist_ok=True)

    # Simulate another thread creating the file between outer and inner checks.
    class _SeedOnEnter:
        def __enter__(self) -> "_SeedOnEnter":
            cache_p.touch()
            return self

        def __exit__(self, *args: object) -> bool:
            return False

    with mock.patch.object(
        res_module, "_path_lock_for_path", return_value=_SeedOnEnter()
    ):
        with mock.patch("recount3.resource.download_to_file") as mock_dl:
            result = _ensure_cached_url(
                url=url,
                cache_root=cfg.cache_dir,
                cfg=cfg,
                chunk_size=1024,
            )
    assert result == cache_p
    mock_dl.assert_not_called()


def test_ensure_cached_url_downloads_missing_file(cfg: Config) -> None:
    """Calls download_to_file when file is absent on both checks."""
    url = "https://example.org/recount3/human/thing3.tsv.gz"
    with mock.patch("recount3.resource.download_to_file") as mock_dl:
        result = _ensure_cached_url(
            url=url,
            cache_root=cfg.cache_dir,
            cfg=cfg,
            chunk_size=512,
        )
    mock_dl.assert_called_once()
    assert mock_dl.call_args.args[0] == url
    assert "thing3.tsv.gz" in result.name


# ===========================================================================
# _read_id_rail_ids
# ===========================================================================


def test_read_id_rail_ids_rail_id_column() -> None:
    """Reads rail_id values from a real compressed junction ID file."""
    ids = _read_id_rail_ids(_JXN_ID_GZ)
    assert isinstance(ids, list)
    assert len(ids) == 1
    assert ids[0] == "985131"


def test_read_id_rail_ids_no_rail_id_uses_first_column(
    tmp_path: Path,
) -> None:
    """Uses the first column when no 'rail_id' header is present."""
    p = tmp_path / "noid.tsv.gz"
    df = pd.DataFrame({"sample_id": ["S1", "S2"]})
    with gzip.open(p, "wt") as fh:
        df.to_csv(fh, sep="\t", index=False)
    ids = _read_id_rail_ids(p)
    assert ids == ["S1", "S2"]


def test_read_id_rail_ids_canonicalizes_widened_numeric_column(
    tmp_path: Path,
) -> None:
    """A blank row must not turn the other rail IDs into "123488.0".

    These strings become the junction matrix's column labels and are
    matched against the rail IDs in the sample metadata, so a float
    rendering here would match nothing and drop every sample.
    """
    p = tmp_path / "ids.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("rail_id\tsample\n123488\ta\n\tb\n123474\tc\n")
    assert _read_id_rail_ids(p) == ["123488", "<NA>", "123474"]


def test_read_id_rail_ids_empty_dataframe_raises(tmp_path: Path) -> None:
    """Raises LoadError when the parsed DataFrame is empty."""
    p = tmp_path / "empty.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("rail_id\n")  # header only, no rows
    with pytest.raises(LoadError, match="parsed empty"):
        _read_id_rail_ids(p)


def test_read_id_rail_ids_empty_list_raises_load_error(
    tmp_path: Path,
) -> None:
    """LoadError when rail_ids list is empty despite a non-empty DataFrame.

    This covers the empty-rail-ID guard, which is unreachable via normal
    I/O (a non-empty DataFrame always yields >= 1 items). It is reached by
    mocking pd.read_csv to return a DataFrame that passes df.empty, and by
    mocking the canonicalizer it is passed through to return a column whose
    tolist() is empty.
    """
    p = tmp_path / "data.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("rail_id\nsome_value\n")

    mock_canonical = mock.MagicMock()
    mock_canonical.astype.return_value.tolist.return_value = []

    mock_df = mock.MagicMock()
    mock_df.empty = False
    mock_df.columns = ["rail_id"]
    mock_df.__getitem__ = mock.Mock(return_value=mock.MagicMock())

    with (
        mock.patch("recount3.resource.pd.read_csv", return_value=mock_df),
        mock.patch(
            "recount3.resource.canonical_identifier_series",
            return_value=mock_canonical,
        ),
    ):
        with pytest.raises(LoadError, match="has no rail IDs"):
            _read_id_rail_ids(p)


def test_read_id_rail_ids_python_engine_fallback(tmp_path: Path) -> None:
    """Python engine is used as fallback when C engine raises."""
    p = tmp_path / "data.tsv.gz"
    df_src = pd.DataFrame({"rail_id": ["999"]})
    with gzip.open(p, "wt") as fh:
        df_src.to_csv(fh, sep="\t", index=False)

    original = pd.read_csv
    call_count = 0

    def _patched_read_csv(*args: object, **kwargs: object) -> pd.DataFrame:
        nonlocal call_count
        call_count += 1
        if call_count == 1:
            raise pd.errors.ParserError("simulated C engine failure")
        return original(*args, **kwargs)  # type: ignore[arg-type]

    with mock.patch("recount3.resource.pd.read_csv", _patched_read_csv):
        ids = _read_id_rail_ids(p)

    assert ids == ["999"]
    assert call_count == 2


# ===========================================================================
# _read_mm_matrix
# ===========================================================================


def test_read_mm_matrix_gz_file() -> None:
    """Reads a .MM.gz file into a CSR sparse array."""
    mat = _read_mm_matrix(_JXN_MM_GZ)
    assert isinstance(mat, scipy.sparse.csr_array)
    assert mat.ndim == 2


def test_read_mm_matrix_plain_file(tmp_path: Path) -> None:
    """Reads a plain (non-compressed) MatrixMarket file."""
    # scipy.io.mmwrite appends .mtx when the name has no .mtx suffix;
    # use .mtx directly so the written file matches the path we read.
    mm_file = tmp_path / "test.mtx"
    scipy.io.mmwrite(str(mm_file), scipy.sparse.eye(4, k=1, format="coo"))
    mat = _read_mm_matrix(mm_file)
    assert isinstance(mat, scipy.sparse.csr_array)
    assert mat.shape == (4, 4)


def test_read_mm_matrix_emits_no_deprecation_warning() -> None:
    """Reading a MatrixMarket file is free of SciPy deprecation warnings."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        mat = _read_mm_matrix(_JXN_MM_GZ)
    assert isinstance(mat, scipy.sparse.csr_array)


def test_read_mm_matrix_without_spmatrix_support(tmp_path: Path) -> None:
    """Older SciPy releases lack ``spmatrix=`` and must still be supported."""
    mm_file = tmp_path / "legacy.mtx"
    scipy.io.mmwrite(str(mm_file), scipy.sparse.eye(3, format="coo"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        with mock.patch.object(res_module, "_MMREAD_KWARGS", {}):
            mat = _read_mm_matrix(mm_file)
    assert isinstance(mat, scipy.sparse.csr_array)
    assert mat.shape == (3, 3)


def test_detect_mmread_kwargs_matches_installed_scipy() -> None:
    """Keyword detection agrees with the installed SciPy signature."""
    supported = "spmatrix" in inspect.signature(scipy.io.mmread).parameters
    expected = {"spmatrix": False} if supported else {}
    assert _detect_mmread_kwargs() == expected


def test_detect_mmread_kwargs_handles_unintrospectable_mmread() -> None:
    """Detection degrades to no keywords when the signature is unavailable."""
    with mock.patch.object(
        inspect, "signature", side_effect=ValueError("no signature")
    ):
        assert _detect_mmread_kwargs() == {}


def test_read_mm_matrix_bad_file_raises_load_error(tmp_path: Path) -> None:
    """LoadError is raised when the file cannot be parsed."""
    bad = tmp_path / "bad.MM.gz"
    with gzip.open(bad, "wb") as fh:
        fh.write(b"not a matrix market file")
    with pytest.raises(LoadError, match="Failed to read MatrixMarket"):
        _read_mm_matrix(bad)


def test_read_mm_matrix_non_2d_raises_load_error(tmp_path: Path) -> None:
    """LoadError when parsed object reports ndim != 2."""

    class _FakeMat:
        ndim = 1

    plain = tmp_path / "test2.mtx"
    scipy.io.mmwrite(str(plain), scipy.sparse.eye(2, format="coo"))
    with mock.patch.object(scipy.io, "mmread", return_value=_FakeMat()):
        with pytest.raises(LoadError, match="not 2-dimensional"):
            _read_mm_matrix(plain)


# ===========================================================================
# R3Resource — construction
# ===========================================================================


def test_post_init_derives_url_from_description(cfg: Config) -> None:
    """url=None causes __post_init__ to derive the URL from description."""
    res = _make("data_sources", cfg, organism="human")
    assert res.url is not None
    assert res.url.startswith("https://example.org/recount3/")
    assert "human/homes_index" in res.url


def test_post_init_preserves_explicit_url(cfg: Config) -> None:
    """An explicitly-provided URL is kept unchanged."""
    desc = R3ResourceDescription(resource_type="data_sources", organism="human")
    res = R3Resource(
        description=desc,
        url="https://custom.host/file.gz",
        config=cfg,
    )
    assert res.url == "https://custom.host/file.gz"


def test_post_init_config_none_calls_default_config(
    tmp_path: Path,
) -> None:
    """default_config() is called when config=None in __post_init__."""
    desc = R3ResourceDescription(resource_type="data_sources", organism="human")
    fc = _fake_cfg(tmp_path)
    with mock.patch("recount3.resource.default_config", return_value=fc):
        res = R3Resource(description=desc)
    assert res.url is not None
    assert res.url.startswith("https://fallback.org/recount3/")


def test_r3resource_defaults_are_none(cfg: Config) -> None:
    """Filepath and _cached_data default to None on construction."""
    res = _make("data_sources", cfg, organism="human")
    assert res.filepath is None
    assert res.get_loaded() is None


# ===========================================================================
# R3Resource — arcname / cache helpers
# ===========================================================================


def test_arcname_has_no_leading_slash(cfg: Config) -> None:
    """Arcname never starts with a slash."""
    res = _make("data_sources", cfg, organism="human")
    assert not res.arcname.startswith("/")
    assert res.arcname == "human/homes_index"


def test_cache_root_uses_config(cfg: Config) -> None:
    """_cache_root() returns cfg.cache_dir when config is set."""
    res = _make("data_sources", cfg, organism="human")
    assert res._cache_root() == cfg.cache_dir


def test_cache_root_calls_default_config(tmp_path: Path) -> None:
    """_cache_root() calls default_config() when config=None."""
    desc = R3ResourceDescription(resource_type="data_sources", organism="human")
    fc = _fake_cfg(tmp_path)
    with mock.patch("recount3.resource.default_config", return_value=fc):
        res = R3Resource(description=desc)
        root = res._cache_root()
    assert root == tmp_path / "dc"


def test_cached_path_filename_contains_url_basename(cfg: Config) -> None:
    """_cached_path() filename includes the URL basename component."""
    res = _make(
        "count_files_gene_or_exon",
        cfg,
        organism="human",
        data_source="sra",
        genomic_unit="gene",
        project="SRP014565",
        annotation_extension="G026",
    )
    cp = res._cached_path()
    assert "sra.gene_sums.SRP014565.G026.gz" in cp.name


class TestResourceEquality:
    """Equality asks which file a resource is, not what it has parsed.

    ``_cached_data`` holds whatever ``load()`` produced. For counts and
    metadata that is a ``pandas.DataFrame``, whose ``==`` is element-wise,
    so including it in the comparison made equality raise
    ``ValueError: The truth value of a DataFrame is ambiguous`` for two
    loaded resources naming the same file.
    """

    @staticmethod
    def _resource(project: str = "P1") -> R3Resource:
        return R3Resource(
            R3ResourceDescription(
                resource_type="count_files_gene_or_exon",
                organism="human",
                data_source="sra",
                project=project,
                genomic_unit="gene",
                annotation_extension="G026",
            )
        )

    def test_two_loaded_resources_for_the_same_file_compare_equal(
        self,
    ) -> None:
        first, second = self._resource(), self._resource()
        first._cached_data = pd.DataFrame({"x": [1, 2]})
        second._cached_data = pd.DataFrame({"x": [3, 4]})
        assert first == second
        assert first in [second]

    def test_loading_does_not_change_equality(self) -> None:
        first, second = self._resource(), self._resource()
        assert first == second
        first._cached_data = pd.DataFrame({"x": [1]})
        assert first == second

    def test_different_files_are_still_unequal(self) -> None:
        assert self._resource("P1") != self._resource("P2")

    def test_a_materialized_path_still_distinguishes_resources(self) -> None:
        first, second = self._resource(), self._resource()
        first.filepath = "/tmp/one.gz"
        assert first != second


class TestEnsureCached:
    """`_cached_path()` computes a path; `ensure_cached()` guarantees a file.

    The distinction matters because callers that read the cached file
    themselves cannot tell a cache miss from a hit by calling
    `_cached_path()`, it returns a path either way and never raises.
    """

    @staticmethod
    def _resource(cfg: Config) -> R3Resource:
        return _make(
            "count_files_gene_or_exon",
            cfg,
            organism="human",
            data_source="sra",
            genomic_unit="gene",
            project="SRP014565",
            annotation_extension="G026",
        )

    def test_returns_the_cached_file_without_downloading(
        self, cfg: Config, tmp_path: Path
    ) -> None:
        res = self._resource(cfg)
        src = tmp_path / "seed.gz"
        src.write_bytes(b"payload")
        seeded = _seed(res, src)

        with mock.patch.object(res_module, "download_to_file") as dl:
            assert res.ensure_cached() == seeded
        dl.assert_not_called()

    def test_downloads_when_the_file_is_absent(self, cfg: Config) -> None:
        res = self._resource(cfg)
        assert not res._cached_path().exists()

        def fake_download(url: str, out_path: Path, **_kw: object) -> None:
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_bytes(b"fetched")

        with mock.patch.object(
            res_module, "download_to_file", side_effect=fake_download
        ) as dl:
            path = res.ensure_cached()
        dl.assert_called_once()
        assert path.read_bytes() == b"fetched"

    def test_raises_instead_of_downloading_when_disabled(
        self, cfg: Config
    ) -> None:
        res = self._resource(cfg)
        with mock.patch.object(res_module, "download_to_file") as dl:
            with pytest.raises(FileNotFoundError):
                res.ensure_cached(download=False)
        dl.assert_not_called()

    def test_a_resource_with_no_usable_path_still_tries_to_download(
        self, cfg: Config, tmp_path: Path
    ) -> None:
        """The only case where `_cached_path()` itself can fail."""
        res = self._resource(cfg)
        real_path = res._cached_path()
        calls = {"n": 0}

        def flaky_cached_path(_self: R3Resource) -> Path:
            calls["n"] += 1
            if calls["n"] == 1:
                raise RuntimeError("malformed resource")
            return real_path

        def fake_download(url: str, out_path: Path, **_kw: object) -> None:
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_bytes(b"fetched")

        with (
            mock.patch.object(R3Resource, "_cached_path", flaky_cached_path),
            mock.patch.object(
                res_module, "download_to_file", side_effect=fake_download
            ),
        ):
            assert res.ensure_cached().read_bytes() == b"fetched"


# ===========================================================================
# R3Resource._ensure_cached
# ===========================================================================


def test_ensure_cached_disable_raises(cfg: Config) -> None:
    """_ensure_cached('disable') raises ValueError."""
    res = _make("data_sources", cfg, organism="human")
    with pytest.raises(ValueError, match="disable.*mode"):
        res._ensure_cached(mode="disable", chunk_size=1024)


def test_ensure_cached_enable_calls_ensure_cached_url(
    cfg: Config,
) -> None:
    """_ensure_cached('enable') delegates to _ensure_cached_url."""
    res = _make("data_sources", cfg, organism="human")
    with mock.patch(
        "recount3.resource._ensure_cached_url",
        return_value=Path("/fake"),
    ) as mock_ecu:
        result = res._ensure_cached(mode="enable", chunk_size=512)
    mock_ecu.assert_called_once()
    assert result == Path("/fake")


def test_ensure_cached_update_calls_download_to_file(
    cfg: Config,
) -> None:
    """_ensure_cached('update') calls download_to_file under the lock."""
    res = _make("data_sources", cfg, organism="human")
    with mock.patch("recount3.resource.download_to_file") as mock_dl:
        res._ensure_cached(mode="update", chunk_size=256)
    mock_dl.assert_called_once()
    assert mock_dl.call_args.kwargs["chunk_size"] == 256


def test_ensure_cached_unknown_mode_raises(cfg: Config) -> None:
    """_ensure_cached raises ValueError for an unknown mode string."""
    res = _make("data_sources", cfg, organism="human")
    with pytest.raises(ValueError, match="Unknown cache mode"):
        res._ensure_cached(
            mode="bogus",  # type: ignore[arg-type]
            chunk_size=1024,
        )


def test_ensure_cached_config_none_calls_default_config(
    tmp_path: Path,
) -> None:
    """_ensure_cached calls default_config() when config=None."""
    desc = R3ResourceDescription(resource_type="data_sources", organism="human")
    fc = _fake_cfg(tmp_path)
    with mock.patch("recount3.resource.default_config", return_value=fc):
        res = R3Resource(description=desc)
        with mock.patch("recount3.resource.download_to_file"):
            res._ensure_cached(mode="update", chunk_size=512)


# ===========================================================================
# R3Resource.download — validation / cache-only
# ===========================================================================


def test_download_invalid_cache_mode_raises(cfg: Config) -> None:
    """download() raises ValueError for an unrecognised cache_mode."""
    res = _make("data_sources", cfg, organism="human")
    with pytest.raises(ValueError, match="Invalid cache_mode"):
        res.download(cache_mode="bogus")  # type: ignore[arg-type]


def test_download_cache_disabled_forces_disable_mode(
    cfg: Config, tmp_path: Path
) -> None:
    """When cfg.cache_disabled=True, cache_mode is overridden to 'disable'."""
    disabled_cfg = Config(
        base_url=cfg.base_url,
        timeout=cfg.timeout,
        insecure_ssl=cfg.insecure_ssl,
        max_retries=cfg.max_retries,
        user_agent=cfg.user_agent,
        cache_dir=cfg.cache_dir,
        cache_disabled=True,
    )
    desc = R3ResourceDescription(resource_type="data_sources", organism="human")
    res = R3Resource(description=desc, config=disabled_cfg)
    dest_dir = tmp_path / "out"
    dest_dir.mkdir()
    with mock.patch("recount3.resource.download_to_file") as mock_dl:
        res.download(str(dest_dir))
    mock_dl.assert_called_once()


def test_download_config_none_calls_default_config(
    tmp_path: Path,
) -> None:
    """download() calls default_config() when config=None."""
    desc = R3ResourceDescription(resource_type="data_sources", organism="human")
    fc = _fake_cfg(tmp_path)
    (tmp_path / "out").mkdir()
    with mock.patch("recount3.resource.default_config", return_value=fc):
        res = R3Resource(description=desc)
        with mock.patch("recount3.resource.download_to_file"):
            with mock.patch("recount3.resource._hardlink_or_copy"):
                res.download(str(tmp_path / "out"))


def test_download_path_none_disable_raises(cfg: Config) -> None:
    """path=None with cache_mode='disable' raises ValueError."""
    res = _make("data_sources", cfg, organism="human")
    with pytest.raises(ValueError, match="cache-only requires cache enabled"):
        res.download(path=None, cache_mode="disable")


def test_download_path_none_enable_returns_none(cfg: Config) -> None:
    """path=None with cache_mode='enable' returns None."""
    res = _make("data_sources", cfg, organism="human")
    with mock.patch.object(
        R3Resource, "_ensure_cached", return_value=Path("/fake")
    ):
        result = res.download(path=None, cache_mode="enable")
    assert result is None


def test_download_path_none_update_returns_none(cfg: Config) -> None:
    """path=None with cache_mode='update' returns None."""
    res = _make("data_sources", cfg, organism="human")
    with mock.patch.object(
        R3Resource, "_ensure_cached", return_value=Path("/fake")
    ):
        result = res.download(path=None, cache_mode="update")
    assert result is None


def test_download_unsupported_extension_raises(
    cfg: Config, tmp_path: Path
) -> None:
    """Path with non-directory, non-.zip extension raises ValueError."""
    res = _make("data_sources", cfg, organism="human")
    with pytest.raises(ValueError, match="directory or a .zip"):
        res.download(str(tmp_path / "out.tar.gz"))


def test_download_chunk_size_explicit_value_used(
    cfg: Config, tmp_path: Path
) -> None:
    """Explicit chunk_size is forwarded; config default is not used."""
    res = _make("data_sources", cfg, organism="human")
    (tmp_path / "out").mkdir()
    fake = tmp_path / "cf"
    fake.touch()
    with mock.patch.object(
        R3Resource, "_ensure_cached", return_value=fake
    ) as mock_ec:
        with mock.patch("recount3.resource._hardlink_or_copy"):
            res.download(str(tmp_path / "out"), chunk_size=128)
    assert mock_ec.call_args.kwargs["chunk_size"] == 128


def test_download_chunk_size_defaults_to_config(
    cfg: Config, tmp_path: Path
) -> None:
    """When chunk_size=None, cfg.chunk_size is forwarded."""
    res = _make("data_sources", cfg, organism="human")
    (tmp_path / "out").mkdir()
    fake = tmp_path / "cf"
    fake.touch()
    with mock.patch.object(
        R3Resource, "_ensure_cached", return_value=fake
    ) as mock_ec:
        with mock.patch("recount3.resource._hardlink_or_copy"):
            res.download(str(tmp_path / "out"))
    assert mock_ec.call_args.kwargs["chunk_size"] == cfg.chunk_size


# ===========================================================================
# R3Resource.download — ZIP archive materialization
# ===========================================================================


def test_download_to_zip_disable_streams_directly(
    cfg: Config, tmp_path: Path
) -> None:
    """ZIP + disable: calls download_stream_to_zip without caching."""
    res = _make("data_sources", cfg, organism="human")
    zip_path = tmp_path / "out.zip"
    with mock.patch("recount3.resource.download_stream_to_zip") as mock_ds:
        result = res.download(str(zip_path), cache_mode="disable")
    mock_ds.assert_called_once()
    assert result is None


def test_download_to_zip_enable_caches_then_writes(
    cfg: Config, tmp_path: Path
) -> None:
    """ZIP + enable: caches the file, then calls write_cached_file_to_zip."""
    res = _make("data_sources", cfg, organism="human")
    zip_path = tmp_path / "out.zip"
    fake_cached = tmp_path / "cached_file"
    fake_cached.touch()
    with mock.patch.object(
        R3Resource, "_ensure_cached", return_value=fake_cached
    ):
        with mock.patch(
            "recount3.resource.write_cached_file_to_zip"
        ) as mock_wc:
            result = res.download(
                str(zip_path), cache_mode="enable", overwrite=True
            )
    mock_wc.assert_called_once_with(
        fake_cached, zip_path, res.arcname, overwrite=True
    )
    assert result is None


def test_download_to_zip_update_caches_then_writes(
    cfg: Config, tmp_path: Path
) -> None:
    """ZIP + update: forced cache refresh, then write_cached_file_to_zip."""
    res = _make("data_sources", cfg, organism="human")
    zip_path = tmp_path / "out.zip"
    fake_cached = tmp_path / "cached_file"
    fake_cached.touch()
    with mock.patch.object(
        R3Resource, "_ensure_cached", return_value=fake_cached
    ):
        with mock.patch(
            "recount3.resource.write_cached_file_to_zip"
        ) as mock_wc:
            result = res.download(str(zip_path), cache_mode="update")
    mock_wc.assert_called_once_with(
        fake_cached, zip_path, res.arcname, overwrite=False
    )
    assert result is None


# ===========================================================================
# R3Resource.download — directory materialization
# ===========================================================================


def test_download_to_dir_disable_calls_download_to_file(
    cfg: Config, tmp_path: Path
) -> None:
    """Directory + disable: streams directly via download_to_file."""
    res = _make("data_sources", cfg, organism="human")
    dest_dir = tmp_path / "out"
    dest_dir.mkdir()
    with mock.patch("recount3.resource.download_to_file") as mock_dl:
        result = res.download(str(dest_dir), cache_mode="disable")
    mock_dl.assert_called_once()
    assert result is not None
    assert result.endswith("homes_index")
    assert res.filepath == result


def test_download_to_dir_enable_hardlinks_from_cache(
    cfg: Config, tmp_path: Path
) -> None:
    """Directory + enable: hardlinks cached file, sets filepath."""
    res = _make("data_sources", cfg, organism="human")
    dest_dir = tmp_path / "out"
    dest_dir.mkdir()
    fake_cached = tmp_path / "cf"
    fake_cached.touch()
    with mock.patch.object(
        R3Resource, "_ensure_cached", return_value=fake_cached
    ):
        with mock.patch("recount3.resource._hardlink_or_copy") as mock_hl:
            result = res.download(str(dest_dir), cache_mode="enable")
    mock_hl.assert_called_once()
    assert result is not None
    assert res.filepath == result


def test_download_to_dir_dest_exists_no_overwrite_returns_early(
    cfg: Config, tmp_path: Path
) -> None:
    """Destination exists + overwrite=False: skips hardlink, returns path."""
    res = _make("data_sources", cfg, organism="human")
    dest_dir = tmp_path / "out"
    dest_dir.mkdir()
    (dest_dir / "homes_index").touch()
    fake_cached = tmp_path / "cf"
    fake_cached.touch()
    with mock.patch.object(
        R3Resource, "_ensure_cached", return_value=fake_cached
    ):
        with mock.patch("recount3.resource._hardlink_or_copy") as mock_hl:
            result = res.download(
                str(dest_dir), cache_mode="enable", overwrite=False
            )
    mock_hl.assert_not_called()
    assert result == str(dest_dir / "homes_index")
    assert res.filepath == result


def test_download_to_dir_dest_exists_overwrite_relinks(
    cfg: Config, tmp_path: Path
) -> None:
    """Destination exists + overwrite=True: hardlinks anyway."""
    res = _make("data_sources", cfg, organism="human")
    dest_dir = tmp_path / "out"
    dest_dir.mkdir()
    (dest_dir / "homes_index").touch()
    fake_cached = tmp_path / "cf"
    fake_cached.touch()
    with mock.patch.object(
        R3Resource, "_ensure_cached", return_value=fake_cached
    ):
        with mock.patch("recount3.resource._hardlink_or_copy") as mock_hl:
            result = res.download(
                str(dest_dir), cache_mode="enable", overwrite=True
            )
    mock_hl.assert_called_once()
    assert result == str(dest_dir / "homes_index")


def test_download_to_dir_with_suffix_uses_is_dir(
    cfg: Config, tmp_path: Path
) -> None:
    """A path that is a directory but has a suffix is treated as a directory."""
    # Exercises the path_p.is_dir() branch in: suffix == "" or is_dir()
    res = _make("data_sources", cfg, organism="human")
    dest_dir = tmp_path / "out.d"  # non-empty suffix but IS a directory
    dest_dir.mkdir()
    fake_cached = tmp_path / "cf"
    fake_cached.touch()
    with mock.patch.object(
        R3Resource, "_ensure_cached", return_value=fake_cached
    ):
        with mock.patch("recount3.resource._hardlink_or_copy"):
            result = res.download(str(dest_dir), cache_mode="enable")
    assert result is not None


# ===========================================================================
# R3Resource.download — os.PathLike destinations
# ===========================================================================


def test_download_to_dir_accepts_pathlike(cfg: Config, tmp_path: Path) -> None:
    """A bare Path directory works like the equivalent str."""
    res = _make("data_sources", cfg, organism="human")
    dest_dir = tmp_path / "out"
    dest_dir.mkdir()
    fake_cached = tmp_path / "cf"
    fake_cached.touch()
    with mock.patch.object(
        R3Resource, "_ensure_cached", return_value=fake_cached
    ):
        with mock.patch("recount3.resource._hardlink_or_copy") as mock_hl:
            result = res.download(dest_dir, cache_mode="enable")
    mock_hl.assert_called_once()
    assert result == str(dest_dir / "homes_index")
    assert isinstance(res.filepath, str)


def test_download_to_zip_accepts_pathlike(cfg: Config, tmp_path: Path) -> None:
    """A bare Path .zip destination works like the equivalent str."""
    res = _make("data_sources", cfg, organism="human")
    zip_path = tmp_path / "out.zip"
    fake_cached = tmp_path / "cached_file"
    fake_cached.touch()
    with mock.patch.object(
        R3Resource, "_ensure_cached", return_value=fake_cached
    ):
        with mock.patch(
            "recount3.resource.write_cached_file_to_zip"
        ) as mock_wc:
            result = res.download(zip_path, cache_mode="enable")
    mock_wc.assert_called_once_with(
        fake_cached, zip_path, res.arcname, overwrite=False
    )
    assert result is None


def test_download_accepts_pathlike_from_the_packages_own_output(
    cfg: Config, tmp_path: Path
) -> None:
    """A Path produced by the package feeds back in without conversion."""
    res = _make("data_sources", cfg, organism="human")
    cached = tmp_path / "cf"
    cached.touch()
    dest_dir = tmp_path / "out"
    dest_dir.mkdir()
    with mock.patch.object(R3Resource, "_ensure_cached", return_value=cached):
        # ensure_cached() returns a Path; its parent is a Path too.
        assert isinstance(res.ensure_cached(), Path)
        with mock.patch("recount3.resource._hardlink_or_copy"):
            result = res.download(dest_dir, cache_mode="enable")
    assert result == str(dest_dir / "homes_index")


def test_filepath_pathlike_is_normalized_to_str(cfg: Config) -> None:
    """A Path passed to the constructor is stored as a str."""
    desc = R3ResourceDescription(resource_type="data_sources", organism="human")
    res = R3Resource(description=desc, config=cfg, filepath=Path("/a/b.gz"))
    assert isinstance(res.filepath, str)
    assert res.filepath == str(Path("/a/b.gz"))
    assert "PosixPath(" not in repr(res)
    assert "WindowsPath(" not in repr(res)


def test_filepath_pathlike_and_str_compare_equal(cfg: Config) -> None:
    """Two resources naming the same file are equal regardless of spelling."""
    desc = R3ResourceDescription(resource_type="data_sources", organism="human")
    as_path = R3Resource(description=desc, config=cfg, filepath=Path("/a/b.gz"))
    as_str = R3Resource(
        description=desc, config=cfg, filepath=str(Path("/a/b.gz"))
    )
    assert as_path == as_str


# ===========================================================================
# R3Resource.load — in-memory cache
# ===========================================================================


def test_load_returns_cached_without_force(cfg: Config) -> None:
    """load() returns _cached_data without disk I/O when force=False."""
    res = _make("data_sources", cfg, organism="human")
    sentinel = object()
    res._cached_data = sentinel
    assert res.load() is sentinel


def test_load_force_bypasses_in_memory_cache(cfg: Config) -> None:
    """load(force=True) re-reads from disk ignoring _cached_data."""
    res = _make(
        "metadata_files",
        cfg,
        organism="human",
        data_source="sra",
        project="SRP014565",
        table_name="recount_qc",
    )
    res._cached_data = "stale"
    _seed(res, _META_MD_GZ)
    result = res.load(force=True)
    assert isinstance(result, pd.DataFrame)
    assert result is not "stale"  # noqa: F632


# ===========================================================================
# R3Resource.load — bigwig_files
# ===========================================================================


@pytest.mark.requires_pybigwig
def test_load_bigwig_returns_open_bigwigfile(cfg: Config) -> None:
    """load() for bigwig_files returns an open BigWigFile instance."""
    res = _make(
        "bigwig_files",
        cfg,
        organism="human",
        data_source="sra",
        project="SRP009615",
        sample="SRR387777",
    )
    _seed(res, _BW_PATH)
    result = res.load()
    assert isinstance(result, BigWigFile)
    assert result.is_open()
    res.clear_loaded()


@pytest.mark.requires_pybigwig
def test_load_bigwig_second_call_returns_cached(cfg: Config) -> None:
    """Second load() without force=True returns the same BigWigFile."""
    res = _make(
        "bigwig_files",
        cfg,
        organism="human",
        data_source="sra",
        project="SRP009615",
        sample="SRR387777",
    )
    _seed(res, _BW_PATH)
    first = res.load()
    second = res.load()
    assert first is second
    res.clear_loaded()


def test_load_bigwig_file_not_found_raises(cfg: Config) -> None:
    """FileNotFoundError when the bigwig cache file is absent post-download."""
    res = _make(
        "bigwig_files",
        cfg,
        organism="human",
        data_source="sra",
        project="SRP009615",
        sample="SRR387777",
    )
    with mock.patch.object(R3Resource, "download"):
        with pytest.raises(FileNotFoundError):
            res.load()


# ===========================================================================
# R3Resource.load — count_files_gene_or_exon
# ===========================================================================


def _make_gene_res(cfg: Config) -> R3Resource:
    return _make(
        "count_files_gene_or_exon",
        cfg,
        organism="human",
        data_source="sra",
        genomic_unit="gene",
        project="SRP014565",
        annotation_extension="G026",
    )


def test_load_gene_counts_returns_dataframe_indexed_by_gene_id(
    cfg: Config,
) -> None:
    """load() for gene counts returns a DataFrame with gene_id index."""
    res = _make_gene_res(cfg)
    _seed(res, _GENE_GZ)
    df = res.load()
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "gene_id"
    assert df.shape[1] >= 1


def test_load_gene_counts_caches_result(cfg: Config) -> None:
    """Second load() on gene counts returns the same DataFrame object."""
    res = _make_gene_res(cfg)
    _seed(res, _GENE_GZ)
    assert res.load() is res.load()


def test_load_gene_counts_file_not_found_raises(cfg: Config) -> None:
    """FileNotFoundError when gene/exon counts file is absent post-download."""
    res = _make_gene_res(cfg)
    with mock.patch.object(R3Resource, "download"):
        with pytest.raises(FileNotFoundError):
            res.load()


def test_load_gene_counts_python_engine_fallback(
    cfg: Config, tmp_path: Path
) -> None:
    """Python engine is used as fallback when C engine raises."""
    p = tmp_path / "fallback.tsv.gz"
    df_src = pd.DataFrame({"gene_id": ["g1", "g2"], "SRR001": [1, 2]})
    with gzip.open(p, "wt") as fh:
        df_src.to_csv(fh, sep="\t", index=False)

    original = pd.read_csv
    call_count = 0

    def _patched(*args: object, **kwargs: object) -> pd.DataFrame:
        nonlocal call_count
        call_count += 1
        if call_count == 1:
            raise pd.errors.ParserError("simulated C fail")
        return original(*args, **kwargs)  # type: ignore[arg-type]

    desc = R3ResourceDescription(
        resource_type="count_files_gene_or_exon",
        organism="human",
        data_source="sra",
        genomic_unit="gene",
        project="SRP014565",
        annotation_extension="G026",
    )
    res = R3Resource(description=desc, config=cfg)
    _seed(res, p)
    with mock.patch("recount3.resource.pd.read_csv", _patched):
        result = res.load()
    assert isinstance(result, pd.DataFrame)
    assert call_count == 2


def test_load_exon_counts_uses_exon_id_index(
    cfg: Config, tmp_path: Path
) -> None:
    """exon_id is used as the index column when gene_id is absent."""
    p = tmp_path / "exon.tsv.gz"
    df_src = pd.DataFrame({"exon_id": ["e1", "e2"], "SRR001": [10, 20]})
    with gzip.open(p, "wt") as fh:
        df_src.to_csv(fh, sep="\t", index=False)
    desc = R3ResourceDescription(
        resource_type="count_files_gene_or_exon",
        organism="human",
        data_source="sra",
        genomic_unit="exon",
        project="SRP014565",
        annotation_extension="G026",
    )
    res = R3Resource(description=desc, config=cfg)
    _seed(res, p)
    result = res.load()
    assert result.index.name == "exon_id"  # type: ignore


def test_load_counts_uses_feature_id_index(cfg: Config, tmp_path: Path) -> None:
    """feature_id is used as index when gene_id and exon_id are absent."""
    p = tmp_path / "feat.tsv.gz"
    df_src = pd.DataFrame({"feature_id": ["f1", "f2"], "SRR001": [5, 15]})
    with gzip.open(p, "wt") as fh:
        df_src.to_csv(fh, sep="\t", index=False)
    desc = R3ResourceDescription(
        resource_type="count_files_gene_or_exon",
        organism="human",
        data_source="sra",
        genomic_unit="gene",
        project="SRP014565",
        annotation_extension="G026",
    )
    res = R3Resource(description=desc, config=cfg)
    _seed(res, p)
    result = res.load()
    assert result.index.name == "feature_id"  # type: ignore


def test_load_counts_uses_first_column_as_fallback_index(
    cfg: Config, tmp_path: Path
) -> None:
    """First column is used as index when no known name is found."""
    p = tmp_path / "other.tsv.gz"
    df_src = pd.DataFrame({"row_name": ["r1", "r2"], "SRR001": [1, 2]})
    with gzip.open(p, "wt") as fh:
        df_src.to_csv(fh, sep="\t", index=False)
    desc = R3ResourceDescription(
        resource_type="count_files_gene_or_exon",
        organism="human",
        data_source="sra",
        genomic_unit="gene",
        project="SRP014565",
        annotation_extension="G026",
    )
    res = R3Resource(description=desc, config=cfg)
    _seed(res, p)
    result = res.load()
    assert result.index.name == "row_name"  # type: ignore


def test_load_counts_empty_matrix_raises_load_error(
    cfg: Config, tmp_path: Path
) -> None:
    """LoadError when parsed counts file is completely empty (no rows)."""
    p = tmp_path / "empty.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("gene_id\n")  # header only, zero data rows
    desc = R3ResourceDescription(
        resource_type="count_files_gene_or_exon",
        organism="human",
        data_source="sra",
        genomic_unit="gene",
        project="SRP014565",
        annotation_extension="G026",
    )
    res = R3Resource(description=desc, config=cfg)
    _seed(res, p)
    with pytest.raises(LoadError, match="empty or 1-column"):
        res.load()


def test_load_counts_single_column_raises_load_error(
    cfg: Config, tmp_path: Path
) -> None:
    """LoadError when parsed counts file has exactly one column."""
    p = tmp_path / "onecol.tsv.gz"
    df_src = pd.DataFrame({"gene_id": ["g1", "g2"]})
    with gzip.open(p, "wt") as fh:
        df_src.to_csv(fh, sep="\t", index=False)
    desc = R3ResourceDescription(
        resource_type="count_files_gene_or_exon",
        organism="human",
        data_source="sra",
        genomic_unit="gene",
        project="SRP014565",
        annotation_extension="G026",
    )
    res = R3Resource(description=desc, config=cfg)
    _seed(res, p)
    with pytest.raises(LoadError, match="empty or 1-column"):
        res.load()


# ===========================================================================
# R3Resource.load — count_files_junctions (MM)
# ===========================================================================


def _make_jxn(cfg: Config, ext: str) -> R3Resource:
    return _make(
        "count_files_junctions",
        cfg,
        organism="human",
        data_source="sra",
        project="SRP014565",
        junction_type="ALL",
        junction_extension=ext,
    )


def _seed_jxn_mm_and_id(res: R3Resource, cfg: Config) -> None:
    """Pre-seed the MM cache file and its ID sidecar."""
    _seed(res, _JXN_MM_GZ)
    id_url = _derive_junction_sidecar_url(res.url or "", "ID")
    id_cache = _cache_path(id_url, cfg.cache_dir)
    id_cache.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(_JXN_ID_GZ, id_cache)


def test_load_junctions_mm_returns_sparse_dataframe(cfg: Config) -> None:
    """load() for MM junctions returns a sparse DataFrame."""
    res = _make_jxn(cfg, "MM")
    _seed_jxn_mm_and_id(res, cfg)
    df = res.load()
    assert isinstance(df, pd.DataFrame)
    assert df.shape[1] == 1  # one sample in the test ID file


def test_load_junctions_mm_file_not_found_raises(cfg: Config) -> None:
    """FileNotFoundError when the MM file is absent after download."""
    res = _make_jxn(cfg, "MM")
    with mock.patch.object(R3Resource, "download"):
        with pytest.raises(FileNotFoundError):
            res.load()


def test_load_junctions_mm_bad_sidecar_url_raises_load_error(
    cfg: Config,
) -> None:
    """LoadError when _derive_junction_sidecar_url raises."""
    res = _make_jxn(cfg, "MM")
    _seed(res, _JXN_MM_GZ)
    with mock.patch(
        "recount3.resource._derive_junction_sidecar_url",
        side_effect=ValueError("bad url"),
    ):
        with pytest.raises(LoadError, match="Cannot derive junction ID URL"):
            res.load()


def test_load_junctions_mm_column_mismatch_raises_load_error(
    cfg: Config, tmp_path: Path
) -> None:
    """LoadError when MM column count does not match ID rail_id count."""
    res = _make_jxn(cfg, "MM")
    _seed(res, _JXN_MM_GZ)
    # Create an ID file with 2 rail_ids, but the MM file has 1 column.
    p = tmp_path / "id2.tsv.gz"
    df_id = pd.DataFrame({"rail_id": ["111", "222"]})
    with gzip.open(p, "wt") as fh:
        df_id.to_csv(fh, sep="\t", index=False)
    id_url = _derive_junction_sidecar_url(res.url or "", "ID")
    id_cache = _cache_path(id_url, cfg.cache_dir)
    id_cache.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(p, id_cache)
    with pytest.raises(LoadError, match="column count does not match"):
        res.load()


def test_load_junctions_mm_undefined_shape_raises_load_error(
    cfg: Config,
) -> None:
    """LoadError when the loaded MM matrix reports shape=None."""
    res = _make_jxn(cfg, "MM")
    _seed_jxn_mm_and_id(res, cfg)

    class _ShapeNoneMat:
        shape = None
        ndim = 2

    with mock.patch(
        "recount3.resource._read_mm_matrix",
        return_value=_ShapeNoneMat(),
    ):
        with pytest.raises(LoadError, match="undefined shape"):
            res.load()


def test_load_junctions_mm_config_none_calls_default_config(
    tmp_path: Path,
) -> None:
    """load() for MM junctions calls default_config() when config=None."""
    fc = _fake_cfg(tmp_path)
    desc = R3ResourceDescription(
        resource_type="count_files_junctions",
        organism="human",
        data_source="sra",
        project="SRP014565",
        junction_type="ALL",
        junction_extension="MM",
    )
    with mock.patch("recount3.resource.default_config", return_value=fc):
        res = R3Resource(description=desc)
        _seed(res, _JXN_MM_GZ)
        id_url = _derive_junction_sidecar_url(res.url or "", "ID")
        id_cache = _cache_path(id_url, fc.cache_dir)
        id_cache.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(_JXN_ID_GZ, id_cache)
        df = res.load()
    assert isinstance(df, pd.DataFrame)


# ===========================================================================
# R3Resource.load — count_files_junctions (ID / RR / unknown)
# ===========================================================================


def test_load_junctions_id_returns_dataframe(cfg: Config) -> None:
    """load() for ID junctions reads the table directly."""
    res = _make_jxn(cfg, "ID")
    _seed(res, _JXN_ID_GZ)
    result = res.load()
    assert isinstance(result, pd.DataFrame)


def test_load_junctions_rr_returns_dataframe(cfg: Config) -> None:
    """load() for RR junctions reads the table directly."""
    res = _make_jxn(cfg, "RR")
    _seed(res, _JXN_RR_GZ)
    result = res.load()
    assert isinstance(result, pd.DataFrame)


def test_load_junctions_unsupported_extension_raises(
    cfg: Config, tmp_path: Path
) -> None:
    """LoadError for an unrecognised junction_extension string."""
    p = tmp_path / "data.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("col\nval\n")
    res = _make_jxn(cfg, "XY")
    _seed(res, p)
    with pytest.raises(LoadError, match="Unsupported junction_extension"):
        res.load()


def test_load_junctions_file_not_found_raises(cfg: Config) -> None:
    """FileNotFoundError when junction cache file is absent post-download."""
    res = _make_jxn(cfg, "MM")
    with mock.patch.object(R3Resource, "download"):
        with pytest.raises(FileNotFoundError):
            res.load()


# ===========================================================================
# R3Resource.load — generic fallback (metadata, tsv.gz, tsv, unknown)
# ===========================================================================


def test_load_generic_md_gz_returns_dataframe(cfg: Config) -> None:
    """load() for metadata (.MD.gz) reads via the generic read_table path."""
    res = _make(
        "metadata_files",
        cfg,
        organism="human",
        data_source="sra",
        project="SRP014565",
        table_name="recount_qc",
    )
    _seed(res, _META_MD_GZ)
    result = res.load()
    assert isinstance(result, pd.DataFrame)
    assert not result.empty


def test_load_generic_tsv_gz_returns_dataframe(
    cfg: Config, tmp_path: Path
) -> None:
    """load() for a .tsv.gz URL uses read_table via the generic path."""
    p = tmp_path / "data.tsv.gz"
    with gzip.open(p, "wt") as fh:
        pd.DataFrame({"a": [1, 2], "b": [3, 4]}).to_csv(
            fh, sep="\t", index=False
        )
    desc = R3ResourceDescription(resource_type="data_sources", organism="human")
    res = R3Resource(
        description=desc,
        url="https://example.org/recount3/human/thing.tsv.gz",
        config=cfg,
    )
    _seed(res, p)
    result = res.load()
    assert isinstance(result, pd.DataFrame)


def test_load_generic_tsv_returns_dataframe(
    cfg: Config, tmp_path: Path
) -> None:
    """load() for a plain .tsv URL uses read_table via the generic path."""
    p = tmp_path / "data.tsv"
    pd.DataFrame({"a": [1, 2]}).to_csv(p, sep="\t", index=False)
    desc = R3ResourceDescription(resource_type="data_sources", organism="human")
    res = R3Resource(
        description=desc,
        url="https://example.org/recount3/human/thing.tsv",
        config=cfg,
    )
    _seed(res, p)
    result = res.load()
    assert isinstance(result, pd.DataFrame)


def test_load_generic_unsupported_extension_raises(
    cfg: Config, tmp_path: Path
) -> None:
    """LoadError for a cached file whose name has an unsupported extension."""
    p = tmp_path / "data.xyz"
    p.write_bytes(b"dummy")
    desc = R3ResourceDescription(resource_type="data_sources", organism="human")
    res = R3Resource(
        description=desc,
        url="https://example.org/recount3/human/thing.xyz",
        config=cfg,
    )
    _seed(res, p)
    with pytest.raises(LoadError, match="Unsupported load"):
        res.load()


def test_load_generic_file_not_found_raises(cfg: Config) -> None:
    """FileNotFoundError when generic cached file is absent post-download."""
    res = _make(
        "metadata_files",
        cfg,
        organism="human",
        data_source="sra",
        project="SRP014565",
        table_name="recount_qc",
    )
    with mock.patch.object(R3Resource, "download"):
        with pytest.raises(FileNotFoundError):
            res.load()


# ===========================================================================
# R3Resource.is_loaded / get_loaded / clear_loaded
# ===========================================================================


def test_is_loaded_false_initially(cfg: Config) -> None:
    """is_loaded() returns False before any successful load call."""
    assert _make("data_sources", cfg, organism="human").is_loaded() is False


def test_is_loaded_true_after_load(cfg: Config) -> None:
    """is_loaded() returns True after a successful load call."""
    res = _make(
        "metadata_files",
        cfg,
        organism="human",
        data_source="sra",
        project="SRP014565",
        table_name="recount_qc",
    )
    _seed(res, _META_MD_GZ)
    res.load()
    assert res.is_loaded() is True


def test_get_loaded_returns_none_before_load(cfg: Config) -> None:
    """get_loaded() returns None when nothing has been loaded."""
    assert _make("data_sources", cfg, organism="human").get_loaded() is None


def test_get_loaded_returns_object_after_load(cfg: Config) -> None:
    """get_loaded() returns the loaded object without triggering I/O."""
    res = _make(
        "metadata_files",
        cfg,
        organism="human",
        data_source="sra",
        project="SRP014565",
        table_name="recount_qc",
    )
    _seed(res, _META_MD_GZ)
    loaded = res.load()
    assert res.get_loaded() is loaded


def test_clear_loaded_noop_when_nothing_cached(cfg: Config) -> None:
    """clear_loaded() is a safe no-op when nothing is in memory."""
    res = _make("data_sources", cfg, organism="human")
    res.clear_loaded()
    assert not res.is_loaded()


def test_clear_loaded_evicts_dataframe(cfg: Config) -> None:
    """clear_loaded() removes a cached DataFrame from memory."""
    res = _make(
        "metadata_files",
        cfg,
        organism="human",
        data_source="sra",
        project="SRP014565",
        table_name="recount_qc",
    )
    _seed(res, _META_MD_GZ)
    res.load()
    assert res.is_loaded()
    res.clear_loaded()
    assert not res.is_loaded()


@pytest.mark.requires_pybigwig
def test_clear_loaded_closes_bigwig(cfg: Config) -> None:
    """clear_loaded() calls close() on a cached BigWigFile."""
    res = _make(
        "bigwig_files",
        cfg,
        organism="human",
        data_source="sra",
        project="SRP009615",
        sample="SRR387777",
    )
    _seed(res, _BW_PATH)
    res.load()
    bw = res.get_loaded()
    assert isinstance(bw, BigWigFile)
    res.clear_loaded()
    assert not res.is_loaded()
    assert not bw.is_open()


def test_clear_loaded_sets_none_even_if_close_raises(cfg: Config) -> None:
    """_cached_data is cleared in the finally block even if close() raises."""
    mock_bw = mock.MagicMock(spec=BigWigFile)
    mock_bw.close.side_effect = OSError("handle gone")
    res = _make("data_sources", cfg, organism="human")
    res._cached_data = mock_bw
    with pytest.raises(OSError, match="handle gone"):
        res.clear_loaded()
    assert not res.is_loaded()


# ===========================================================================
# R3Resource.__repr__
# ===========================================================================


def test_repr_contains_expected_fields(cfg: Config) -> None:
    """__repr__ includes class name, url=, arcname=, and filepath=."""
    res = _make("data_sources", cfg, organism="human")
    r = repr(res)
    assert r.startswith("R3Resource(")
    assert "url=" in r
    assert "arcname=" in r
    assert "filepath=" in r
    assert "human/homes_index" in r


def test_repr_reflects_filepath_after_set(cfg: Config) -> None:
    """__repr__ shows an updated filepath once download() sets it."""
    res = _make("data_sources", cfg, organism="human")
    res.filepath = "/some/path/homes_index"
    assert "/some/path/homes_index" in repr(res)


# ===========================================================================
# build_url
# ===========================================================================


def test_build_url_with_explicit_config(cfg: Config) -> None:
    """Joins the description's url_path onto the supplied config's base_url."""
    url = res_module.build_url(
        "bigwig_files",
        config=cfg,
        organism="human",
        data_source="sra",
        project="SRP009615",
        sample="SRR387777",
    )
    expected = R3ResourceDescription(
        resource_type="bigwig_files",
        organism="human",
        data_source="sra",
        project="SRP009615",
        sample="SRR387777",
    ).url_path()
    assert url == cfg.base_url + expected
    assert url.endswith("sra.base_sums.SRP009615_SRR387777.ALL.bw")


def test_build_url_uses_default_config_when_omitted(tmp_path: Path) -> None:
    """Falls back to default_config() when no config is given."""
    fallback = _fake_cfg(tmp_path)
    with mock.patch("recount3.resource.default_config", return_value=fallback):
        url = res_module.build_url(
            "bigwig_files",
            organism="human",
            data_source="sra",
            project="SRP009615",
            sample="SRR387777",
        )
    assert url.startswith("https://fallback.org/recount3/")


def test_build_url_invalid_resource_type_raises(cfg: Config) -> None:
    """Propagates the description factory's error for an unknown type."""
    with pytest.raises((KeyError, ValueError)):
        res_module.build_url("not_a_real_type", config=cfg, organism="human")


# ===========================================================================
# R3Resource.from_mapping
# ===========================================================================


def test_from_mapping_builds_configured_resource(cfg: Config) -> None:
    """Rehydrates a resource with the right description, config, and url."""
    res = R3Resource.from_mapping(
        {"resource_type": "data_sources", "organism": "human"},
        config=cfg,
    )
    assert isinstance(res, R3Resource)
    assert res.config is cfg
    assert res.description.organism == "human"
    assert res.url == cfg.base_url + res.description.url_path()


def test_from_mapping_ignores_derived_url_and_arcname_keys(cfg: Config) -> None:
    """Stale serialized url/arcname keys are dropped; url is recomputed."""
    res = R3Resource.from_mapping(
        {
            "resource_type": "data_sources",
            "organism": "human",
            "url": "STALE",
            "arcname": "STALE",
        },
        config=cfg,
    )
    assert res.url != "STALE"
    assert res.url == cfg.base_url + res.description.url_path()


def test_from_mapping_defaults_config_when_omitted(tmp_path: Path) -> None:
    """With no config, the url derives from default_config()'s base_url."""
    fallback = _fake_cfg(tmp_path)
    with mock.patch("recount3.resource.default_config", return_value=fallback):
        res = R3Resource.from_mapping(
            {"resource_type": "data_sources", "organism": "human"}
        )
    assert res.config is None
    assert res.url.startswith("https://fallback.org/recount3/")


def test_read_metadata_table_maps_r_logicals_to_boolean(
    tmp_path: Path,
) -> None:
    """R writes logicals as TRUE/FALSE/T/F; they must round-trip as booleans.

    Anything else stays as written: a text column with blanks keeps them, and
    ``NA`` becomes a missing value rather than the string.
    """
    path = tmp_path / "metadata.tsv"
    path.write_text(
        "sample\tpaired\tflag\tnote\n"
        "S1\tTRUE\tT\tok\n"
        "S2\tFALSE\tF\t\n"
        "S3\tNA\tTRUE\tx\n",
        encoding="utf-8",
    )

    frame = res_module._read_metadata_table(path)

    assert frame["paired"].dtype == "boolean"
    assert list(frame["paired"]) == [True, False, pd.NA]
    assert frame["flag"].dtype == "boolean"
    assert list(frame["flag"]) == [True, False, True]
    assert list(frame["note"]) == ["ok", "", "x"]
    assert list(frame["sample"]) == ["S1", "S2", "S3"]


def _fetch(
    cfg: Config, url: str = _CACHE_URL, *, refresh: bool = False
) -> Path:
    """Drive :func:`_ensure_cached_url` with a chunk size small enough to loop."""
    return _ensure_cached_url(
        url=url,
        cache_root=cfg.cache_dir,
        cfg=cfg,
        chunk_size=32,
        refresh=refresh,
    )


def _counts_resource(cfg: Config, project: str = "SRP001") -> R3Resource:
    """Return a gene-counts resource wired to the temporary cache."""
    return R3Resource(
        R3ResourceDescription(
            resource_type="count_files_gene_or_exon",
            organism="human",
            data_source="sra",
            project=project,
            genomic_unit="gene",
            annotation_extension="G026",
        ),
        config=cfg,
    )


@pytest.mark.parametrize("refresh", [False, True])
@pytest.mark.parametrize("backend", ["filesystem", "pybiocfilecache"])
def test_ensure_cached_url_lets_distinct_destinations_overlap(
    cfg: Config, refresh: bool, backend: str
) -> None:
    """Every transfer enters before any can finish, with either backend.

    The destination lock must not serialize unrelated URLs, so a barrier that
    only trips once all four transfers are in flight has to be reachable.
    """
    if backend == "pybiocfilecache":
        pytest.importorskip("pybiocfilecache")
    cfg = dataclasses.replace(cfg, cache_backend=backend)
    barrier = threading.Barrier(4, timeout=5)

    def transfer(url: str, path: Path, **_: object) -> None:
        barrier.wait()
        path.write_bytes(url.encode())

    with mock.patch.object(res_module, "download_to_file", transfer):
        with futures.ThreadPoolExecutor(4) as pool:
            pending = [
                pool.submit(
                    _fetch,
                    cfg,
                    f"https://example.org/{index}",
                    refresh=refresh,
                )
                for index in range(4)
            ]
            paths = [future.result(timeout=10) for future in pending]

    assert len(set(paths)) == 4
    assert all(path.read_bytes() for path in paths)


@pytest.mark.parametrize("refresh", [False, True])
@pytest.mark.parametrize("fail", [False, True])
def test_ensure_cached_url_same_destination_waits_and_recovers(
    cfg: Config, refresh: bool, fail: bool
) -> None:
    """A confirmed waiter deduplicates misses or performs its own refresh.

    The second caller is observed blocking on the destination lock while the
    first transfer is still running. It then downloads only when it has to: a
    forced refresh always transfers, and a failed first transfer leaves nothing
    to reuse.
    """
    entered = threading.Event()
    waiting = threading.Event()
    release = threading.Event()
    calls: list[str] = []
    real_lock = _path_lock_for_path

    @contextlib.contextmanager
    def observed_lock(path: Path):
        if entered.is_set():
            waiting.set()
        with real_lock(path):
            yield

    def transfer(url: str, path: Path, **_: object) -> None:
        calls.append(url)
        if len(calls) == 1:
            entered.set()
            assert release.wait(5)
            if fail:
                raise DownloadError("controlled failure")
        path.write_bytes(str(len(calls)).encode())

    with mock.patch.object(res_module, "_path_lock_for_path", observed_lock):
        with mock.patch.object(res_module, "download_to_file", transfer):
            with futures.ThreadPoolExecutor(2) as pool:
                first = pool.submit(_fetch, cfg, refresh=refresh)
                try:
                    assert entered.wait(5)
                    second = pool.submit(_fetch, cfg, refresh=refresh)
                    assert waiting.wait(5)
                    assert len(calls) == 1
                finally:
                    release.set()
                if fail:
                    with pytest.raises(DownloadError):
                        first.result(timeout=5)
                else:
                    first.result(timeout=5)
                result = second.result(timeout=5)

    assert len(calls) == (2 if refresh or fail else 1)
    assert result.read_bytes() == str(len(calls)).encode()


def test_ensure_cached_url_failed_refresh_keeps_the_old_payload(
    cfg: Config,
) -> None:
    """A broken stream preserves the old file and cleans its temporary file."""
    path = _cache_path(_CACHE_URL, cfg.cache_dir)
    path.write_bytes(b"old")

    class BrokenStream(io.BytesIO):
        """Fails once the first chunk has been handed over."""

        def read(self, *args: object) -> bytes:
            if self.tell():
                raise OSError("broken transfer")
            return super().read(*args)

    with mock.patch(
        "recount3._utils.http_open", return_value=BrokenStream(b"new")
    ):
        with pytest.raises(DownloadError):
            _fetch(cfg, refresh=True)

    assert path.read_bytes() == b"old"
    assert list(cfg.cache_dir.iterdir()) == [path]

    with mock.patch(
        "recount3._utils.http_open", return_value=io.BytesIO(b"new")
    ):
        _fetch(cfg, refresh=True)
    assert path.read_bytes() == b"new"


def test_ensure_cached_url_adopts_a_native_payload_into_the_registry(
    cfg: Config,
) -> None:
    """An existing file is registered without a transfer and then re-used."""
    module = pytest.importorskip("pybiocfilecache")
    cfg = dataclasses.replace(cfg, cache_backend="pybiocfilecache")
    path = _cache_path(_CACHE_URL, cfg.cache_dir)
    path.write_bytes(b"native")

    with mock.patch.object(res_module, "download_to_file") as download:
        assert _fetch(cfg) == path
        assert _fetch(cfg) == path
        download.assert_not_called()

    with module.BiocFileCache(cfg.cache_dir) as registry:
        record = registry.get(rname=_CACHE_URL)
    assert Path(record["rpath"]) == path
    assert record["etag"]


def test_ensure_cached_url_refresh_rechecksums_the_registry_entry(
    cfg: Config,
) -> None:
    """A refreshed payload replaces the recorded checksum, not the row."""
    module = pytest.importorskip("pybiocfilecache")
    cfg = dataclasses.replace(cfg, cache_backend="pybiocfilecache")
    with mock.patch(
        "recount3._utils.http_open", return_value=io.BytesIO(b"old")
    ):
        path = _fetch(cfg)
    with module.BiocFileCache(cfg.cache_dir) as registry:
        checksum = registry.get(rname=_CACHE_URL)["etag"]

    with mock.patch(
        "recount3._utils.http_open", return_value=io.BytesIO(b"new")
    ):
        _fetch(cfg, refresh=True)

    assert path.read_bytes() == b"new"
    with module.BiocFileCache(cfg.cache_dir) as registry:
        assert registry.get(rname=_CACHE_URL)["etag"] != checksum
        assert len(registry) == 1


def test_ensure_cached_url_recovers_after_the_cache_is_emptied(
    cfg: Config,
) -> None:
    """Clearing the cache leaves the registry usable for the next transfer."""
    pytest.importorskip("pybiocfilecache")
    cfg = dataclasses.replace(cfg, cache_backend="pybiocfilecache")
    with mock.patch(
        "recount3._utils.http_open", return_value=io.BytesIO(b"first")
    ):
        assert _fetch(cfg).read_bytes() == b"first"

    assert recount3_cache_rm(config=cfg)
    assert recount3_cache_rm(config=cfg) == []
    assert (cfg.cache_dir / "BiocFileCache.sqlite").exists()

    with mock.patch(
        "recount3._utils.http_open", return_value=io.BytesIO(b"again")
    ):
        assert _fetch(cfg).read_bytes() == b"again"


def test_ensure_cached_url_registry_failure_keeps_the_payload(
    cfg: Config,
) -> None:
    """A registry error never destroys a completed download."""
    pytest.importorskip("pybiocfilecache")
    cfg = dataclasses.replace(cfg, cache_backend="pybiocfilecache")
    with mock.patch(
        "recount3._utils.http_open", return_value=io.BytesIO(b"good")
    ):
        with mock.patch.object(
            res_module, "_register_biocfilecache", side_effect=OSError
        ):
            with pytest.raises(OSError):
                _fetch(cfg)

    with mock.patch.object(res_module, "download_to_file") as download:
        assert _fetch(cfg).read_bytes() == b"good"
        download.assert_not_called()


def test_ensure_cached_url_repairs_a_failed_refresh_on_the_next_hit(
    cfg: Config,
) -> None:
    """A completed refresh with a failed DB update recovers without network."""
    module = pytest.importorskip("pybiocfilecache")
    cfg = dataclasses.replace(cfg, cache_backend="pybiocfilecache")
    with mock.patch(
        "recount3._utils.http_open", return_value=io.BytesIO(b"old")
    ):
        _fetch(cfg)
    with mock.patch.object(
        res_module, "_register_biocfilecache", side_effect=OSError
    ):
        with mock.patch(
            "recount3._utils.http_open", return_value=io.BytesIO(b"new")
        ):
            with pytest.raises(OSError):
                _fetch(cfg, refresh=True)

    with mock.patch.object(res_module, "download_to_file") as download:
        assert _fetch(cfg).read_bytes() == b"new"
        download.assert_not_called()

    with module.BiocFileCache(cfg.cache_dir) as registry:
        assert registry.get(rname=_CACHE_URL)["etag"] == (
            hashlib.md5(b"new").hexdigest()
        )


def test_resource_reuses_and_repairs_an_r_created_relative_entry(
    cfg: Config, tmp_path: Path
) -> None:
    """R-style relative entries are reused and atomically repaired in place."""
    module = pytest.importorskip("pybiocfilecache")
    cfg = dataclasses.replace(cfg, cache_backend="pybiocfilecache")
    res = _counts_resource(cfg)
    source = tmp_path / "from-r.tsv.gz"
    source.write_bytes(b"R bytes")
    with module.BiocFileCache(cfg.cache_dir) as cache:
        row = cache.add(rname=res.url, fpath=source, rtype="relative")
        registered = Path(row["rpath"])

    assert res._cached_path() == registered
    with mock.patch.object(res_module, "download_to_file") as download:
        assert res.ensure_cached(download=False) == registered
        res.download()
        download.assert_not_called()

    registered.unlink()
    assert registered in recount3_cache_files(cfg)
    with mock.patch(
        "recount3._utils.http_open", return_value=io.BytesIO(b"fixed")
    ):
        assert res.ensure_cached().read_bytes() == b"fixed"

    with module.BiocFileCache(cfg.cache_dir) as cache:
        assert cache.get(rname=res.url)["rtype"] == "relative"


def test_download_to_a_new_directory_creates_the_parent(
    cfg: Config, tmp_path: Path
) -> None:
    """Cached directory downloads work with a new destination hierarchy."""
    res = _counts_resource(cfg)
    res._cached_path().write_bytes(b"cached")

    result = res.download(str(tmp_path / "new" / "nested"))

    assert Path(result).read_bytes() == b"cached"


def test_optional_backend_reports_the_extra_and_disable_bypasses_it(
    cfg: Config, tmp_path: Path
) -> None:
    """The extra is requested lazily and is unnecessary for direct downloads."""
    cfg = dataclasses.replace(cfg, cache_backend="pybiocfilecache")
    path = _cache_path(_CACHE_URL, cfg.cache_dir)
    path.write_bytes(b"cached")
    import_optional_module.cache_clear()
    try:
        with mock.patch(
            "recount3._utils.importlib.import_module",
            side_effect=ModuleNotFoundError("pybiocfilecache"),
        ):
            # Even a warm cache needs the registry to resolve the destination.
            with pytest.raises(
                ImportError, match=r"recount3\[pybiocfilecache\]"
            ):
                _fetch(cfg)
            path.unlink()
            with mock.patch.object(res_module, "download_to_file") as download:
                with pytest.raises(
                    ImportError, match=r"recount3\[pybiocfilecache\]"
                ):
                    _fetch(cfg)
                download.assert_not_called()

            res = _counts_resource(cfg)
            with mock.patch.object(res_module, "download_to_file") as download:
                res.download(str(tmp_path / "direct"), cache_mode="disable")
                download.assert_called_once()
    finally:
        import_optional_module.cache_clear()
