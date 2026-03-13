# pylint: disable=redefined-outer-name

from __future__ import annotations

import gzip
import shutil
from pathlib import Path
from unittest import mock

import pandas as pd
import pytest
import scipy.io
import scipy.sparse

import recount3.resource as res_module
from recount3._bigwig import BigWigFile
from recount3._descriptions import R3ResourceDescription
from recount3._utils import _cache_path, _derive_junction_sidecar_url
from recount3.config import Config
from recount3.errors import LoadError
from recount3.resource import (
    R3Resource,
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
    dest = resource._cached_path()  # pylint: disable=protected-access
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
        def __enter__(self) -> _SeedOnEnter:
            cache_p.touch()
            return self

        def __exit__(self, *args: object) -> bool:
            return False

    with mock.patch.object(res_module, "_FILE_LOCK", _SeedOnEnter()):
        with mock.patch(
            "recount3.resource.download_to_file"
        ) as mock_dl:
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
    with mock.patch(
        "recount3.resource.download_to_file"
    ) as mock_dl:
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

    This covers the guard at line 115 which is unreachable via normal I/O
    (a non-empty DataFrame always yields ≥1 items). It is reached by mocking
    pd.read_csv to return a DataFrame that passes df.empty but whose column
    series returns [] from tolist().
    """
    p = tmp_path / "data.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("rail_id\nsome_value\n")

    mock_series = mock.MagicMock()
    mock_series.astype.return_value.tolist.return_value = []

    mock_df = mock.MagicMock()
    mock_df.empty = False
    mock_df.columns = ["rail_id"]
    mock_df.__getitem__ = mock.Mock(return_value=mock_series)

    with mock.patch(
        "recount3.resource.pd.read_csv", return_value=mock_df
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
    """Reads a .MM.gz file into a CSR sparse matrix."""
    mat = _read_mm_matrix(_JXN_MM_GZ)
    assert isinstance(mat, scipy.sparse.csr_matrix)
    assert mat.ndim == 2


def test_read_mm_matrix_plain_file(tmp_path: Path) -> None:
    """Reads a plain (non-compressed) MatrixMarket file."""
    # scipy.io.mmwrite appends .mtx when the name has no .mtx suffix;
    # use .mtx directly so the written file matches the path we read.
    mm_file = tmp_path / "test.mtx"
    scipy.io.mmwrite(str(mm_file), scipy.sparse.eye(4, k=1, format="coo"))
    mat = _read_mm_matrix(mm_file)
    assert isinstance(mat, scipy.sparse.csr_matrix)
    assert mat.shape == (4, 4)


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
    desc = R3ResourceDescription(
        resource_type="data_sources", organism="human"
    )
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
    desc = R3ResourceDescription(
        resource_type="data_sources", organism="human"
    )
    fc = _fake_cfg(tmp_path)
    with mock.patch(
        "recount3.resource.default_config", return_value=fc
    ):
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
    assert (
        res._cache_root()  # pylint: disable=protected-access
        == cfg.cache_dir
    )


def test_cache_root_calls_default_config(tmp_path: Path) -> None:
    """_cache_root() calls default_config() when config=None."""
    desc = R3ResourceDescription(
        resource_type="data_sources", organism="human"
    )
    fc = _fake_cfg(tmp_path)
    with mock.patch(
        "recount3.resource.default_config", return_value=fc
    ):
        res = R3Resource(description=desc)
        root = res._cache_root()  # pylint: disable=protected-access
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
    cp = res._cached_path()  # pylint: disable=protected-access
    assert "sra.gene_sums.SRP014565.G026.gz" in cp.name


# ===========================================================================
# R3Resource._ensure_cached
# ===========================================================================


def test_ensure_cached_disable_raises(cfg: Config) -> None:
    """_ensure_cached('disable') raises ValueError."""
    res = _make("data_sources", cfg, organism="human")
    with pytest.raises(ValueError, match="disable.*mode"):
        res._ensure_cached(  # pylint: disable=protected-access
            mode="disable", chunk_size=1024
        )


def test_ensure_cached_enable_calls_ensure_cached_url(
    cfg: Config,
) -> None:
    """_ensure_cached('enable') delegates to _ensure_cached_url."""
    res = _make("data_sources", cfg, organism="human")
    with mock.patch(
        "recount3.resource._ensure_cached_url",
        return_value=Path("/fake"),
    ) as mock_ecu:
        result = res._ensure_cached(  # pylint: disable=protected-access
            mode="enable", chunk_size=512
        )
    mock_ecu.assert_called_once()
    assert result == Path("/fake")


def test_ensure_cached_update_calls_download_to_file(
    cfg: Config,
) -> None:
    """_ensure_cached('update') calls download_to_file under the lock."""
    res = _make("data_sources", cfg, organism="human")
    with mock.patch(
        "recount3.resource.download_to_file"
    ) as mock_dl:
        res._ensure_cached(  # pylint: disable=protected-access
            mode="update", chunk_size=256
        )
    mock_dl.assert_called_once()
    assert mock_dl.call_args.kwargs["chunk_size"] == 256


def test_ensure_cached_unknown_mode_raises(cfg: Config) -> None:
    """_ensure_cached raises ValueError for an unknown mode string."""
    res = _make("data_sources", cfg, organism="human")
    with pytest.raises(ValueError, match="Unknown cache mode"):
        res._ensure_cached(  # pylint: disable=protected-access
            mode="bogus",  # type: ignore[arg-type]
            chunk_size=1024,
        )


def test_ensure_cached_config_none_calls_default_config(
    tmp_path: Path,
) -> None:
    """_ensure_cached calls default_config() when config=None."""
    desc = R3ResourceDescription(
        resource_type="data_sources", organism="human"
    )
    fc = _fake_cfg(tmp_path)
    with mock.patch(
        "recount3.resource.default_config", return_value=fc
    ):
        res = R3Resource(description=desc)
        with mock.patch("recount3.resource.download_to_file"):
            res._ensure_cached(  # pylint: disable=protected-access
                mode="update", chunk_size=512
            )


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
    desc = R3ResourceDescription(
        resource_type="data_sources", organism="human"
    )
    res = R3Resource(description=desc, config=disabled_cfg)
    dest_dir = tmp_path / "out"
    dest_dir.mkdir()
    with mock.patch(
        "recount3.resource.download_to_file"
    ) as mock_dl:
        res.download(str(dest_dir))
    mock_dl.assert_called_once()


def test_download_config_none_calls_default_config(
    tmp_path: Path,
) -> None:
    """download() calls default_config() when config=None."""
    desc = R3ResourceDescription(
        resource_type="data_sources", organism="human"
    )
    fc = _fake_cfg(tmp_path)
    (tmp_path / "out").mkdir()
    with mock.patch(
        "recount3.resource.default_config", return_value=fc
    ):
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
    with mock.patch(
        "recount3.resource.download_stream_to_zip"
    ) as mock_ds:
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
    with mock.patch(
        "recount3.resource.download_to_file"
    ) as mock_dl:
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
        with mock.patch(
            "recount3.resource._hardlink_or_copy"
        ) as mock_hl:
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
        with mock.patch(
            "recount3.resource._hardlink_or_copy"
        ) as mock_hl:
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
        with mock.patch(
            "recount3.resource._hardlink_or_copy"
        ) as mock_hl:
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
# R3Resource.load — in-memory cache
# ===========================================================================


def test_load_returns_cached_without_force(cfg: Config) -> None:
    """load() returns _cached_data without disk I/O when force=False."""
    res = _make("data_sources", cfg, organism="human")
    sentinel = object()
    res._cached_data = sentinel  # pylint: disable=protected-access
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
    res._cached_data = "stale"  # pylint: disable=protected-access
    _seed(res, _META_MD_GZ)
    result = res.load(force=True)
    assert isinstance(result, pd.DataFrame)
    assert result is not "stale"  # noqa: F632


# ===========================================================================
# R3Resource.load — bigwig_files
# ===========================================================================


def test_load_bigwig_returns_open_bigwigfile(cfg: Config) -> None:
    """load() for bigwig_files returns an open BigWigFile instance."""
    pytest.importorskip("pyBigWig")
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


def test_load_bigwig_second_call_returns_cached(cfg: Config) -> None:
    """Second load() without force=True returns the same BigWigFile."""
    pytest.importorskip("pyBigWig")
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


def test_load_counts_uses_feature_id_index(
    cfg: Config, tmp_path: Path
) -> None:
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


def _seed_jxn_mm_and_id(
    res: R3Resource, cfg: Config
) -> None:
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
    with mock.patch(
        "recount3.resource.default_config", return_value=fc
    ):
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
    desc = R3ResourceDescription(
        resource_type="data_sources", organism="human"
    )
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
    desc = R3ResourceDescription(
        resource_type="data_sources", organism="human"
    )
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
    desc = R3ResourceDescription(
        resource_type="data_sources", organism="human"
    )
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
    assert (
        _make("data_sources", cfg, organism="human").get_loaded() is None
    )


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


def test_clear_loaded_closes_bigwig(cfg: Config) -> None:
    """clear_loaded() calls close() on a cached BigWigFile."""
    pytest.importorskip("pyBigWig")
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
    res._cached_data = mock_bw  # pylint: disable=protected-access
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
