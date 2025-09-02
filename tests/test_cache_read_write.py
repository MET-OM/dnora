import dnora as dn
import os
import shutil
import glob
import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
import stat


def handle_remove_readonly(func, path, exc_info):
    """Clear the read-only bit and reattempt the removal."""
    os.chmod(path, stat.S_IWRITE)  # Change to writable
    func(path)


def test_write_wind():
    """Write wind data to cache"""
    grid = dn.grid.Grid(lon=(4, 11), lat=(60, 65))
    grid.set_spacing(dlon=1, dlat=1)

    model = dn.modelrun.ModelRun(
        grid, start_time="2020-01-31 00:00", end_time="2020-02-02 01:00"
    )
    wind_cache = Path("wind_cache")
    if wind_cache.exists():
        shutil.rmtree(wind_cache, onerror=handle_remove_readonly)

    # Automatically uses expansion factor 1.2, so tiles 55-60 also included
    model.import_wind(
        dn.read.generic.ConstantData(debug_cache=True), u=1, v=2, write_cache=True
    )
    assert sorted(wind_cache.glob("constantdata/*")) == [
        wind_cache / "constantdata/2020"
    ]
    assert set(wind_cache.glob("constantdata/2020/*")) == {
        wind_cache / "constantdata/2020/01",
        wind_cache / "constantdata/2020/02",
    }

    jan_files = list(wind_cache.glob("constantdata/2020/01/*"))
    days = 1
    lat_tiles = 3  # 55-60, 60-65, 65-70
    lon_tiles = 3  # 0-5, 5-10, 10-15
    assert len(jan_files) == days * lat_tiles * lon_tiles

    feb_files = list(wind_cache.glob("constantdata/2020/02/*"))
    days = 2
    lat_tiles = 3  # 55-60, 60-65, 65-70
    lon_tiles = 3  # 0-5, 5-10, 10-15
    assert len(feb_files) == days * lat_tiles * lon_tiles


def test_content_of_cached_wind_files():
    """Test that the cached files actually have the data they should."""
    wind_cache = Path("wind_cache")
    all_files = list(wind_cache.glob("constantdata/2020/01/*")) + list(
        wind_cache.glob("constantdata/2020/02/*")
    )
    ds = xr.open_mfdataset(all_files)
    np.testing.assert_array_almost_equal(ds.lon.values, np.arange(0, 15, 1))
    np.testing.assert_array_almost_equal(ds.lat.values, np.arange(55, 70, 1))
    times = pd.to_datetime(ds.time.values)
    assert len(times) == 3 * 24
    assert times[0].strftime("%Y-%m-%d %H:%M") == "2020-01-31 00:00"
    assert times[-1].strftime("%Y-%m-%d %H:%M") == "2020-02-02 23:00"


def test_write_spectra():
    """Write spectral data to cache"""

    grid = dn.grid.Grid(lon=(4, 11), lat=(60, 65))
    grid.set_spacing(dlon=1, dlat=1)

    model = dn.modelrun.ModelRun(
        grid, start_time="2020-01-31 00:00", end_time="2020-02-02 01:00"
    )

    spectra_cache = Path("spectra_cache")
    if spectra_cache.exists():
        shutil.rmtree(spectra_cache, onerror=handle_remove_readonly)

    # Automatically uses expansion_factor=1.5, so:
    # lat tiles 65-70 and 70-75 included
    # lon tiles -5-0 included
    # However, no data for 70-75 of -5-0
    model.import_spectra(
        dn.read.generic.ConstantData(debug_cache=True), write_cache=True
    )

    # Cross-platform assertion for the constantdata directory
    assert sorted(spectra_cache.glob("constantdata/*")) == [
        spectra_cache / "constantdata/2020"
    ]

    # Cross-platform assertion for the subdirectories under 2020
    assert set(spectra_cache.glob("constantdata/2020/*")) == {
        spectra_cache / "constantdata/2020/01",
        spectra_cache / "constantdata/2020/02",
    }

    # Check the number of files in January (constantdata/2020/01)
    jan_files = list(spectra_cache.glob("constantdata/2020/01/*"))
    days = 1
    lat_tiles = 3  # 55-60, 60-65, 65-70
    lon_tiles = 3  # 0-5, 5-10, 10-15
    assert len(jan_files) == days * lat_tiles * lon_tiles

    # Check the number of files in February (constantdata/2020/02)
    feb_files = list(spectra_cache.glob("constantdata/2020/02/*"))
    days = 2
    lat_tiles = 3  # 55-60, 60-65, 65-70
    lon_tiles = 3  # 0-5, 5-10, 10-15
    assert len(feb_files) == days * lat_tiles * lon_tiles


def test_content_of_cached_spectral_files():
    """Test that the cached files actually have the data they should. Can't use open_mfdataset since the common cariable is inds"""
    grid = dn.grid.Grid(lon=(0, 14), lat=(55, 69))
    grid.set_spacing(dlon=1, dlat=1)

    model = dn.modelrun.ModelRun(
        grid, start_time="2020-01-31 00:00", end_time="2020-02-02 23:00"
    )

    model.import_spectra(
        dn.read.generic.ConstantData(debug_cache=True),
        read_cache=True,
        point_picker=dn.pick.Area(),
        expansion_factor=1,
    )

    ds = model.spectra().ds()
    np.testing.assert_array_almost_equal(np.unique(ds.lon.values), np.arange(0, 15, 1))
    np.testing.assert_array_almost_equal(np.unique(ds.lat.values), np.arange(55, 70, 1))

    times = pd.to_datetime(ds.time.values)

    assert len(times) == 3 * 24
    assert times[0].strftime("%Y-%m-%d %H:%M") == "2020-01-31 00:00"
    assert times[-1].strftime("%Y-%m-%d %H:%M") == "2020-02-02 23:00"


def test_read_subset_of_cached_wind_data():
    """Test that we can read in a subset of the cached wind data"""
    grid = dn.grid.Grid(lon=(0, 11), lat=(60, 65))
    grid.set_spacing(dlon=1, dlat=1)

    model = dn.modelrun.ModelRun(
        grid, start_time="2020-01-31 15:00", end_time="2020-02-02 14:00"
    )

    model.import_wind(
        dn.read.generic.ConstantData(debug_cache=True),
        read_cache=True,
        expansion_factor=1,
    )

    ds = model.wind().ds()
    np.testing.assert_array_almost_equal(ds.lon.values, np.arange(0, 12, 1))
    np.testing.assert_array_almost_equal(ds.lat.values, np.arange(60, 66, 1))

    times = pd.to_datetime(ds.time.values)

    assert len(times) == 2 * 24
    assert times[0].strftime("%Y-%m-%d %H:%M") == "2020-01-31 15:00"
    assert times[-1].strftime("%Y-%m-%d %H:%M") == "2020-02-02 14:00"

    # We should automatically use an expansion factor 1.2
    model.import_wind(
        dn.read.generic.ConstantData(debug_cache=True),
        read_cache=True,
        dont_patch=True,
    )

    ds = model.wind().ds()
    np.testing.assert_array_almost_equal(np.unique(ds.lon.values), np.arange(0, 13, 1))
    np.testing.assert_array_almost_equal(np.unique(ds.lat.values), np.arange(60, 66, 1))

    times = pd.to_datetime(ds.time.values)

    assert len(times) == 2 * 24
    assert times[0].strftime("%Y-%m-%d %H:%M") == "2020-01-31 15:00"
    assert times[-1].strftime("%Y-%m-%d %H:%M") == "2020-02-02 14:00"

    if os.path.isdir("wind_cache"):
        shutil.rmtree("wind_cache", onerror=handle_remove_readonly)


def test_read_subset_of_cached_spectral_data():
    """Test that we can read in a subset of the cached wind data"""
    grid = dn.grid.Grid(lon=(0, 11), lat=(60, 65))
    grid.set_spacing(dlon=1, dlat=1)

    model = dn.modelrun.ModelRun(
        grid, start_time="2020-01-31 15:00", end_time="2020-02-02 14:00"
    )

    model.import_spectra(
        dn.read.generic.ConstantData(debug_cache=True),
        read_cache=True,
        point_picker=dn.pick.Area(),
        expansion_factor=1,
    )

    ds = model.spectra().ds()
    np.testing.assert_array_almost_equal(np.unique(ds.lon.values), np.arange(0, 12, 1))
    np.testing.assert_array_almost_equal(np.unique(ds.lat.values), np.arange(60, 66, 1))

    times = pd.to_datetime(ds.time.values)

    assert len(times) == 2 * 24

    assert times[0].strftime("%Y-%m-%d %H:%M") == "2020-01-31 15:00"
    assert times[-1].strftime("%Y-%m-%d %H:%M") == "2020-02-02 14:00"

    # We should automatically use an expansion factor 1.5
    model.import_spectra(
        dn.read.generic.ConstantData(debug_cache=True),
        read_cache=True,
        point_picker=dn.pick.Area(),
        dont_patch=True,
    )

    ds = model.spectra().ds()

    np.testing.assert_array_almost_equal(np.unique(ds.lon.values), np.arange(0, 14, 1))
    np.testing.assert_array_almost_equal(np.unique(ds.lat.values), np.arange(59, 67, 1))

    times = pd.to_datetime(ds.time.values)

    assert len(times) == 2 * 24

    assert times[0].strftime("%Y-%m-%d %H:%M") == "2020-01-31 15:00"
    assert times[-1].strftime("%Y-%m-%d %H:%M") == "2020-02-02 14:00"

    if os.path.isdir("spectra_cache"):
        shutil.rmtree("spectra_cache", onerror=handle_remove_readonly)
