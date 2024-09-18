import dnora as dn
import os
import shutil
import glob
import xarray as xr
import numpy as np
import pandas as pd


def test_write_wind():
    """Write wind data to cache"""
    grid = dn.grid.Grid(lon=(4, 11), lat=(60, 65))
    grid.set_spacing(dlon=1, dlat=1)

    model = dn.modelrun.ModelRun(
        grid, start_time="2020-01-31 00:00", end_time="2020-02-02 01:00"
    )
    if os.path.isdir("wind_cache"):
        shutil.rmtree("wind_cache")

    model.import_wind(
        dn.read.generic.ConstantData(debug_cache=True), u=1, v=2, write_cache=True
    )
    assert glob.glob("wind_cache/constantdata/*") == ["wind_cache/constantdata/2020"]
    assert set(glob.glob("wind_cache/constantdata/2020/*")) == {
        "wind_cache/constantdata/2020/01",
        "wind_cache/constantdata/2020/02",
    }
    jan_files = glob.glob("wind_cache/constantdata/2020/01/*")
    days = 1
    lat_tiles = 2  # 60-65, 65-70
    lon_tiles = 3  # 0-5, 5-10, 10-15
    assert len(jan_files) == days * lat_tiles * lon_tiles

    feb_files = glob.glob("wind_cache/constantdata/2020/02/*")
    days = 2
    lat_tiles = 2  # 60-65, 65-70
    lon_tiles = 3  # 0-5, 5-10, 10-15
    assert len(feb_files) == days * lat_tiles * lon_tiles


def test_content_of_cached_wind_files():
    """Test that the cached files actually have the data they should."""
    all_files = glob.glob("wind_cache/constantdata/2020/01/*") + glob.glob(
        "wind_cache/constantdata/2020/02/*"
    )
    ds = xr.open_mfdataset(all_files)
    np.testing.assert_array_almost_equal(ds.lon.values, np.arange(0, 15, 1))
    np.testing.assert_array_almost_equal(ds.lat.values, np.arange(60, 70, 1))
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
    if os.path.isdir("spectra_cache"):
        shutil.rmtree("spectra_cache")

    model.import_spectra(
        dn.read.generic.ConstantData(debug_cache=True), write_cache=True
    )
    assert glob.glob("spectra_cache/constantdata/*") == [
        "spectra_cache/constantdata/2020"
    ]
    assert set(glob.glob("spectra_cache/constantdata/2020/*")) == {
        "spectra_cache/constantdata/2020/01",
        "spectra_cache/constantdata/2020/02",
    }
    jan_files = glob.glob("spectra_cache/constantdata/2020/01/*")
    days = 1
    lat_tiles = 2  # 60-65, 65-70
    lon_tiles = 3  # 0-5, 5-10, 10-15
    assert len(jan_files) == days * lat_tiles * lon_tiles

    feb_files = glob.glob("spectra_cache/constantdata/2020/02/*")
    days = 2
    lat_tiles = 2  # 60-65, 65-70
    lon_tiles = 3  # 0-5, 5-10, 10-15
    assert len(feb_files) == days * lat_tiles * lon_tiles


def test_content_of_cached_spectral_files():
    """Test that the cached files actually have the data they should. Can't use open_mfdataset since the common cariable is inds"""
    grid = dn.grid.Grid(lon=(0, 14), lat=(60, 69))
    grid.set_spacing(dlon=1, dlat=1)

    model = dn.modelrun.ModelRun(
        grid, start_time="2020-01-31 00:00", end_time="2020-02-02 23:00"
    )

    model.import_spectra(
        dn.read.generic.ConstantData(debug_cache=True),
        read_cache=True,
        point_picker=dn.pick.Area(),
    )

    ds = model.spectra().ds()
    np.testing.assert_array_almost_equal(np.unique(ds.lon.values), np.arange(0, 15, 1))
    np.testing.assert_array_almost_equal(np.unique(ds.lat.values), np.arange(60, 70, 1))

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
    )

    ds = model.wind().ds()
    np.testing.assert_array_almost_equal(np.unique(ds.lon.values), np.arange(0, 13, 1))
    np.testing.assert_array_almost_equal(np.unique(ds.lat.values), np.arange(60, 66, 1))

    times = pd.to_datetime(ds.time.values)

    assert len(times) == 2 * 24
    assert times[0].strftime("%Y-%m-%d %H:%M") == "2020-01-31 15:00"
    assert times[-1].strftime("%Y-%m-%d %H:%M") == "2020-02-02 14:00"

    if os.path.isdir("wind_cache"):
        shutil.rmtree("wind_cache")


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

    # We should automatically use an expansion factor
    model.import_spectra(
        dn.read.generic.ConstantData(debug_cache=True),
        read_cache=True,
        point_picker=dn.pick.Area(),
    )

    ds = model.spectra().ds()
    np.testing.assert_array_almost_equal(np.unique(ds.lon.values), np.arange(0, 14, 1))
    np.testing.assert_array_almost_equal(np.unique(ds.lat.values), np.arange(60, 67, 1))

    times = pd.to_datetime(ds.time.values)

    assert len(times) == 2 * 24

    assert times[0].strftime("%Y-%m-%d %H:%M") == "2020-01-31 15:00"
    assert times[-1].strftime("%Y-%m-%d %H:%M") == "2020-02-02 14:00"

    if os.path.isdir("spectra_cache"):
        shutil.rmtree("spectra_cache")
