from dnora.read.generic import ConstantData
from dnora.grid import Grid, TriGrid
from dnora.type_manager.dnora_types import DnoraDataType
import numpy as np
import pandas as pd


def test_gridded_data_on_gridded_grid_cartesian():
    grid = Grid(x=(0, 10), y=(10, 20))
    grid.set_spacing(nx=20, ny=30)
    grid.utm.set((34, "V"))
    reader = ConstantData()
    start_time = "2020-01-01 00:00"
    end_time = "2020-01-01 23:00"

    coord_dict, data_dict, meta_dict = reader(
        DnoraDataType.GRID,
        grid,
        start_time,
        end_time,
        source=None,
        folder=None,
    )
    assert len(coord_dict.keys()) == 2
    assert coord_dict.get("lon") is None
    assert coord_dict.get("lat") is None
    np.testing.assert_almost_equal(coord_dict.get("x"), grid.x())
    np.testing.assert_almost_equal(coord_dict.get("y"), grid.y())

    assert len(data_dict.keys()) == 1
    assert data_dict.get("topo").shape == grid.size()
    np.testing.assert_almost_equal(data_dict.get("topo"), 1.0)

    assert meta_dict == {"zone_number": 34, "zone_letter": "V"}


def test_gridded_data_on_gridded_grid_sperical():
    grid = Grid(lon=(0, 10), lat=(10, 20))
    grid.set_spacing(nx=20, ny=30)
    reader = ConstantData(vars={"u": 2.0, "v": 5.0})
    start_time = "2020-01-01 00:00"
    hours = 24
    end_time = pd.to_datetime(start_time) + pd.Timedelta(hours=hours - 1)

    coord_dict, data_dict, meta_dict = reader(
        DnoraDataType.WIND,
        grid,
        start_time,
        end_time,
        source=None,
        folder=None,
    )
    assert len(coord_dict.keys()) == 3
    assert coord_dict.get("x") is None
    assert coord_dict.get("y") is None
    np.testing.assert_almost_equal(coord_dict.get("lon"), grid.lon())
    np.testing.assert_almost_equal(coord_dict.get("lat"), grid.lat())

    assert list(coord_dict.get("time").strftime("%Y-%m-%d %h%m")) == list(
        pd.date_range(start_time, end_time, freq="h").strftime("%Y-%m-%d %h%m")
    )

    assert len(data_dict.keys()) == 2
    assert data_dict.get("u").shape == (hours, grid.ny(), grid.nx())
    assert data_dict.get("v").shape == (hours, grid.ny(), grid.nx())
    np.testing.assert_almost_equal(data_dict.get("u"), 2.0)
    np.testing.assert_almost_equal(data_dict.get("v"), 5.0)
    assert meta_dict == {}


def test_point_data_on_gridded_grid_cartesian():
    grid = Grid(x=(0, 10), y=(10, 20))
    grid.set_spacing(nx=5, ny=10)
    reader = ConstantData(vars={"spec": 2.0}, peaks={"freq": None, "dirs": None})
    start_time = "2020-01-01 00:00"
    hours = 24
    end_time = pd.to_datetime(start_time) + pd.Timedelta(hours=hours - 1)

    coord_dict, data_dict, meta_dict = reader(
        DnoraDataType.SPECTRA1D,
        grid,
        start_time,
        end_time,
        source=None,
        folder=None,
    )
    assert len(coord_dict.keys()) == 4
    assert coord_dict.get("lon") is None
    assert coord_dict.get("lat") is None
    x, y = grid.xy()
    np.testing.assert_almost_equal(coord_dict.get("x"), x)
    np.testing.assert_almost_equal(coord_dict.get("y"), y)

    assert list(coord_dict.get("time").strftime("%Y-%m-%d %h%m")) == list(
        pd.date_range(start_time, end_time, freq="h").strftime("%Y-%m-%d %h%m")
    )
    assert set(data_dict.keys()) == {"spr", "dirm", "spec"}
    assert data_dict.get("spec").shape == (hours, grid.ny() * grid.nx(), 10)
    np.testing.assert_almost_equal(data_dict.get("spec"), 2.0)
    # np.testing.assert_almost_equal(data_dict.get("v"), 5.0)

    assert meta_dict == {}


def test_point_data_on_gridded_grid_spherical():
    grid = Grid(lon=(0, 10), lat=(10, 20))
    grid.set_spacing(nx=5, ny=10)
    reader = ConstantData(vars={"spec": 2.0}, peaks={"freq": None, "dirs": None})
    start_time = "2020-01-01 00:00"
    hours = 24
    end_time = pd.to_datetime(start_time) + pd.Timedelta(hours=hours - 1)

    coord_dict, data_dict, meta_dict = reader(
        DnoraDataType.SPECTRA1D,
        grid,
        start_time,
        end_time,
        source=None,
        folder=None,
    )
    assert len(coord_dict.keys()) == 4
    assert coord_dict.get("x") is None
    assert coord_dict.get("y") is None
    lon, lat = grid.lonlat()
    np.testing.assert_almost_equal(coord_dict.get("lon"), lon)
    np.testing.assert_almost_equal(coord_dict.get("lat"), lat)

    assert list(coord_dict.get("time").strftime("%Y-%m-%d %h%m")) == list(
        pd.date_range(start_time, end_time, freq="h").strftime("%Y-%m-%d %h%m")
    )
    assert set(data_dict.keys()) == {"spr", "dirm", "spec"}
    assert data_dict.get("spec").shape == (hours, grid.ny() * grid.nx(), 10)
    np.testing.assert_almost_equal(data_dict.get("spec"), 2.0)
    # np.testing.assert_almost_equal(data_dict.get("v"), 5.0)

    assert meta_dict == {}


def test_gridded_data_on_point_grid_cartesian():
    grid = TriGrid(x=np.linspace(0, 10, 11), y=np.linspace(10, 20, 11))
    grid.utm.set((32, "W"))
    reader = ConstantData(vars={"eta": -5.0})
    start_time = "2020-01-01 00:00"
    end_time = "2020-01-01 23:00"

    coord_dict, data_dict, meta_dict = reader(
        DnoraDataType.WATERLEVEL,
        grid,
        start_time,
        end_time,
        source=None,
        folder=None,
    )
    assert set(coord_dict.keys()) == {"time", "x", "y"}
    assert coord_dict.get("lon") is None
    assert coord_dict.get("lat") is None

    np.testing.assert_almost_equal(coord_dict.get("x"), np.linspace(0, 10, 4))
    np.testing.assert_almost_equal(coord_dict.get("y"), np.linspace(10, 20, 4))

    assert len(data_dict.keys()) == 1
    assert data_dict.get("eta").shape == (len(coord_dict["time"]), 4, 4)
    np.testing.assert_almost_equal(data_dict.get("eta"), -5.0)

    assert meta_dict == {"zone_number": 32, "zone_letter": "W"}


def test_point_data_on_point_grid_sperical():
    grid = TriGrid(lon=np.linspace(0, 10, 11), lat=np.linspace(10, 20, 11))
    reader = ConstantData(
        vars={"spec": 0.1},
        coords={"freq": np.linspace(0.2, 0.5, 4)},
        peaks={"freq": None, "dirs": None},
    )
    start_time = "2020-01-01 00:00"
    end_time = "2020-01-01 23:00"

    coord_dict, data_dict, meta_dict = reader(
        DnoraDataType.SPECTRA,
        grid,
        start_time,
        end_time,
        source=None,
        folder=None,
    )
    assert set(coord_dict.keys()) == {"time", "lon", "lat", "freq", "dirs"}
    assert coord_dict.get("x") is None
    assert coord_dict.get("y") is None

    np.testing.assert_almost_equal(coord_dict.get("lon"), np.linspace(0, 10, 11))
    np.testing.assert_almost_equal(coord_dict.get("lat"), np.linspace(10, 20, 11))
    np.testing.assert_almost_equal(coord_dict.get("freq"), np.linspace(0.2, 0.5, 4))
    np.testing.assert_almost_equal(coord_dict.get("dirs"), np.linspace(0, 350, 36))
    assert len(data_dict.keys()) == 1
    assert data_dict.get("spec").shape == (len(coord_dict["time"]), 11, 4, 36)
    np.testing.assert_almost_equal(data_dict.get("spec"), 0.1)

    assert meta_dict == {}


def test_gridded_data_on_gridded_grid_force_spherical():
    grid = Grid(x=(0, 10), y=(10, 20))
    grid.utm.set((33, "W"))
    grid.set_spacing(nx=20, ny=30)
    reader = ConstantData()
    start_time = "2020-01-01 00:00"
    end_time = "2020-01-01 23:00"

    coord_dict, data_dict, meta_dict = reader(
        DnoraDataType.GRID,
        grid,
        start_time,
        end_time,
        source=None,
        folder=None,
        force_type="spherical",
    )
    assert len(coord_dict.keys()) == 2
    assert coord_dict.get("x") is None
    assert coord_dict.get("y") is None
    np.testing.assert_almost_equal(coord_dict.get("lon"), grid.lon())
    np.testing.assert_almost_equal(coord_dict.get("lat"), grid.lat())

    assert len(data_dict.keys()) == 1
    assert data_dict.get("topo").shape == grid.size()
    np.testing.assert_almost_equal(data_dict.get("topo"), 1.0)

    assert meta_dict == {}


def test_point_data_on_gridded_grid_force_cartesian():
    grid = Grid(lon=(0, 10), lat=(10, 20))
    grid.set_spacing(nx=5, ny=10)
    reader = ConstantData(vars={"spec": 2.0}, peaks={"freq": None, "dirs": None})
    start_time = "2020-01-01 00:00"
    hours = 24
    end_time = pd.to_datetime(start_time) + pd.Timedelta(hours=hours - 1)

    coord_dict, data_dict, meta_dict = reader(
        DnoraDataType.SPECTRA1D,
        grid,
        start_time,
        end_time,
        source=None,
        folder=None,
        force_type="cartesian",
    )
    assert len(coord_dict.keys()) == 4
    assert coord_dict.get("lon") is None
    assert coord_dict.get("lat") is None
    x, y = grid.xy()
    np.testing.assert_almost_equal(coord_dict.get("x"), x)
    np.testing.assert_almost_equal(coord_dict.get("y"), y)

    assert list(coord_dict.get("time").strftime("%Y-%m-%d %h%m")) == list(
        pd.date_range(start_time, end_time, freq="h").strftime("%Y-%m-%d %h%m")
    )
    assert set(data_dict.keys()) == {"spr", "dirm", "spec"}
    assert data_dict.get("spec").shape == (hours, grid.ny() * grid.nx(), 10)
    np.testing.assert_almost_equal(data_dict.get("spec"), 2.0)
    # np.testing.assert_almost_equal(data_dict.get("v"), 5.0)

    assert meta_dict == {
        "zone_number": grid.utm.zone()[0],
        "zone_letter": grid.utm.zone()[1],
    }
