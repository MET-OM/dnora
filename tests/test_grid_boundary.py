from dnora.grid import Grid
from dnora.grid import mask
import numpy as np


def test_boundary_lonlat():
    grid = Grid(lon=(1, 2), lat=(0, 3))
    grid.set_spacing(nx=10, ny=5)

    grid.set_boundary_points(mask.LonLat(lon=1, lat=0))
    assert grid.boundary_mask()[0, 0]

    lon, lat = grid.boundary_points()
    np.testing.assert_array_almost_equal(lon, np.array([1]))
    np.testing.assert_array_almost_equal(lat, np.array([0]))


def test_boundary_all():
    grid = Grid(lon=(1, 2), lat=(0, 3))
    grid.set_spacing(nx=10, ny=5)

    grid.set_boundary_points(mask.All())
    assert np.all(grid.boundary_mask())

    lon, lat = grid.boundary_points()
    all_lon, all_lat = grid.lonlat()
    np.testing.assert_array_almost_equal(lon, all_lon)
    np.testing.assert_array_almost_equal(lat, all_lat)


def test_boundary_clear():
    grid = Grid(lon=(1, 2), lat=(0, 3))
    grid.set_spacing(nx=10, ny=5)

    grid.set_boundary_points(mask.All())
    assert np.all(grid.boundary_mask())

    grid.set_boundary_points(mask.Clear())
    assert np.all(np.logical_not(grid.boundary_mask()))

    lon, lat = grid.boundary_points()
    np.testing.assert_array_almost_equal(lon, np.array([]))
    np.testing.assert_array_almost_equal(lat, np.array([]))


def test_boundary_edges():
    grid = Grid(lon=(1, 2), lat=(0, 3))
    grid.set_spacing(nx=10, ny=5)

    grid.set_boundary_points(mask.Edges(edges=["w", "s"]))

    assert np.all(grid.boundary_mask()[0, :])
    assert np.all(grid.boundary_mask()[:, 0])

    lon, lat = grid.boundary_points()
    all_lon, all_lat = grid.lonlat()
    bnd_mask = np.logical_or(all_lon < 1.001, all_lat < 0.001)

    all_lon = all_lon[bnd_mask]
    all_lat = all_lat[bnd_mask]

    np.testing.assert_array_almost_equal(lon, all_lon)
    np.testing.assert_array_almost_equal(lat, all_lat)


def test_boundary_midpoint():
    grid = Grid(lon=(1, 2), lat=(0, 3))
    grid.set_spacing(nx=10, ny=5)

    grid.set_boundary_points(mask.MidPoint(edges=["w", "s"]))

    assert grid.boundary_mask()[0, 4]
    assert not grid.boundary_mask()[1, 4]
    assert not grid.boundary_mask()[0, 3]
    assert not grid.boundary_mask()[0, 5]
    assert grid.boundary_mask()[2, 0]
    assert not grid.boundary_mask()[2, 1]
    assert not grid.boundary_mask()[1, 0]
    assert not grid.boundary_mask()[3, 0]

    lon, lat = grid.boundary_points()
    np.testing.assert_array_almost_equal(lon, np.array([grid.lon()[4], 1]))
    np.testing.assert_array_almost_equal(lat, np.array([0, 1.5]))
