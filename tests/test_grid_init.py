from dnora.grid import Grid
import numpy as np


def test_init_trivial():
    grid = Grid(lon=(1, 2), lat=(0, 3))
    assert grid.nx() == 2
    assert grid.ny() == 2
    assert grid.size() == (2, 2)
    np.testing.assert_array_almost_equal(grid.lon(), np.array([1, 2]))
    np.testing.assert_array_almost_equal(grid.lat(), np.array([0, 3]))
    assert grid.inds() is None
    np.testing.assert_array_almost_equal(
        grid.lonlat(), (np.array([1, 2, 1, 2]), np.array([0, 0, 3, 3]))
    )


def test_init_one_point_in_lat():
    grid = Grid(lon=(3, 5), lat=(0, 0))
    assert grid.nx() == 2
    assert grid.ny() == 1
    assert grid.size() == (1, 2)
    np.testing.assert_array_almost_equal(grid.lon(), np.array([3, 5]))
    np.testing.assert_array_almost_equal(grid.lat(), np.array([0]))
    assert grid.inds() is None
    np.testing.assert_array_almost_equal(
        grid.lonlat(), (np.array([3, 5]), np.array([0, 0]))
    )


def test_init_one_point_in_lon():
    grid = Grid(lon=(0, 0), lat=(3, 5))
    assert grid.nx() == 1
    assert grid.ny() == 2
    assert grid.size() == (2, 1)
    np.testing.assert_array_almost_equal(grid.lat(), np.array([3, 5]))
    np.testing.assert_array_almost_equal(grid.lon(), np.array([0]))
    assert grid.inds() is None
    np.testing.assert_array_almost_equal(
        grid.lonlat(), (np.array([0, 0]), np.array([3, 5]))
    )


def test_init_one_point():
    grid = Grid(lon=0, lat=3)
    assert grid.nx() == 1
    assert grid.ny() == 1
    assert grid.size() == (1, 1)
    np.testing.assert_array_almost_equal(grid.lat(), np.array([3]))
    np.testing.assert_array_almost_equal(grid.lon(), np.array([0]))
    assert grid.inds() is None
    np.testing.assert_array_almost_equal(grid.lonlat(), (np.array([0]), np.array([3])))
