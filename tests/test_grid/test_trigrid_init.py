from dnora.grid import TriGrid
import numpy as np


def test_init_trivial():
    grid = TriGrid(lon=(1, 2), lat=(0, 3), ntriang=range(5), corner=range(3))
    assert grid.nx() == 2
    assert grid.ny() == 2
    assert grid.size(coord_group="grid") == (2,)
    np.testing.assert_array_almost_equal(grid.lon(), np.array([1, 2]))
    np.testing.assert_array_almost_equal(grid.lat(), np.array([0, 3]))
    np.testing.assert_array_almost_equal(grid.inds(), np.array([0, 1]))
    np.testing.assert_array_almost_equal(grid.ntriang(), np.array([0, 1, 2, 3, 4]))
    np.testing.assert_array_almost_equal(grid.corner(), np.array([0, 1, 2]))
    np.testing.assert_array_almost_equal(
        grid.lonlat(), (np.array([1, 2]), np.array([0, 3]))
    )


def test_init_one_point_in_lat():
    grid = TriGrid(lon=(3, 5), lat=(0, 0), ntriang=range(5), corner=range(3))
    assert grid.nx() == 2
    assert grid.ny() == 2
    assert grid.size(coord_group="grid") == (2,)
    np.testing.assert_array_almost_equal(grid.lon(), np.array([3, 5]))
    np.testing.assert_array_almost_equal(grid.lat(), np.array([0, 0]))
    np.testing.assert_array_almost_equal(grid.inds(), np.array([0, 1]))
    np.testing.assert_array_almost_equal(
        grid.lonlat(), (np.array([3, 5]), np.array([0, 0]))
    )


def test_init_one_point_in_lon():
    grid = TriGrid(lon=(0, 0), lat=(3, 5), ntriang=range(5), corner=range(3))
    assert grid.nx() == 2
    assert grid.ny() == 2
    assert grid.size(coord_group="grid") == (2,)
    np.testing.assert_array_almost_equal(grid.lat(), np.array([3, 5]))
    np.testing.assert_array_almost_equal(grid.lon(), np.array([0, 0]))
    np.testing.assert_array_almost_equal(grid.inds(), np.array([0, 1]))
    np.testing.assert_array_almost_equal(
        grid.lonlat(), (np.array([0, 0]), np.array([3, 5]))
    )


def test_init_one_point():
    grid = TriGrid(lon=0, lat=3, ntriang=range(5), corner=range(3))
    assert grid.nx() == 1
    assert grid.ny() == 1
    assert grid.size(coord_group="grid") == (1,)
    np.testing.assert_array_almost_equal(grid.lat(), np.array([3]))
    np.testing.assert_array_almost_equal(grid.lon(), np.array([0]))
    np.testing.assert_array_almost_equal(grid.inds(), np.array([0]))
    np.testing.assert_array_almost_equal(grid.lonlat(), (np.array([0]), np.array([3])))


def test_init_long():
    grid = TriGrid(
        lon=[0, 2, 4, 5, 6, 6],
        lat=[3, 1, 2, 3, 4, 5],
        ntriang=range(5),
        corner=range(3),
    )
    assert grid.nx() == 6
    assert grid.ny() == 6
    assert grid.size(coord_group="grid") == (6,)
    np.testing.assert_array_almost_equal(grid.lat(), np.array([3, 1, 2, 3, 4, 5]))
    np.testing.assert_array_almost_equal(grid.lon(), np.array([0, 2, 4, 5, 6, 6]))
    np.testing.assert_array_almost_equal(grid.inds(), np.array([0, 1, 2, 3, 4, 5]))
    np.testing.assert_array_almost_equal(
        grid.lonlat(), (np.array([0, 2, 4, 5, 6, 6]), np.array([3, 1, 2, 3, 4, 5]))
    )
