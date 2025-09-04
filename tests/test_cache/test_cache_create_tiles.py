import dnora as dn
from dnora.cacher.tiles import create_tiles
import numpy as np


def test_create_tiles():
    area = dn.grid.Grid(lon=(0, 1), lat=(50, 52))
    start_time = "2020-01-01 00:00"
    end_time = "2020-01-02 01:00"
    lons, lats, times = create_tiles(area, start_time, end_time, expansion_factor=1)
    assert times == ["2020-01-01", "2020-01-02"]
    np.testing.assert_almost_equal(lons[0][0], 0)
    np.testing.assert_almost_equal(lons[0][1], 5)
    np.testing.assert_almost_equal(lats[0][0], 50)
    np.testing.assert_almost_equal(lats[0][1], 55)


def test_create_tiles_2x2():
    area = dn.grid.Grid(lon=(0, 6), lat=(50, 56))
    start_time = "2020-01-01 00:00"
    end_time = "2020-01-02 01:00"
    lons, lats, times = create_tiles(area, start_time, end_time, expansion_factor=1)
    assert times == [
        "2020-01-01",
        "2020-01-01",
        "2020-01-01",
        "2020-01-01",
        "2020-01-02",
        "2020-01-02",
        "2020-01-02",
        "2020-01-02",
    ]
    for n in range(len(lons)):
        if n % 2 == 0:
            np.testing.assert_almost_equal(lons[n][0], 0)
            np.testing.assert_almost_equal(lons[n][1], 5)

        else:
            np.testing.assert_almost_equal(lons[n][0], 5)
            np.testing.assert_almost_equal(lons[n][1], 10)

        if n in [0, 1, 4, 5]:
            np.testing.assert_almost_equal(lats[n][0], 50)
            np.testing.assert_almost_equal(lats[n][1], 55)
        else:
            np.testing.assert_almost_equal(lats[n][0], 55)
            np.testing.assert_almost_equal(lats[n][1], 60)


def test_create_tiles_single_point_grid():
    area = dn.grid.Grid(lon=(0), lat=(50))
    start_time = "2020-01-01 00:00"
    end_time = "2020-01-02 01:00"
    lons, lats, times = create_tiles(area, start_time, end_time, expansion_factor=1)
    assert times == ["2020-01-01", "2020-01-02"]
    np.testing.assert_almost_equal(lons[0][0], 0)
    np.testing.assert_almost_equal(lons[0][1], 5)
    np.testing.assert_almost_equal(lats[0][0], 50)
    np.testing.assert_almost_equal(lats[0][1], 55)
