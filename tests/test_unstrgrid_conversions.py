from dnora.grd import UnstrGrid
import numpy as np
import utm
def test_utm_conversion():
    lon = np.array([3,4,5])
    lat = np.array([60,60.5,61])
    x, y, zone, letter = utm.from_latlon(lat, lon)

    grid = UnstrGrid(lon=lon, lat=lat)
    assert grid.nx() == 3
    assert grid.ny() == 3
    assert grid.size() == (3,)
    np.testing.assert_array_almost_equal(grid.lon(),lon)
    np.testing.assert_array_almost_equal(grid.lat(),lat)
    np.testing.assert_array_almost_equal(grid.inds(),np.array([0,1,2]))
    np.testing.assert_array_almost_equal(grid.lonlat(),(lon ,lat))
    np.testing.assert_array_almost_equal(grid.lon(strict=True),lon)
    np.testing.assert_array_almost_equal(grid.lat(strict=True),lat)
    np.testing.assert_array_almost_equal(grid.lonlat(strict=True),(lon ,lat))
    assert grid.x(strict=True) is None
    assert grid.y(strict=True) is None
    assert np.all(grid.xy(strict=True) == (None, None))

    grid.set_utm(zone, letter)
    assert grid.utm()[0] == zone
    assert grid.utm()[1] == letter

    np.testing.assert_array_almost_equal(grid.x(), x)
    np.testing.assert_array_almost_equal(grid.y(), y)
    np.testing.assert_array_almost_equal(grid.xy(), (x, y))



    grid2 = UnstrGrid(x=x, y=y)
    grid2.set_utm(zone, letter)
    assert grid.utm()[0] == zone
    assert grid.utm()[1] == letter

    assert grid2.nx() == 3
    assert grid2.ny() == 3
    assert grid2.size() == (3,)

    assert grid2.lon(strict=True) is None
    assert grid2.lat(strict=True) is None
    assert np.all(grid2.lonlat(strict=True) == (None, None))

    np.testing.assert_array_almost_equal(grid2.x(strict=True), x)
    np.testing.assert_array_almost_equal(grid2.y(strict=True), y)
    np.testing.assert_array_almost_equal(grid2.xy(strict=True), (x, y))

    np.testing.assert_array_almost_equal(grid2.lon(), lon, decimal=5)
    np.testing.assert_array_almost_equal(grid2.lat(), lat, decimal=5)
    np.testing.assert_array_almost_equal(grid2.lonlat(), (lon, lat), decimal=5)
