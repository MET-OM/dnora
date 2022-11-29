from dnora.grd import Grid
import numpy as np
def test_nx_ny_cartesian():
    grid = Grid(x=(-2,2), y=(-3,3))
    grid.set_spacing(nx=5, ny=7)
    assert grid.nx() == 5
    assert grid.ny() == 7
    assert grid.size() == (7,5)
    np.testing.assert_array_almost_equal(grid.x(),np.array([-2,-1,0,1,2]))
    np.testing.assert_array_almost_equal(grid.y(),np.array([-3,-2,-1,0,1,2,3]))

def test_nx_ny_spherical():
    grid = Grid(lon=(-2,2), lat=(0,3))
    grid.set_spacing(nx=5, ny=4)
    assert grid.nx() == 5
    assert grid.ny() == 4
    assert grid.size() == (4,5)
    np.testing.assert_array_almost_equal(grid.lon(),np.array([-2,-1,0,1,2]))
    np.testing.assert_array_almost_equal(grid.lat(),np.array([0,1,2,3]))

def test_dx_dy_cartesian():
    grid = Grid(x=(-1,1), y=(-3,3))
    grid.set_spacing(dx=0.5, dy=3)
    assert grid.nx() == 5
    assert grid.ny() == 3
    assert grid.size() == (3,5)
    np.testing.assert_array_almost_equal(grid.x(),np.array([-1,-0.5,0,0.5,1]))
    np.testing.assert_array_almost_equal(grid.y(),np.array([-3,0,3]))

def test_dm_cartesian():
    grid = Grid(x=(-1,1), y=(-2,2))
    grid.set_spacing(dm=0.5)
    assert grid.nx() == 5
    assert grid.ny() == 9
    assert grid.size() == (9,5)
    np.testing.assert_array_almost_equal(grid.x(),np.array([-1,-0.5,0,0.5,1]))
    np.testing.assert_array_almost_equal(grid.y(),np.array([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]))

def test_dlon_dlat_spherical():
    grid = Grid(lon=(-1,1), lat=(0,3))
    grid.set_spacing(dlon=0.5, dlat=1.5)
    assert grid.nx() == 5
    assert grid.ny() == 3
    assert grid.size() == (3,5)
    np.testing.assert_array_almost_equal(grid.lon(),np.array([-1,-0.5,0,0.5,1]))
    np.testing.assert_array_almost_equal(grid.lat(),np.array([0,1.5,3]))

def test_dx_dy_spherical():
    grid = Grid(lon=(4,5), lat=(60,61))
    grid.set_spacing(dx=1110, dy=1110)
    breakpoint()
    np.testing.assert_array_almost_equal(grid.dlat(),0.01,decimal=3)
    np.testing.assert_array_almost_equal(grid.dlon(),0.02,decimal=3)
    np.testing.assert_array_almost_equal(grid.dy(),1115,decimal=0)
    np.testing.assert_array_almost_equal(grid.dx(),1110,decimal=0)
    assert grid.nx() == 1/0.02
    assert grid.ny() == 1/0.01

def test_dlon_dlat_cartesian():
    grid = Grid(x=(0,150_000),y=(6_700_000,6_800_000))
    grid.set_spacing(dlon=0.02, dlat=0.01)
    np.testing.assert_array_almost_equal(grid.dlat(),0.01,decimal=3)
    np.testing.assert_array_almost_equal(grid.dlon(),0.02,decimal=3)
    np.testing.assert_array_almost_equal(grid.dy(),1123,decimal=0)
    np.testing.assert_array_almost_equal(grid.dx(),1103,decimal=0)
    assert grid.nx() == 137
    assert grid.ny() == 90

def test_dlon_dlat_spherical_floating():
    grid = Grid(lon=(-1,0.999), lat=(0,2.999))
    grid.set_spacing(dlon=0.5, dlat=1.5, floating_edge=True)
    assert grid.nx() == 5
    assert grid.ny() == 3
    assert grid.size() == (3,5)
    np.testing.assert_array_almost_equal(grid.lon(),np.array([-1,-0.5,0,0.5,1]))
    np.testing.assert_array_almost_equal(grid.lat(),np.array([0,1.5,3]))

def test_dx_dy_cartesian_floating():
    grid = Grid(x=(-1,0.999), y=(-3,2.999))
    grid.set_spacing(dx=0.5, dy=3, floating_edge=True)
    assert grid.nx() == 5
    assert grid.ny() == 3
    assert grid.size() == (3,5)
    np.testing.assert_array_almost_equal(grid.x(),np.array([-1,-0.5,0,0.5,1]))
    np.testing.assert_array_almost_equal(grid.y(),np.array([-3,0,3]))
