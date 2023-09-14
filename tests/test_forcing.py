from dnora import grd, mdl, wnd
import numpy as np

def test_import_constant_forcing_one_point():
    grid = grd.Grid(lon=5, lat=60)
    model = mdl.ModelRun(grid, start_time = '2020-01-01 00:00', end_time='2020-01-02 00:00')
    
    model.import_forcing(wnd.read.ConstantForcing(u=1., v=2.))
    
    assert model.forcing().size() == (25,1,1)

    np.testing.assert_almost_equal(np.mean(model.forcing().u()), 1)
    np.testing.assert_almost_equal(np.mean(model.forcing().v()), 2)
    np.testing.assert_almost_equal(np.mean(model.forcing().magnitude()), (2**2+1**2)**0.5)


def test_import_constant_forcing():
    grid = grd.Grid(lon=(5,6), lat=(60,61))
    grid.set_spacing(nx=5,ny=10)

    model = mdl.ModelRun(grid, start_time = '2020-01-01 00:00', end_time='2020-01-02 00:00')
    
    model.import_forcing(wnd.read.ConstantForcing(u=1., v=2.))
    
    assert model.forcing().size() == (25,10,5)

    np.testing.assert_almost_equal(np.mean(model.forcing().u()), 1)
    np.testing.assert_almost_equal(np.mean(model.forcing().v()), 2)
    np.testing.assert_almost_equal(np.mean(model.forcing().magnitude()), (2**2+1**2)**0.5)