from dnora import grd, mdl, wlv
import numpy as np

def test_import_constant_waterlevel_one_point():
    grid = grd.Grid(lon=5, lat=60)
    model = mdl.ModelRun(grid, start_time = '2020-01-01 00:00', end_time='2020-01-02 00:00')
    
    model.import_waterlevel(wlv.read.ConstantWaterLevel(waterlevel=1.))
    
    assert model.waterlevel().size() == (25,1,1)

    np.testing.assert_almost_equal(np.mean(model.waterlevel().waterlevel()), 1)


def test_import_constant_waterlevel():
    grid = grd.Grid(lon=(5,6), lat=(60,61))
    grid.set_spacing(nx=5,ny=10)

    model = mdl.ModelRun(grid, start_time = '2020-01-01 00:00', end_time='2020-01-02 00:00')
    
    model.import_waterlevel(wlv.read.ConstantWaterLevel(waterlevel=1.))
    
    assert model.waterlevel().size() == (25,10,5)

    np.testing.assert_almost_equal(np.mean(model.waterlevel().waterlevel()), 1)
    