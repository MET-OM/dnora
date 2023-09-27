from dnora import grd, mdl, bnd, wsr, pick
import numpy as np

def test_import_constant_waveseries_one_point():
    grid = grd.Grid(lon=5, lat=60)
    model = mdl.ModelRun(grid, start_time = '2020-01-01 00:00', end_time='2020-01-02 00:00')
    
    model.import_boundary(bnd.read.ConstantBoundary(), point_picker=pick.TrivialPicker())
    model.boundary_to_waveseries()
  
    assert model.waveseries().size() == (25,1)

    np.testing.assert_almost_equal(np.mean(model.waveseries().hs()), 1.94162591)