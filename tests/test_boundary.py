from dnora import grd, mdl, bnd, pick
import numpy as np

def test_import_constant_boundary_one_point():
    grid = grd.Grid(lon=5, lat=60)
    model = mdl.ModelRun(grid, start_time = '2020-01-01 00:00', end_time='2020-01-02 00:00')
    
    model.import_boundary(bnd.read.ConstantBoundary(), point_picker=pick.TrivialPicker())
  
    assert model.boundary().size() == (25,1,10,24)
    
    np.testing.assert_almost_equal(np.max(model.boundary().spec()), 1)
    np.testing.assert_almost_equal(np.min(model.boundary().spec()), 0)
    

def test_import_constant_boundary():
    grid = grd.Grid(lon=(5,6), lat=(60,61))
    grid.set_spacing(nx=5,ny=10)

    model = mdl.ModelRun(grid, start_time = '2020-01-01 00:00', end_time='2020-01-02 00:00')
    
    model.import_boundary(bnd.read.ConstantBoundary(), point_picker=pick.TrivialPicker())
  
    assert model.boundary().size() == (25,50,10,24)
    
    np.testing.assert_almost_equal(np.max(model.boundary().spec()), 1)
    np.testing.assert_almost_equal(np.min(model.boundary().spec()), 0)