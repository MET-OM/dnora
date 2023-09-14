from dnora import grd, mdl, bnd
import numpy as np

def test_import_constant_spectra_one_point():
    grid = grd.Grid(lon=5, lat=60)
    model = mdl.ModelRun(grid, start_time = '2020-01-01 00:00', end_time='2020-01-02 00:00')
    
    model.import_boundary(bnd.read.ConstantBoundary(), point_picker=bnd.pick.TrivialPicker())
    model.boundary_to_spectra()
  
    assert model.spectra().size() == (25,1,10)

    np.testing.assert_almost_equal(np.max(model.spectra().spec()), np.max(model.boundary().spec())*(model.boundary().dd())*np.pi/180)

def test_import_constant_spectra():
    grid = grd.Grid(lon=(5,6), lat=(60,61))
    grid.set_spacing(nx=5,ny=10)

    model = mdl.ModelRun(grid, start_time = '2020-01-01 00:00', end_time='2020-01-02 00:00')
    
    model.import_boundary(bnd.read.ConstantBoundary(), point_picker=bnd.pick.TrivialPicker())
    model.boundary_to_spectra()
    assert model.spectra().size() == (25,50,10)
    
    np.testing.assert_almost_equal(np.max(model.spectra().spec()), np.max(model.boundary().spec())*(model.boundary().dd())*np.pi/180)
