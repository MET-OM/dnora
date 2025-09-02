import dnora as dn
import xarray as xr
import numpy as np
import pathlib

def test_read_wind():
    grid = dn.grid.Grid(lon=(0, 10), lat=(50,60))
    grid.set_spacing(nx=10, ny=10)

    model = dn.modelrun.Constant(grid, year=2020,month=1, day=1)
    model.import_wind(u=5.0, v=6.0)
    model.wind().ds().to_netcdf('test_wind.nc')
    
    model2 = dn.modelrun.ModelRun(grid, year=2020,month=1, day=1)
    model2.import_wind(dn.read.generic.Netcdf(), filename='test_wind.nc')

    np.testing.assert_array_almost_equal(model.wind().lon(), model2.wind().lon())
    np.testing.assert_array_almost_equal(model.wind().lat(), model2.wind().lat())
    np.testing.assert_array_almost_equal(model.wind().u(), model2.wind().u())
    np.testing.assert_array_almost_equal(model.wind().v(), model2.wind().v())

    pathlib.Path.unlink(pathlib.Path('test_wind.nc'))

def test_read_spectra():
    grid = dn.grid.Grid(lon=(0, 10), lat=(50,60))
    grid.set_spacing(nx=10, ny=10)

    model = dn.modelrun.Constant(grid, year=2020,month=1, day=1)
    model.import_spectra(point_picker=dn.pick.Trivial())
    model.spectra().ds().to_netcdf('test_spectra.nc')
    
    model2 = dn.modelrun.ModelRun(grid, year=2020,month=1, day=1)
    model2.import_spectra(dn.read.generic.PointNetcdf(), filename='test_spectra.nc',point_picker=dn.pick.Trivial())

    np.testing.assert_array_almost_equal(model.spectra().lon(), model2.spectra().lon())
    np.testing.assert_array_almost_equal(model.spectra().lat(), model2.spectra().lat())
    np.testing.assert_array_almost_equal(model.spectra().spec(), model2.spectra().spec())
    
    pathlib.Path.unlink(pathlib.Path('test_spectra.nc'))