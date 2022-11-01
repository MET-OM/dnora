from dnora import grd
from dnora import mdl
from dnora import bnd
from dnora import wnd
from dnora import spc
from dnora import wsr
from dnora.bnd.conventions import SpectralConvention
import pandas as pd
import numpy as np

def get_freq_and_dir_vector(math=False):
	f = np.linspace(0.1,1,10) #np.loadtxt('data/freq.test')
	if math:
		D = np.mod(np.linspace(90.,-255.,24),360)
	else:
		D = np.linspace(0.,345.,24)
	return f, D

def find_north_from_dirs(dirs, north):
    return np.where(dirs==north)[0][0]

def find_north_from_spec(spec):
    spec = spec[0,0,0,:]
    return np.where(spec==1)[0][0]

def test_import_export():
    grid = grd.Grid(lon=(5,6), lat=(60,61))
    start_time = '2020-01-31 00:00:00'
    end_time = '2020-02-01 00:00:00'
    model = mdl.ModelRun(grid=grid, start_time=start_time, end_time=end_time)

    model.import_forcing(wnd.read.ConstantForcing())
    model.import_boundary(bnd.read.ConstantBoundary(grid), point_picker=bnd.pick.TrivialPicker())
    model.boundary_to_spectra()
    model.spectra_to_waveseries()
    # model.export_grid(grd.write.Null())
    # model.export_forcing(wnd.write.Null())
    # model.export_boundary(bnd.write.Null())
    # model.export_spectra(spc.write.Null())
    # model.export_waveseries(wsr.write.Null())



def test_boundary_processing():
    grid = grd.Grid(lon=(5,6), lat=(60,61))
    start_time = '2020-01-31 00:00:00'
    end_time = '2020-02-01 00:00:00'
    model = mdl.ModelRun(grid=grid, start_time=start_time, end_time=end_time)
    model.import_boundary(bnd.read.ConstantBoundary(grid, spectral_convention=SpectralConvention.OCEAN), point_picker=bnd.pick.TrivialPicker())
    model.boundary_to_spectra()
    model.spectra_to_waveseries()

    f, D = get_freq_and_dir_vector(math=False)

    np.testing.assert_array_almost_equal(model.boundary().dirs(),D)
    np.testing.assert_array_almost_equal(model.boundary().freq(),f)
    np.testing.assert_array_almost_equal(find_north_from_dirs(D,0), find_north_from_spec(model.boundary().spec()))

    model.boundary()._set_convention(SpectralConvention.MET)
    np.testing.assert_array_almost_equal(model.boundary().dirs(),D)
    np.testing.assert_array_almost_equal(model.boundary().freq(),f)
    np.testing.assert_array_almost_equal(find_north_from_dirs(D,180), find_north_from_spec(model.boundary().spec()))

    model.boundary()._set_convention(SpectralConvention.MATH)
    np.testing.assert_array_almost_equal(model.boundary().dirs(),D)
    np.testing.assert_array_almost_equal(model.boundary().freq(),f)
    np.testing.assert_array_almost_equal(find_north_from_dirs(D,90), find_north_from_spec(model.boundary().spec()))

    f, D = get_freq_and_dir_vector(math=True)

    model.boundary()._set_convention(SpectralConvention.WW3)
    np.testing.assert_array_almost_equal(model.boundary().dirs(),D)
    np.testing.assert_array_almost_equal(model.boundary().freq(),f)
    np.testing.assert_array_almost_equal(find_north_from_dirs(D,0), find_north_from_spec(model.boundary().spec()))

    model.boundary()._set_convention(SpectralConvention.MATHVEC)
    np.testing.assert_array_almost_equal(model.boundary().dirs(),D)
    np.testing.assert_array_almost_equal(model.boundary().freq(),f)
    np.testing.assert_array_almost_equal(find_north_from_dirs(D,90), find_north_from_spec(model.boundary().spec()))
