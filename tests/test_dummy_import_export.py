from dnora import grd, mdl, bnd, wnd, spc, wsr, pick, exp, wlv
from dnora.bnd.conventions import SpectralConvention
import pandas as pd
import numpy as np


def get_freq_and_dir_vector(math=False):
    f = np.linspace(0.1, 1, 10)  # np.loadtxt('data/freq.test')
    if math:
        D = np.mod(np.linspace(90.0, -255.0, 24), 360)
    else:
        D = np.linspace(0.0, 345.0, 24)
    return f, D


def find_north_from_dirs(dirs, north):
    return np.where(dirs == north)[0][0]


def find_north_from_spec(spec):
    spec = spec[0, 0, 0, :]
    return np.where(spec == 1)[0][0]


def test_import_export():
    grid = grd.Grid(lon=(5, 6), lat=(60, 61))
    start_time = "2020-01-31 00:00:00"
    end_time = "2020-02-01 00:00:00"
    model = mdl.ModelRun(grid=grid, start_time=start_time, end_time=end_time)

    model.import_forcing(wnd.read.ConstantForcing())
    model.import_boundary(
        bnd.read.ConstantBoundary(), point_picker=pick.TrivialPicker()
    )
    model.boundary_to_spectra()
    model.spectra_to_waveseries()
    model.import_waterlevel(wlv.read.ConstantWaterLevel())

    for obj_type in ["Forcing", "Boundary", "Spectra", "WaveSeries"]:
        assert model[obj_type] is not None

    exporter = exp.NullExporter(model)
    exporter.export_boundary()
    exporter.export_forcing()
    exporter.export_spectra()
    exporter.export_waveseries()
    exporter.export_waterlevel()
    exporter.write_input_file()


def test_conventions():
    grid = grd.Grid(lon=(5, 6), lat=(60, 61))
    start_time = "2020-01-31 00:00:00"
    end_time = "2020-02-01 00:00:00"
    model = mdl.ModelRun(grid=grid, start_time=start_time, end_time=end_time)
    # Import constant spectra in oceanic convention with one component going north
    model.import_boundary(
        bnd.read.ConstantBoundary(spectral_convention=SpectralConvention.OCEAN),
        point_picker=pick.TrivialPicker(),
    )

    assert model.boundary().convention() == SpectralConvention.OCEAN

    f, D = get_freq_and_dir_vector(math=False)

    np.testing.assert_array_almost_equal(model.boundary().dirs(), D)
    np.testing.assert_array_almost_equal(model.boundary().freq(), f)
    np.testing.assert_array_almost_equal(
        find_north_from_dirs(D, 0), find_north_from_spec(model.boundary().spec())
    )

    # Check 1D spectra convention
    model.boundary_to_spectra()
    assert model.spectra().convention() == SpectralConvention.OCEAN
    mdir = int(np.median(model.spectra().mdir()))
    assert mdir == 0
    model.spectra_to_waveseries()  # Converts sepctra to MET before feeding into WaveSeries
    mdir = int(np.median(model.waveseries().dirm()))
    assert mdir == 180

    # Meteorological convention
    model.boundary()._set_convention(SpectralConvention.MET)
    assert model.boundary().convention() == SpectralConvention.MET
    np.testing.assert_array_almost_equal(model.boundary().dirs(), D)
    np.testing.assert_array_almost_equal(model.boundary().freq(), f)
    np.testing.assert_array_almost_equal(
        find_north_from_dirs(D, 180), find_north_from_spec(model.boundary().spec())
    )

    model.boundary_to_spectra()
    assert model.spectra().convention() == SpectralConvention.MET
    mdir = int(np.median(model.spectra().mdir()))
    assert mdir == 180
    model.spectra_to_waveseries()
    mdir = int(np.median(model.waveseries().dirm()))
    assert mdir == 180

    # Mathematical convention (directional vector still starts from 0!)
    model.boundary()._set_convention(SpectralConvention.MATH)
    assert model.boundary().convention() == SpectralConvention.MATH
    np.testing.assert_array_almost_equal(model.boundary().dirs(), D)
    np.testing.assert_array_almost_equal(model.boundary().freq(), f)
    np.testing.assert_array_almost_equal(
        find_north_from_dirs(D, 90), find_north_from_spec(model.boundary().spec())
    )

    model.boundary_to_spectra()
    assert model.spectra().convention() == SpectralConvention.MATH
    mdir = int(np.median(model.spectra().mdir()))
    assert mdir == 90
    model.spectra_to_waveseries()
    mdir = int(np.median(model.waveseries().dirm()))
    assert mdir == 180

    f, D = get_freq_and_dir_vector(math=True)
    # WW3 convention (Oceanic, but vector starts from 90 downwards)
    model.boundary()._set_convention(SpectralConvention.WW3)
    assert model.boundary().convention() == SpectralConvention.WW3
    np.testing.assert_array_almost_equal(model.boundary().dirs(), D)
    np.testing.assert_array_almost_equal(model.boundary().freq(), f)
    np.testing.assert_array_almost_equal(
        find_north_from_dirs(D, 0), find_north_from_spec(model.boundary().spec())
    )

    model.boundary_to_spectra()
    assert model.spectra().convention() == SpectralConvention.OCEAN
    mdir = int(np.median(model.spectra().mdir()))
    assert mdir == 0
    model.spectra_to_waveseries()
    mdir = int(np.median(model.waveseries().dirm()))
    assert mdir == 180

    # MATHVEC convention (MAthematical and vector starts from 90 downwards)
    model.boundary()._set_convention(SpectralConvention.MATHVEC)
    assert model.boundary().convention() == SpectralConvention.MATHVEC
    np.testing.assert_array_almost_equal(model.boundary().dirs(), D)
    np.testing.assert_array_almost_equal(model.boundary().freq(), f)
    np.testing.assert_array_almost_equal(
        find_north_from_dirs(D, 90), find_north_from_spec(model.boundary().spec())
    )

    model.boundary_to_spectra()
    assert model.spectra().convention() == SpectralConvention.MATH
    mdir = int(np.median(model.spectra().mdir()))
    assert mdir == 90
    mdir = int(np.median(model.waveseries().dirm()))
    assert mdir == 180
