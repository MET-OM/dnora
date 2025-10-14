from dnora.grid import Grid
from dnora import (
    modelrun,
    pick,
    export,
)
from dnora.type_manager.spectral_conventions import SpectralConvention
import numpy as np

from dnora.read.generic import ConstantData
from dnora.type_manager.dnora_types import DnoraDataType


def get_freq_and_dir_vector(math=False):
    f = np.linspace(0.1, 1, 10)  # np.loadtxt('data/freq.test')
    if math:
        D = np.mod(np.linspace(90.0, -260.0, 36), 360)
    else:
        D = np.linspace(0.0, 350.0, 36)
    return f, D


def find_north_from_dirs(dirs, north):
    return np.where(dirs == north)[0][0]


def find_north_from_spec(spec):
    spec = spec[0, 0, 0, :]
    return np.where(spec == 1)[0][0]


def test_import_export():
    grid = Grid(lon=(5, 6), lat=(60, 61))
    start_time = "2020-01-31 00:00:00"
    end_time = "2020-02-01 00:00:00"
    model = modelrun.ModelRun(grid=grid, start_time=start_time, end_time=end_time)

    model.import_wind(ConstantData())
    model.import_spectra(ConstantData(), point_picker=pick.Trivial())
    model.spectra_to_waveseries()
    model.import_waterlevel(ConstantData())

    for obj_type in [
        DnoraDataType.WIND,
        DnoraDataType.SPECTRA,
        DnoraDataType.SPECTRA1D,
        DnoraDataType.WAVESERIES,
        DnoraDataType.WATERLEVEL,
    ]:
        assert model[obj_type] is not None

    exporter = export.NullExporter(model)
    exporter.export_spectra()
    exporter.export_wind()
    exporter.export_spectra1d()
    exporter.export_waveseries()
    exporter.export_waterlevel()


def test_conventions():
    grid = Grid(lon=(5, 6), lat=(60, 61))
    start_time = "2020-01-31 00:00:00"
    end_time = "2020-02-01 00:00:00"
    model = modelrun.ModelRun(grid=grid, start_time=start_time, end_time=end_time)
    # Import constant spectra in oceanic convention with one component going north
    model.import_spectra(
        ConstantData(peaks={"freq": None}, convention=SpectralConvention.OCEAN),
        point_picker=pick.Trivial(),
    )

    assert model.spectra().convention() == SpectralConvention.OCEAN

    f, D = get_freq_and_dir_vector(math=False)

    np.testing.assert_array_almost_equal(model.spectra().dirs(), D)
    np.testing.assert_array_almost_equal(model.spectra().freq(), f)
    np.testing.assert_array_almost_equal(
        find_north_from_dirs(D, 0),
        find_north_from_spec(model.spectra().spec(dask=False)),
    )

    # Check 1D spectra convention
    model.spectra_to_1d()
    mdir = int(np.median(model.spectra1d().dirm(dask=False)))
    assert mdir == 180

    model.spectra_to_waveseries()  # Converts sepctra to MET before feeding into WaveSeries

    mdir = int(np.median(model.spectra1d().dirm(dask=False)))
    assert mdir == 180
    mdir = int(np.median(model.waveseries().dirm(dask=False)))
    assert mdir == 180

    # Meteorological convention
    model.spectra().set_convention(SpectralConvention.MET)
    assert model.spectra().convention() == SpectralConvention.MET
    np.testing.assert_array_almost_equal(model.spectra().dirs(), D)
    np.testing.assert_array_almost_equal(model.spectra().freq(), f)
    np.testing.assert_array_almost_equal(
        find_north_from_dirs(D, 180),
        find_north_from_spec(model.spectra().spec(dask=False)),
    )

    model.spectra_to_1d()
    mdir = int(np.median(model.spectra1d().dirm(dask=False)))
    assert mdir == 180
    model.spectra_to_waveseries()
    mdir = int(np.median(model.waveseries().dirm(dask=False)))
    assert mdir == 180

    # Mathematical convention (directional vector still starts from 0!)
    model.spectra().set_convention(SpectralConvention.MATH)
    assert model.spectra().convention() == SpectralConvention.MATH
    np.testing.assert_array_almost_equal(model.spectra().dirs(), D)
    np.testing.assert_array_almost_equal(model.spectra().freq(), f)
    np.testing.assert_array_almost_equal(
        find_north_from_dirs(D, 90),
        find_north_from_spec(model.spectra().spec(dask=False)),
    )

    model.spectra_to_1d()
    mdir = int(np.median(model.spectra1d().dirm(dask=False)))
    assert mdir == 180
    model.spectra_to_waveseries()
    mdir = int(np.median(model.waveseries().dirm(dask=False)))
    assert mdir == 180

    f, D = get_freq_and_dir_vector(math=True)
    # WW3 convention (Oceanic, but vector starts from 90 downwards)
    model.spectra().set_convention(SpectralConvention.WW3)
    assert model.spectra().convention() == SpectralConvention.WW3
    np.testing.assert_array_almost_equal(model.spectra().dirs(), D)
    np.testing.assert_array_almost_equal(model.spectra().freq(), f)
    np.testing.assert_array_almost_equal(
        find_north_from_dirs(D, 0),
        find_north_from_spec(model.spectra().spec(dask=False)),
    )

    model.spectra_to_1d()
    mdir = int(np.median(model.spectra1d().dirm(dask=False)))
    assert mdir == 180
    model.spectra_to_waveseries()
    mdir = int(np.median(model.waveseries().dirm(dask=False)))
    assert mdir == 180

    # MATHVEC convention (MAthematical and vector starts from 90 downwards)
    model.spectra().set_convention(SpectralConvention.MATHVEC)
    assert model.spectra().convention() == SpectralConvention.MATHVEC
    np.testing.assert_array_almost_equal(model.spectra().dirs(), D)
    np.testing.assert_array_almost_equal(model.spectra().freq(), f)
    np.testing.assert_array_almost_equal(
        find_north_from_dirs(D, 90),
        find_north_from_spec(model.spectra().spec(dask=False)),
    )

    model.spectra_to_1d()
    mdir = int(np.median(model.spectra1d().dirm(dask=False)))
    assert mdir == 180
    mdir = int(np.median(model.waveseries().dirm(dask=False)))
    assert mdir == 180
