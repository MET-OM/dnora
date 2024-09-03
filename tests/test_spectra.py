from dnora.grid import Grid
from dnora import modelrun, pick
import numpy as np

from dnora.read.generic import ConstantData
from dnora.type_manager.spectral_conventions import SpectralConvention


def test_import_constant_spectra_one_point():
    grid = Grid(lon=5, lat=60)
    model = modelrun.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_spectra(ConstantData(), point_picker=pick.Trivial())
    model.spectra_to_1d()

    assert model.spectra1d().size() == (25, 1, 10)
    np.testing.assert_array_almost_equal(
        model.spectra1d().spec(time="2020-01-01 00:00", dask=True),
        np.sum(
            model.spectra().spec(time="2020-01-01 00:00", dask=True)
            * model.spectra().dd(angular=True),
            axis=1,
        ),
    )

    assert model.spectra().convention() == SpectralConvention.OCEAN


def test_import_constant_spectra_on_gridded():
    grid = Grid(lon=(5, 6), lat=(60, 61))
    grid.set_spacing(nx=5, ny=10)

    model = modelrun.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_spectra(ConstantData(), point_picker=pick.Trivial())
    model.spectra_to_1d()
    assert model.spectra1d().size() == (25, 50, 10)

    np.testing.assert_array_almost_equal(
        model.spectra1d().spec(time="2020-01-01 00:00"),
        np.sum(
            model.spectra().spec(time="2020-01-01 00:00")
            * model.spectra().dd(angular=True),
            axis=2,
        ),
    )

    assert model.spectra().convention() == SpectralConvention.OCEAN
