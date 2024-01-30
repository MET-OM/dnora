from dnora.grid import Grid
from dnora import modelrun, pick
import numpy as np

from dnora.readers.generic_readers import ConstantPointData


def test_import_constant_spectra_one_point():
    grid = Grid(lon=5, lat=60)
    model = modelrun.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_spectra(ConstantPointData(), point_picker=pick.TrivialPicker())
    model.spectra_to_1d()

    assert model.spectra1d().size() == (25, 1, 10)

    np.testing.assert_almost_equal(
        np.max(model.spectra1d().spec()),
        np.max(model.spectra().spec()) * (model.spectra().dd()) * np.pi / 180,
    )


def test_import_constant_spectra():
    grid = Grid(lon=(5, 6), lat=(60, 61))
    grid.set_spacing(nx=5, ny=10)

    model = modelrun.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_spectra(ConstantPointData(), point_picker=pick.TrivialPicker())
    model.spectra_to_1d()
    assert model.spectra1d().size() == (25, 50, 10)

    np.testing.assert_almost_equal(
        np.max(model.spectra1d().spec()),
        np.max(model.spectra().spec()) * (model.spectra().dd()) * np.pi / 180,
    )
