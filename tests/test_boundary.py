from dnora.grid import Grid
from dnora import modelrun, spectra, pick
import numpy as np
from dnora.read.generic import ConstantData


def test_import_constant_boundary_one_point():
    grid = Grid(lon=5, lat=60)
    model = modelrun.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_spectra(ConstantData(), point_picker=pick.Trivial())

    assert model.spectra().size() == (25, 1, 10, 36)

    np.testing.assert_almost_equal(np.max(model.spectra().spec()), 1)
    np.testing.assert_almost_equal(np.min(model.spectra().spec()), 0)


def test_import_constant_boundary():
    grid = Grid(lon=(5, 6), lat=(60, 61))
    grid.set_spacing(nx=5, ny=10)

    model = modelrun.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_spectra(ConstantData(), point_picker=pick.Trivial())

    assert model.spectra().size() == (25, 50, 10, 36)

    np.testing.assert_almost_equal(np.max(model.spectra().spec()), 1)
    np.testing.assert_almost_equal(np.min(model.spectra().spec()), 0)
