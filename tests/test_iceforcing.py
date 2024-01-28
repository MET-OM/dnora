from dnora import grd, mdl, ice
import numpy as np


def test_import_constant_ice_one_point():
    grid = grd.Grid(lon=5, lat=60)
    model = mdl.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_ice(
        ice.read.ConstantIce(concentration=0.4, thickness=0.5)
    )

    assert model.ice().size() == (25, 1, 1)

    np.testing.assert_almost_equal(np.mean(model.ice().concentration()), 0.4)
    np.testing.assert_almost_equal(np.mean(model.ice().thickness()), 0.5)


def test_import_constant_ice():
    grid = grd.Grid(lon=(5, 6), lat=(60, 61))
    grid.set_spacing(nx=5, ny=10)

    model = mdl.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_ice(
        ice.read.ConstantIce(concentration=0.4, thickness=0.5)
    )

    assert model.ice().size() == (25, 10, 5)

    np.testing.assert_almost_equal(np.mean(model.ice().concentration()), 0.4)
    np.testing.assert_almost_equal(np.mean(model.ice().thickness()), 0.5)
