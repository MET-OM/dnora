from dnora import grd, mdl, ocr
import numpy as np


def test_import_constant_oceancurrent_one_point():
    grid = grd.Grid(lon=5, lat=60)
    model = mdl.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_oceancurrent(ocr.read.ConstantOceanCurrent(u=1, v=2))

    assert model.oceancurrent().size() == (25, 1, 1)

    np.testing.assert_almost_equal(np.mean(model.oceancurrent().u()), 1)
    np.testing.assert_almost_equal(np.mean(model.oceancurrent().v()), 2)
    np.testing.assert_almost_equal(
        np.mean(model.oceancurrent().magnitude()), (2**2 + 1**2) ** 0.5
    )

    np.testing.assert_almost_equal(
        np.mean(model.oceancurrent().direction()),
        90 - np.rad2deg(np.arctan2(2, 1)) + 180,
    )


def test_import_constant_waterlevel():
    grid = grd.Grid(lon=(5, 6), lat=(60, 61))
    grid.set_spacing(nx=5, ny=10)

    model = mdl.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_oceancurrent(ocr.read.ConstantOceanCurrent(u=1, v=2))

    assert model.oceancurrent().size() == (25, 10, 5)

    np.testing.assert_almost_equal(np.mean(model.oceancurrent().u()), 1)
    np.testing.assert_almost_equal(np.mean(model.oceancurrent().v()), 2)
    np.testing.assert_almost_equal(
        np.mean(model.oceancurrent().magnitude()), (2**2 + 1**2) ** 0.5
    )
