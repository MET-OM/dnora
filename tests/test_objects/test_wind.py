import dnora as dn
import numpy as np

from dnora.read.generic import ConstantData


def test_import_constant_wind_one_point():
    grid = dn.grid.Grid(lon=5, lat=60)
    model = dn.modelrun.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_wind(ConstantData(), u=1, v=2)

    assert model.wind().size() == (25, 1, 1)

    np.testing.assert_almost_equal(np.mean(model.wind().u()), 1)
    np.testing.assert_almost_equal(np.mean(model.wind().v()), 2)
    np.testing.assert_almost_equal(np.mean(model.wind().mag()), (2**2 + 1**2) ** 0.5)

    assert "from" in model.wind().meta.get("dir").get("standard_name")
    np.testing.assert_almost_equal(
        np.mean(model.wind().dir()),
        90 - np.rad2deg(np.arctan2(2, 1)) + 180,
    )


def test_import_constant_wind():
    grid = dn.grid.Grid(lon=(5, 6), lat=(60, 61))
    grid.set_spacing(nx=5, ny=10)

    model = dn.modelrun.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_wind(ConstantData(), u=1, v=2)

    assert model.wind().size() == (25, 10, 5)

    np.testing.assert_almost_equal(np.mean(model.wind().u()), 1)
    np.testing.assert_almost_equal(np.mean(model.wind().v()), 2)
    np.testing.assert_almost_equal(np.mean(model.wind().mag()), (2**2 + 1**2) ** 0.5)

    np.testing.assert_almost_equal(
        np.mean(model.wind().dir()),
        90 - np.rad2deg(np.arctan2(2, 1)) + 180,
    )
