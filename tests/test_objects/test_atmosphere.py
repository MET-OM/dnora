import dnora as dn
import numpy as np


def test_import_constant_ocean_one_point():
    grid = dn.grid.Grid(lon=5, lat=60)
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_atmosphere(t2m=10, r=30)

    assert model.atmosphere().size() == (25, 1, 1)

    np.testing.assert_almost_equal(np.mean(model.atmosphere().t2m()), 10)
    np.testing.assert_almost_equal(np.mean(model.atmosphere().r()), 30)


def test_import_constant_ocean():
    grid = dn.grid.Grid(lon=(5, 6), lat=(50, 60))
    grid.set_spacing(nx=5, ny=10)
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_atmosphere(t2m=10, r=30)

    assert model.atmosphere().size() == (25, 10, 5)
    np.testing.assert_almost_equal(np.mean(model.atmosphere().t2m()), 10)
    np.testing.assert_almost_equal(np.mean(model.atmosphere().r()), 30)
