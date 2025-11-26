import dnora as dn
import numpy as np


def test_import_constant_wavegrid_one_point():
    grid = dn.grid.Grid(lon=5, lat=60)
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_wavegrid(hs=4)

    assert model.wavegrid().size() == (25, 1, 1)

    np.testing.assert_almost_equal(np.mean(model.wavegrid().hs()), 4)


def test_import_constant_wavegrid():
    grid = dn.grid.Grid(lon=(5, 6), lat=(50, 60))
    grid.set_spacing(nx=5, ny=10)
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_wavegrid(hs=4)

    assert model.wavegrid().size() == (25, 10, 5)

    np.testing.assert_almost_equal(np.mean(model.wavegrid().hs()), 4)
