from dnora.grid import Grid
from dnora import modelrun
from dnora.read.generic import ConstantData
import numpy as np


def test_import_constant_waterlevel_one_point():
    grid = Grid(lon=5, lat=60)
    model = modelrun.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_waterlevel(ConstantData(), eta=1.0)

    assert model.waterlevel().size() == (25, 1, 1)

    np.testing.assert_almost_equal(np.mean(model.waterlevel().eta()), 1)


def test_import_constant_waterlevel():
    grid = Grid(lon=(5, 6), lat=(60, 61))
    grid.set_spacing(nx=5, ny=10)

    model = modelrun.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_waterlevel(ConstantData(), eta=1.0)

    assert model.waterlevel().size() == (25, 10, 5)

    np.testing.assert_almost_equal(np.mean(model.waterlevel().eta()), 1)
