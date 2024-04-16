from dnora import grid, modelrun
import numpy as np
from dnora.readers.generic_readers import ConstantData


def test_import_constant_forcing_one_point():
    area = grid.Grid(lon=5, lat=60)
    model = modelrun.ModelRun(
        area, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_wind(ConstantData(), u=1.0, v=2.0)

    assert model.wind().size() == (25, 1, 1)

    np.testing.assert_almost_equal(np.mean(model.wind().u()), 1)
    np.testing.assert_almost_equal(np.mean(model.wind().v()), 2)
    np.testing.assert_almost_equal(
        np.mean(model.wind().magnitude()), (2**2 + 1**2) ** 0.5
    )


def test_import_constant_forcing():
    area = grid.Grid(lon=(5, 6), lat=(60, 61))
    area.set_spacing(ny=10, nx=5)
    model = modelrun.ModelRun(
        area, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_wind(ConstantData(), u=1.0, v=2.0)

    model.import_wind(ConstantData(), u=1.0, v=2.0)

    assert model.wind().size() == (25, 10, 5)

    np.testing.assert_almost_equal(np.mean(model.wind().u()), 1)
    np.testing.assert_almost_equal(np.mean(model.wind().v()), 2)
    np.testing.assert_almost_equal(
        np.mean(model.wind().magnitude()), (2**2 + 1**2) ** 0.5
    )
