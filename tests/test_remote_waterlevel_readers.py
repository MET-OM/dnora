import dnora as dn
import pytest
import pandas as pd
import numpy as np
from dnora.read.wind.metno import get_meps_urls, MEPS
from dnora.type_manager.data_sources import DataSource


def grid_is_covered(grid, skeleton):
    assert skeleton.lon()[0] < grid.lon()[0]
    assert skeleton.lon()[-1] >= grid.lon()[-1]
    assert skeleton.lat()[0] < grid.lat()[0]
    assert skeleton.lat()[-1] >= grid.lat()[-1]


@pytest.fixture(scope="session")
def grid():
    return dn.grid.Grid(lon=(10, 14), lat=(60, 61))


@pytest.fixture(scope="session")
def timevec():
    return pd.date_range("2022-04-01 00:00:00", "2022-04-01 23:00:00", freq="1h")


@pytest.mark.slow
@pytest.mark.remote
def test_era5(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2022, month=4, day=1)
    model.import_waterlevel(dn.read.waterlevel.ec.GTSM_ERA5())

    assert np.all(model.waterlevel().time() == timevec)

    grid_is_covered(grid, model.waterlevel())
