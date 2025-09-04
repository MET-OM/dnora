import dnora as dn
import pytest
import pandas as pd
import numpy as np
from dnora.read.wind.metno import get_meps_urls, MEPS
from dnora.type_manager.data_sources import DataSource
import os, shutil

def grid_is_covered(grid, skeleton):
    assert skeleton.lon()[0] < grid.lon()[0]
    assert skeleton.lon()[-1] >= grid.lon()[-1]
    assert skeleton.lat()[0] < grid.lat()[0]
    assert skeleton.lat()[-1] >= grid.lat()[-1]


@pytest.fixture(scope="session")
def grid():
    return dn.grid.Grid(lon=(4, 6), lat=(58, 61))


@pytest.fixture(scope="session")
def timevec():
    return pd.date_range("2023-04-01 00:00:00", "2023-04-01 23:00:00", freq="1h")
def cleanup():
    if os.path.isdir("dnora_waterlevel_temp"):
        shutil.rmtree("dnora_waterlevel_temp")


@pytest.mark.slow
@pytest.mark.remote
def test_era5(grid, timevec):
    cleanup()
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_waterlevel(dn.read.waterlevel.ec.GTSM_ERA5())

    assert np.all(model.waterlevel().time() == timevec)

    grid_is_covered(grid, model.waterlevel())
    cleanup()


@pytest.mark.remote
def test_cmems(grid, timevec):
    cleanup()
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_waterlevel(dn.read.waterlevel.cmems.Global())
    assert np.all(model.waterlevel().time() == timevec)

    grid_is_covered(grid, model.waterlevel())
    cleanup()