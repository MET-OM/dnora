import dnora as dn
import pytest
import pandas as pd
import numpy as np


@pytest.fixture(scope="session")
def grid():
    return dn.grid.Grid(lon=(10, 14), lat=(60, 61))


@pytest.fixture(scope="session")
def timevec():
    return pd.date_range("2022-04-01 00:00:00", "2022-04-01 23:00:00", freq="1h")


@pytest.mark.remote
def test_nora3(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2022, month=4, day=1)
    model.import_wind(dn.wind.read_metno.NORA3(), program="fimex")
    assert np.all(model.wind().time() == timevec)


@pytest.mark.remote
def test_mywave3km(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2022, month=4, day=1)
    model.import_wind(dn.wind.read_metno.MyWave3km())
    assert np.all(model.wind().time() == timevec)


@pytest.mark.remote
def test_meps(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2022, month=4, day=1)
    model.import_wind(dn.wind.read_metno.MEPS())
    assert np.all(model.wind().time() == timevec)


@pytest.mark.remote
def test_nora3fp(grid, timevec):
    """NORA3 reader reading the original hourly files"""
    model = dn.modelrun.ModelRun(grid, year=2022, month=4, day=1)
    model.import_wind(dn.wind.read_metno.NORA3_fp())
    assert np.all(model.wind().time() == timevec)


@pytest.mark.remote
def test_era5(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2022, month=4, day=1)
    model.import_wind(dn.wind.read_ec.ERA5())
    assert np.all(model.wind().time() == timevec)
