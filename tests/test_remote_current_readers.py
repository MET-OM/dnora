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
def test_norkyst800(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2022, month=4, day=1)
    model.import_current(dn.read.current.metno.NorKyst800(), program="pyfimex")

    assert np.all(model.current().time() == timevec)


@pytest.mark.remote
def test_norfjords160(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2022, month=4, day=1)
    model.import_current(dn.read.current.metno.NorFjords160(), program="pyfimex")

    assert np.all(model.current().time() == timevec)
