import dnora as dn
import pytest
import pandas as pd
import numpy as np
from dnora.type_manager.data_sources import DataSource


@pytest.fixture(scope="session")
def grid():
    return dn.grid.Grid(lon=(10, 14), lat=(70, 71))


@pytest.fixture(scope="session")
def timevec():
    return pd.date_range("2022-01-01 00:00:00", "2022-01-02 23:00:00", freq="1h")


@pytest.fixture(scope="session")
def timevec_late2022():
    return pd.date_range("2022-10-01 00:00:00", "2022-10-02 23:00:00", freq="1h")


@pytest.fixture(scope="session")
def timevec2024():
    return pd.date_range("2024-01-05 00:00:00", "2024-01-05 23:00:00", freq="1h")


@pytest.mark.remote
def test_nora3(grid, timevec_late2022):
    model = dn.modelrun.ModelRun(
        grid, start_time="2022-10-01 00:00", end_time="2022-10-02 23:00"
    )
    model.import_ice(dn.read.ice.metno.NORA3())
    assert np.all(model.ice().time() == timevec_late2022)


@pytest.mark.remote
def test_barents25(grid, timevec2024):
    model = dn.modelrun.ModelRun(
        grid, start_time="2024-01-05 00:00", end_time="2024-01-05 23:00"
    )
    model.import_ice(dn.read.ice.metno.Barents25())
    assert np.all(model.ice().time() == timevec2024)


@pytest.mark.remote
def test_barents_patch(grid):
    """There are missing files (03 is missing the 18Z-folder) so we can use this to test patching"""
    model = dn.modelrun.ModelRun(
        grid, start_time="2024-01-03 00:00", end_time="2024-01-03 23:00"
    )
    model.import_ice(dn.read.ice.metno.Barents25())
    timevec = pd.date_range("2024-01-03 00:00:00", "2024-01-03 23:00:00", freq="1h")
    assert np.all(model.ice().time() == timevec)
