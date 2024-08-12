import dnora as dn
import pytest
import pandas as pd
import numpy as np
from dnora.dnora_type_manager.data_sources import DataSource


@pytest.fixture(scope="session")
def grid():
    return dn.grid.Grid(lon=(10, 14), lat=(70, 71))


@pytest.fixture(scope="session")
def timevec():
    return pd.date_range("2022-01-01 00:00:00", "2022-01-02 23:00:00", freq="1h")


@pytest.mark.remote
def test_nora3(grid, timevec):
    model = dn.modelrun.ModelRun(
        grid, start_time="2022-01-01 00:00", end_time="2022-01-02 23:00"
    )
    model.import_ice(dn.ice.read_metno.NORA3())
    assert np.all(model.ice().time() == timevec)
