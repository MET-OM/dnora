import dnora as dn
import pytest
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from dnora.type_manager.data_sources import DataSource


@pytest.fixture(scope="session")
def grid():
    return dn.grid.Grid(lon=6, lat=60)


@pytest.fixture(scope="session")
def timevec():
    return pd.date_range("2023-04-01 00:00:00", "2023-04-01 23:00:00", freq="1h")


@pytest.mark.remote
def test_nora3(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(dn.read.spectra.metno.NORA3())
    assert np.all(model.spectra().time() == timevec)
    assert model.spectra().spec(strict=True) is not None

@pytest.mark.remote
def test_wam4km(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(dn.read.spectra.metno.WAM4km())
    assert np.all(model.spectra().time() == timevec)
    assert model.spectra().spec(strict=True) is not None

@pytest.mark.remote
def test_ww3_4km(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(dn.read.spectra.metno.WW3_4km())
    assert np.all(model.spectra().time() == timevec)
    assert model.spectra().spec(strict=True) is not None

@pytest.mark.remote
def test_norac(timevec):
    grid = dn.grid.Grid(lon=9.834990, lat=63.571444)
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(dn.read.spectra.metno.NORAC())
    assert np.all(model.spectra().time() == timevec)
    assert model.spectra().spec(strict=True) is not None
