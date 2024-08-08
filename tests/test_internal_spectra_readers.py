import dnora as dn
import pytest
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from dnora.dnora_type_manager.data_sources import DataSource


@pytest.fixture(scope="session")
def grid():
    return dn.grid.Grid(lon=6, lat=60)


@pytest.fixture(scope="session")
def timevec():
    return pd.date_range("2023-04-01 00:00:00", "2023-04-01 23:00:00", freq="1h")


@pytest.mark.internal
def test_wam3(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(dn.spectra.read_metno.WAM3(), source=DataSource.IMMUTABLE)
    assert np.all(model.spectra().time() == timevec)


@pytest.mark.internal
def test_wam800(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(dn.spectra.read_metno.WAM800(), source=DataSource.IMMUTABLE)
    assert np.all(model.spectra().time() == timevec)


@pytest.mark.internal
def test_nora3(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(dn.spectra.read_metno.NORA3(), source=DataSource.INTERNAL)
    assert np.all(model.spectra().time() == timevec)


@pytest.mark.internal
def test_ww3_4km(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(dn.spectra.read_metno.WW3_4km(), source=DataSource.IMMUTABLE)
    assert np.all(model.spectra().time() == timevec)