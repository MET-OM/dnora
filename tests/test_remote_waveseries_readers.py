import dnora as dn
import pytest
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from dnora.type_manager.data_sources import DataSource
import geo_parameters as gp


@pytest.mark.remote
def test_e39():
    model = dn.modelrun.ModelRun(
        start_time="2018-01-01 00:00", end_time="2018-05-01 00:00"
    )
    model.import_waveseries(dn.read.waveseries.metno.E39(), tile="D")
    assert model.waveseries().time()[0] == pd.Timestamp("2018-01-01 00:00")
    assert model.waveseries().time()[-1] == pd.Timestamp("2018-05-01 00:00")


@pytest.mark.remote
def test_import_wave_from_norac_spectra():
    grid = dn.grid.Grid(lon=9.834990, lat=63.571444)
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=10)
    model.import_waveseries(
        dn.read.spectra.metno.NORAC(),
    )
    timevec = pd.date_range("2023-04-10 00:00:00", "2023-04-10 23:00:00", freq="1h")
    assert np.all(model.waveseries().time() == timevec)
    assert model.waveseries().ff() is not None
    assert model.waveseries().dd() is not None
    assert model.waveseries().current() is not None
    assert model.waveseries().current_dir() is not None
    assert model.waveseries().depth() is not None


@pytest.mark.remote
def test_import_wave_from_ww3_4km_spectra():
    grid = dn.grid.Grid(lon=9.834990, lat=63.571444)
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=10)
    model.import_waveseries(
        dn.read.spectra.metno.WW3_4km(),
    )
    timevec = pd.date_range("2023-04-10 00:00:00", "2023-04-10 23:00:00", freq="1h")
    assert np.all(model.waveseries().time() == timevec)
    assert model.waveseries().ff() is not None
    assert model.waveseries().dd() is not None
    assert model.waveseries().hs() is not None
    assert model.waveseries().tp() is not None
    assert model.waveseries().dirp() is not None


@pytest.mark.remote
def test_import_wave_from_wam800_spectra():
    grid = dn.grid.Grid(lon=9.834990, lat=63.571444)
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=10)
    model.import_waveseries(
        dn.read.spectra.metno.WAM800(),
    )
    timevec = pd.date_range("2023-04-10 00:00:00", "2023-04-10 23:00:00", freq="1h")
    assert np.all(model.waveseries().time() == timevec)
    assert model.waveseries().ff() is not None
    assert model.waveseries().dd() is not None
    assert model.waveseries().hs() is not None
    assert model.waveseries().tp() is not None
    assert model.waveseries().dirp() is not None
    # assert model.waveseries().dirm_swell() is not None
    # assert model.waveseries().dirm_sea() is not None
    assert model.waveseries().depth() is not None


@pytest.mark.remote
def test_import_wave_from_wam4km_spectra():
    grid = dn.grid.Grid(lon=6, lat=60)
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=10)
    model.import_waveseries(
        dn.read.spectra.metno.WAM4km(),
    )
    timevec = pd.date_range("2023-04-10 00:00:00", "2023-04-10 23:00:00", freq="1h")
    assert np.all(model.waveseries().time() == timevec)
    assert model.waveseries().ff() is not None
    assert model.waveseries().dd() is not None
    assert model.waveseries().hs() is not None
    assert model.waveseries().tp() is not None
    assert model.waveseries().dirp() is not None
    # assert model.waveseries().dirm_swell() is not None
    # assert model.waveseries().dirm_sea() is not None
    assert model.waveseries().depth() is not None


@pytest.mark.remote
def test_norac_wave():
    grid = dn.grid.Grid(lon=9.834990, lat=63.571444)
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=10)
    model.import_waveseries(dn.read.waveseries.metno.NORAC())
    timevec = pd.date_range("2023-04-10 00:00:00", "2023-04-10 23:00:00", freq="1h")
    assert np.all(model.waveseries().time() == timevec)
    assert model.waveseries().hs() is not None
