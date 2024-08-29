from dnora import modelrun, pick
from dnora.read.generic import ConstantData
import numpy as np
import dnora as dn
import pandas as pd


def test_import_constant_waveseries_one_point():
    grid = dn.grid.Grid(lon=5, lat=60)
    model = modelrun.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_spectra(ConstantData(), point_picker=pick.Trivial())
    model.spectra_to_waveseries()

    assert model.waveseries().size() == (25, 1)
    manual_hs = 4 * np.sqrt(
        np.sum(model.spectra().isel(time=1).spec())
        * model.spectra().dd(angular=True)
        * model.spectra().df()
    )
    np.testing.assert_almost_equal(np.mean(model.waveseries().hs()), manual_hs)


def test_import_wave_from_norac_spectra():
    grid = dn.grid.Grid(lon=9.834990, lat=63.571444)
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=10)
    model.import_waveseries(
        dn.spectra.read_metno.NORAC(),
    )
    timevec = pd.date_range("2023-04-10 00:00:00", "2023-04-10 23:00:00", freq="1h")
    assert np.all(model.waveseries().time() == timevec)
    assert model.waveseries().wind() is not None
    assert model.waveseries().wind_dir() is not None
    assert model.waveseries().current() is not None
    assert model.waveseries().current_dir() is not None
    assert model.waveseries().depth() is not None


def test_import_wave_from_ww3_4km_spectra():
    grid = dn.grid.Grid(lon=9.834990, lat=63.571444)
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=10)
    model.import_waveseries(
        dn.spectra.read_metno.WW3_4km(),
    )
    timevec = pd.date_range("2023-04-10 00:00:00", "2023-04-10 23:00:00", freq="1h")
    assert np.all(model.waveseries().time() == timevec)
    assert model.waveseries().wind() is not None
    assert model.waveseries().wind_dir() is not None
    assert model.waveseries().hs() is not None
    assert model.waveseries().tp() is not None
    assert model.waveseries().dirp() is not None


def test_import_wave_from_wam800_spectra():
    grid = dn.grid.Grid(lon=9.834990, lat=63.571444)
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=10)
    model.import_waveseries(
        dn.spectra.read_metno.WAM800(),
    )
    timevec = pd.date_range("2023-04-10 00:00:00", "2023-04-10 23:00:00", freq="1h")
    assert np.all(model.waveseries().time() == timevec)
    assert model.waveseries().wind() is not None
    assert model.waveseries().wind_dir() is not None
    assert model.waveseries().hs() is not None
    assert model.waveseries().tp() is not None
    assert model.waveseries().dirp() is not None
    assert model.waveseries().dirm_swell() is not None
    assert model.waveseries().dirm_sea() is not None
    assert model.waveseries().depth() is not None


def test_import_wave_from_wam4km_spectra():
    grid = dn.grid.Grid(lon=6, lat=60)
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=10)
    model.import_waveseries(
        dn.spectra.read_metno.WAM4km(),
    )
    timevec = pd.date_range("2023-04-10 00:00:00", "2023-04-10 23:00:00", freq="1h")
    assert np.all(model.waveseries().time() == timevec)
    assert model.waveseries().wind() is not None
    assert model.waveseries().wind_dir() is not None
    assert model.waveseries().hs() is not None
    assert model.waveseries().tp() is not None
    assert model.waveseries().dirp() is not None
    assert model.waveseries().dirm_swell() is not None
    assert model.waveseries().dirm_sea() is not None
    assert model.waveseries().depth() is not None


def test_norac():
    pass
