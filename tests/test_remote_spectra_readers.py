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


@pytest.mark.remote
def test_nora3(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(dn.spectra.read_metno.NORA3())
    assert np.all(model.spectra().time() == timevec)


@pytest.mark.remote
def test_nora3_generic(timevec):
    grid = dn.grid.Grid(lon=9.834990, lat=63.571444)
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(
        dn.spectra.read.WAM(stride=24, hours_per_file=24),
        folder="https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/%Y/%m",
        filename="SPC%Y%m%d00.nc",
    )
    assert np.all(model.spectra().time() == timevec)


@pytest.mark.remote
def test_wam4km(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(dn.spectra.read_metno.WAM4km())
    assert np.all(model.spectra().time() == timevec)


@pytest.mark.remote
def test_wam4km_generic(timevec):
    grid = dn.grid.Grid(lon=9.834990, lat=63.571444)
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(
        dn.spectra.read.WAM(stride=6, hours_per_file=73),
        source=DataSource.REMOTE,
        folder="https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/%Y/%m/%d",
        filename="MyWave_wam4_SPC_%Y%m%dT%HZ.nc",
    )
    assert np.all(model.spectra().time() == timevec)


@pytest.mark.remote
def test_ww3_4km(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(dn.spectra.read_metno.WW3_4km())
    assert np.all(model.spectra().time() == timevec)


@pytest.mark.remote
def test_ww3_4km_generic(timevec):
    grid = dn.grid.Grid(lon=9.834990, lat=63.571444)
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(
        dn.spectra.read.WAM(stride=6, hours_per_file=73),
        source=DataSource.REMOTE,
        folder="https://thredds.met.no/thredds/dodsC/ww3_4km_archive_files/%Y/%m/%d",
        filename="ww3_4km_POI_SPC_%Y%m%dT%HZ.nc",
    )
    assert np.all(model.spectra().time() == timevec)


@pytest.mark.remote
def test_norac(timevec):
    grid = dn.grid.Grid(lon=9.834990, lat=63.571444)
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(dn.spectra.read_metno.NORAC())
    assert np.all(model.spectra().time() == timevec)


@pytest.mark.remote
def test_norac_generic(timevec):
    grid = dn.grid.Grid(lon=9.834990, lat=63.571444)
    model = dn.modelrun.ModelRun(grid, year=2023, month=4, day=1)
    model.import_spectra(
        dn.spectra.read.WW3(),
        source=DataSource.REMOTE,
        folder="https://thredds.met.no/thredds/dodsC/norac_wave/spec/",
    )
    assert np.all(model.spectra().time() == timevec)


# For some reason the lon/lat is not the same as in the immutable directory and we can't find any points?
# @pytest.mark.remote
# def test_wam3(grid, timevec):
#     today = datetime.today() - timedelta(days=2)
#     model = dn.modelrun.ModelRun(
#         grid, year=today.year, month=today.month, day=today.day
#     )
#     model.import_spectra(dn.spectra.read_metno.WAM3(), source="remote")
#     assert np.all(model.spectra().time() == timevec)
