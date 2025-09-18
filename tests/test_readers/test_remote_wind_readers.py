import dnora as dn
import pytest
import pandas as pd
import numpy as np
from dnora.read.wind.metno import get_meps_urls, MEPS
from dnora.type_manager.data_sources import DataSource


def grid_is_covered(grid, skeleton):
    assert skeleton.lon()[0] < grid.lon()[0]
    assert skeleton.lon()[-1] >= grid.lon()[-1]
    assert skeleton.lat()[0] < grid.lat()[0]
    assert skeleton.lat()[-1] >= grid.lat()[-1]


@pytest.fixture(scope="session")
def grid():
    return dn.grid.Grid(lon=(10, 14), lat=(60, 61))


@pytest.fixture(scope="session")
def timevec():
    return pd.date_range("2022-04-01 00:00:00", "2022-04-01 23:00:00", freq="1h")


@pytest.fixture(scope="session")
def timevec202508():
    return pd.date_range("2025-08-01 00:00:00", "2025-08-01 23:00:00", freq="1h")


@pytest.fixture(scope="session")
def timevec202505():
    return pd.date_range("2025-05-01 00:00:00", "2025-05-01 23:00:00", freq="1h")


@pytest.mark.remote
def test_nora3(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2022, month=4, day=1)
    model.import_wind(dn.read.wind.metno.NORA3(), program="pyfimex")
    assert np.all(model.wind().time() == timevec)


@pytest.mark.remote
def test_nora3_2025(grid, timevec202505):
    model = dn.modelrun.ModelRun(grid, year=2025, month=5, day=1)
    model.import_wind(dn.read.wind.metno.NORA3(), program="pyfimex")
    assert np.all(model.wind().time() == timevec202505)


@pytest.mark.remote
def test_mywave3km(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2022, month=4, day=1)
    model.import_wind(dn.read.wind.metno.MyWave3km())
    assert np.all(model.wind().time() == timevec)


@pytest.mark.remote
def test_meps(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2022, month=4, day=1)
    model.import_wind(dn.read.wind.metno.MEPS())
    assert np.all(model.wind().time() == timevec)


@pytest.mark.remote
def test_meps_subset(grid):
    model = dn.modelrun.ModelRun(grid, year=2020, month=2, day=3)
    model.import_wind(dn.read.wind.metno.MEPS())
    timevec = pd.date_range("2020-02-03 00:00:00", "2020-02-03 23:00:00", freq="1h")
    assert np.all(model.wind().time() == timevec)


def test_meps_url():
    file_times = pd.date_range("2020-02-04 00:00:00", "2020-02-04 23:00:00", freq="6h")
    urls = get_meps_urls(
        MEPS.product_configuration.default_folders[DataSource.REMOTE],
        MEPS.product_configuration.filename,
        file_times,
    )
    assert len(urls) == 4
    assert "subset" in urls[0]
    assert "subset" in urls[1]
    assert "det" in urls[2]
    assert "det" in urls[3]


@pytest.mark.remote
def test_meps_det_subset2025(grid, timevec202508):
    model = dn.modelrun.ModelRun(grid, year=2025, month=8, day=1)
    model.import_wind(dn.read.wind.metno.MEPS())
    assert np.all(model.wind().time() == timevec202508)


@pytest.mark.remote
def test_meps_det_subset(grid):
    model = dn.modelrun.ModelRun(grid, year=2020, month=2, day=4)
    model.import_wind(dn.read.wind.metno.MEPS())
    timevec = pd.date_range("2020-02-04 00:00:00", "2020-02-04 23:00:00", freq="1h")
    assert np.all(model.wind().time() == timevec)


@pytest.mark.remote
def test_meps_old_archive(grid):
    start_time = "2019-12-31 03:00:00"
    end_time = "2020-01-01 03:00:00"
    model = dn.modelrun.ModelRun(grid, start_time=start_time, end_time=end_time)
    model.import_wind(dn.read.wind.metno.MEPS())
    timevec = pd.date_range(start_time, end_time, freq="1h")
    assert np.all(model.wind().time() == timevec)


@pytest.mark.remote
def test_nora3fp(grid, timevec):
    """NORA3 reader reading the original hourly files"""
    model = dn.modelrun.ModelRun(grid, year=2022, month=4, day=1)
    model.import_wind(dn.read.wind.metno.NORA3_fp())
    assert np.all(model.wind().time() == timevec)


@pytest.mark.slow
@pytest.mark.remote
def test_era5(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2022, month=4, day=1)
    model.import_wind(dn.read.wind.ec.ERA5())

    assert np.all(model.wind().time() == timevec)

    grid_is_covered(grid, model.wind())
