import dnora as dn
import pytest
import pandas as pd
import numpy as np
import os, shutil


def grid_is_covered(grid, skeleton):
    assert skeleton.lon()[0] < grid.lon()[0]
    assert skeleton.lon()[-1] >= grid.lon()[-1]
    assert skeleton.lat()[0] < grid.lat()[0]
    assert skeleton.lat()[-1] >= grid.lat()[-1]
@pytest.fixture(scope="session")
def grid():
    return dn.grid.Grid(lon=(10, 14), lat=(60, 61))

@pytest.fixture(scope="session")
def grid2():
    return dn.grid.Grid(lon=(4, 6), lat=(58, 61))


@pytest.fixture(scope="session")
def grid160():
    return dn.grid.Grid(lon=(4.8, 5.2), lat=(59.9, 60.1))


@pytest.fixture(scope="session")
def timevec():
    return pd.date_range("2022-04-01 00:00:00", "2022-04-01 23:00:00", freq="1h")


@pytest.fixture(scope="session")
def timevec2017():
    return pd.date_range("2017-04-01 00:00:00", "2017-04-01 23:00:00", freq="1h")

def cleanup():
    if os.path.isdir("dnora_current_temp"):
        shutil.rmtree("dnora_current_temp")
@pytest.mark.remote
def test_norkyst800(grid, timevec):
    model = dn.modelrun.ModelRun(grid, year=2022, month=4, day=1)
    model.import_current(dn.read.current.metno.NorKyst800(), program="pyfimex")

    assert np.all(model.current().time() == timevec)
    assert model.current().u(strict=True) is not None
    assert model.current().v(strict=True) is not None


@pytest.mark.remote
def test_norkyst800_2017(grid, timevec2017):
    model = dn.modelrun.ModelRun(grid, year=2017, month=4, day=1)
    model.import_current(dn.read.current.metno.NorKyst800(), program="pyfimex")

    assert np.all(model.current().time() == timevec2017)
    assert model.current().u(strict=True) is not None
    assert model.current().v(strict=True) is not None

@pytest.mark.internal
@pytest.mark.remote
def test_norfjords160(grid160, timevec):
    model = dn.modelrun.ModelRun(grid160, year=2022, month=4, day=1)
    model.import_current(dn.read.current.metno.NorFjords160(), program="pyfimex")

    assert np.all(model.current().time() == timevec)
    assert model.current().u(strict=True) is not None
    assert model.current().v(strict=True) is not None

@pytest.mark.remote
def test_cmems_global(grid, timevec):
    cleanup()
    model = dn.modelrun.ModelRun(grid, year=2022, month=4, day=1)
    model.import_current(dn.read.current.cmems.Global(), expansion_factor=1.5)

    assert np.all(model.current().time() == timevec)

    grid_is_covered(grid, model.current())
    cleanup()


@pytest.mark.remote
def test_cmems_europe(grid2, timevec):
    cleanup()
    model = dn.modelrun.ModelRun(grid2, year=2022, month=4, day=1)
    model.import_current(dn.read.current.cmems.EuropeNW())
    assert np.all(model.current().time() == timevec)

    grid_is_covered(grid2, model.current())
    cleanup()