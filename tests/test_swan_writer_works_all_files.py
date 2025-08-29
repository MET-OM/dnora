import dnora as dn
import pytest

from pathlib import Path
import os
import shutil


@pytest.fixture(scope="session")
def model():
    grid = dn.grid.Grid(lon=(6, 7), lat=(60, 61))
    grid.set_spacing(dlon=0.1, dlat=0.1)
    model = dn.modelrun.Constant(grid, year=2020, month=1, day=1)

    model.import_ice()
    model.import_current()
    model.import_waterlevel()
    return model


def cleanup():
    if os.path.isdir("GridName_SWAN"):
        shutil.rmtree("GridName_SWAN")


def test_swan_writer_wind(model):
    model.import_wind()
    exp = dn.export.SWAN(model)
    fn = "swan_test_wind.asc"
    exp.export_wind(filename=fn)
    filepath = Path(f"GridName_SWAN/{fn}")
    assert filepath.is_file()
    assert os.path.getsize(filepath) == 29808
    cleanup()


def test_swan_writer_current(model):
    model.import_current()
    exp = dn.export.SWAN(model)
    fn = "swan_test_current.asc"
    exp.export_current(filename=fn)
    filepath = Path(f"GridName_SWAN/{fn}")
    assert filepath.is_file()
    assert os.path.getsize(filepath) == 29808
    cleanup()


def test_swan_writer_waterlevel(model):
    model.import_waterlevel()
    exp = dn.export.SWAN(model)
    fn = "swan_test_waterlevel.asc"
    exp.export_waterlevel(filename=fn)
    filepath = Path(f"GridName_SWAN/{fn}")
    assert filepath.is_file()
    assert os.path.getsize(filepath) == 14904
    cleanup()


def test_swan_writer_ice(model):
    model.import_ice(sit=10)
    exp = dn.export.SWAN(model)
    fn = "swan_test_ice.asc"
    exp.export_ice(filename=fn)
    filepath = Path(f"GridName_SWAN/{fn}")
    assert filepath.is_file()
    assert os.path.getsize(filepath) == 14904

    fn = "swan_test_sit.asc"
    exp.export_ice(filename=fn, data_vars=["sit"])
    filepath = Path(f"GridName_SWAN/{fn}")
    assert filepath.is_file()
    assert os.path.getsize(filepath) == 17808
    cleanup()
