import dnora as dn
import pytest


@pytest.fixture(scope="session")
def grid():
    return dn.grid.Grid(lon=(10, 14), lat=(60, 61))


@pytest.mark.online
def test_nora3(grid):
    model = dn.modelrun.NORA3(grid, year=2022, month=4, day=1)
    model.import_wind()
