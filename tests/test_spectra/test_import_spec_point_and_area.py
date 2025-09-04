import dnora as dn
import pytest
import numpy as np

@pytest.mark.remote
def test_read_one_point():
    grid = dn.grid.Grid(lon=1,lat=61)
    model = dn.modelrun.NORA3(grid, start_time='2020-01-01 00:00', end_time='2020-01-01 01:00') 
    model.import_spectra()
    assert len(model.spectra().inds()) == 1

@pytest.mark.remote
def test_read_one_point_discarded():
    grid = dn.grid.Grid(lon=1,lat=61)
    model = dn.modelrun.NORA3(grid, start_time='2020-01-01 00:00', end_time='2020-01-01 01:00') 
    model.import_spectra(max_dist=10)
    assert model.spectra() is None

@pytest.mark.remote
def test_read_area():
    grid = dn.grid.Grid(lon=(1,2),lat=(61,62))
    model = dn.modelrun.NORA3(grid, start_time='2020-01-01 00:00', end_time='2020-01-01 01:00') 
    model.import_spectra(expansion_factor=1)
    assert np.max(model.spectra().lon()) < 2
    assert np.min(model.spectra().lon()) > 1
    assert np.max(model.spectra().lat()) < 62
    assert np.min(model.spectra().lat()) > 11
