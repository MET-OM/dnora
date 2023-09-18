from dnora import grd, mdl, bnd
import numpy as np

def test_import_constant_spectra_one_point():
    grid = grd.Grid(lon=5, lat=60)
    model = mdl.ModelRun(grid, start_time = '2020-01-01 00:00', end_time='2020-01-02 00:00')
    
    model.set_spectral_grid()
