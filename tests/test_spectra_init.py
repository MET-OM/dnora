from dnora.grd import Grid
from dnora.spc import Spectra
def test_spectra_one_point():
    grid = Grid(lon=5, lat=60)
    spectra = Spectra(grid=grid)
    spectra.convention
    print(spectra)

def test_spectra_one_area():
    grid = Grid(lon=(5,6), lat=(60,61))
    spectra = Spectra(grid=grid)
    spectra.convention
    print(spectra)

def test_spectra_gridded():
    grid = Grid(lon=(5,6), lat=(60,61))
    grid.set_spacing(nx=10, ny=10)
    spectra = Spectra(grid=grid)
    spectra.convention
    print(spectra)

def test_spectra_one_point_cartesian():
    grid = Grid(x=5, y=3)
    spectra = Spectra(grid=grid)
    spectra.convention
    print(spectra)

def test_spectra_one_area_cartesian():
    grid = Grid(x=(5,6), y=(3,4))
    spectra = Spectra(grid=grid)
    spectra.convention
    print(spectra)

def test_spectra_gridded_cartesian():
    grid = Grid(x=(5,6), y=(3,4))
    grid.set_spacing(nx=10, ny=10)
    spectra = Spectra(grid=grid)
    spectra.convention
    print(spectra)
