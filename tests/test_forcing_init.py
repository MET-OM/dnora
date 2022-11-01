from dnora.grd import Grid
from dnora.wnd import Forcing
def test_forcing_one_point():
    grid = Grid(lon=5, lat=60)
    forcing = Forcing(grid=grid)
    print(forcing)

def test_forcing_one_area():
    grid = Grid(lon=(5,6), lat=(60,61))
    forcing = Forcing(grid=grid)
    print(forcing)

def test_forcing_gridded():
    grid = Grid(lon=(5,6), lat=(60,61))
    grid.set_spacing(nx=10, ny=10)
    forcing = Forcing(grid=grid)
    print(forcing)

def test_forcing_one_point_cartesian():
    grid = Grid(x=5, y=3)
    forcing = Forcing(grid=grid)
    print(forcing)

def test_forcing_one_area_cartesian():
    grid = Grid(x=(5,6), y=(3,4))
    forcing = Forcing(grid=grid)
    print(forcing)

def test_forcing_gridded_cartesian():
    grid = Grid(x=(5,6), y=(3,4))
    grid.set_spacing(nx=10, ny=10)
    forcing = Forcing(grid=grid)
    print(forcing)
