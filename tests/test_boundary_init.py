from dnora.grd import Grid
from dnora.bnd import Boundary
def test_boundary_one_point():
    grid = Grid(lon=5, lat=60)
    boundary = Boundary(grid=grid)
    boundary.convention
    print(boundary)

def test_boundary_one_area():
    grid = Grid(lon=(5,6), lat=(60,61))
    boundary = Boundary(grid=grid)
    boundary.convention
    print(boundary)

def test_boundary_gridded():
    grid = Grid(lon=(5,6), lat=(60,61))
    grid.set_spacing(nx=10, ny=10)
    boundary = Boundary(grid=grid)
    boundary.convention
    print(boundary)

def test_boundary_one_point_cartesian():
    grid = Grid(x=5, y=3)
    boundary = Boundary(grid=grid)
    boundary.convention
    print(boundary)

def test_boundary_one_area_cartesian():
    grid = Grid(x=(5,6), y=(3,4))
    boundary = Boundary(grid=grid)
    boundary.convention
    print(boundary)

def test_boundary_gridded_cartesian():
    grid = Grid(x=(5,6), y=(3,4))
    grid.set_spacing(nx=10, ny=10)
    boundary = Boundary(grid=grid)
    boundary.convention
    print(boundary)
