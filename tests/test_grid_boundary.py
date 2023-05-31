from dnora.grd import Grid
from dnora.grd.mask import SetAll, ClearBoundary
import numpy as np
def test_import_empty():
    grid = Grid(lon=(1,2), lat=(0,3))
    grid.set_spacing(nx=10, ny=5)
    assert np.all(np.logical_not(grid.boundary_mask()))

    grid.set_mask(SetAll())
    assert grid.boundary_nx() == grid.nx()
    assert grid.boundary_ny() == grid.ny()
    assert np.all(grid.boundary_mask())
    assert len(grid.boundary_points()[0]) == grid.nx()*grid.ny()
    assert len(grid.boundary_points()[1]) == grid.nx()*grid.ny()


    grid.set_mask(ClearBoundary())
    assert np.all(np.logical_not(grid.boundary_mask()))
    assert grid.boundary_nx() == 0
    assert grid.boundary_nx() == 0
    assert len(grid.boundary_points()[0]) == 0
    assert len(grid.boundary_points()[1]) == 0
