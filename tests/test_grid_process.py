from dnora.grd import Grid
from dnora.grd.process import SetConstantDepth, SetMinDepth, SetMaxDepth
import numpy as np
def test_constant_depth():
    grid = Grid(lon=(1,2), lat=(0,3))
    grid.set_spacing(nx=10, ny=10)
    grid.process_grid(SetConstantDepth(depth=10))
    np.testing.assert_array_almost_equal(grid.topo().mean(),10)

def test_min_depth():
    grid = Grid(lon=(1,2), lat=(0,3))
    grid.set_spacing(nx=10, ny=10)
    grid.process_grid(SetMinDepth(depth=1000.))
    np.testing.assert_array_almost_equal(grid.topo().mean(),1000.)
    grid.process_grid(SetMinDepth(depth=10.))
    np.testing.assert_array_almost_equal(grid.topo().mean(),1000.)

def test_min_depth():
    grid = Grid(lon=(1,2), lat=(0,3))
    grid.set_spacing(nx=10, ny=10)
    grid.process_grid(SetMaxDepth(depth=100.))
    np.testing.assert_array_almost_equal(grid.topo().mean(),100.)
    grid.process_grid(SetMaxDepth(depth=1000.))
    np.testing.assert_array_almost_equal(grid.topo().mean(),100.)
