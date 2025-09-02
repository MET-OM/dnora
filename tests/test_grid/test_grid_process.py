from dnora.grid import Grid
from dnora.grid.process import SetConstantDepth, SetMinDepth, SetMaxDepth
import numpy as np


def test_constant_depth():
    grid = Grid(lon=(1, 2), lat=(0, 3))
    grid.set_spacing(nx=10, ny=10)
    grid.process_grid(SetConstantDepth(depth=10))
    np.testing.assert_array_almost_equal(grid.topo().mean(), 10)


def test_min_depth():
    grid = Grid(lon=(1, 2), lat=(0, 3))
    grid.set_spacing(nx=10, ny=10)
    grid.process_grid(SetMinDepth(min_depth=1000.0))
    np.testing.assert_array_almost_equal(grid.topo().mean(), 1000.0)
    grid.process_grid(SetMinDepth(), min_depth=10.0)
    np.testing.assert_array_almost_equal(grid.topo().mean(), 1000.0)


def test_min_depth():
    grid = Grid(lon=(1, 2), lat=(0, 3))
    grid.set_spacing(nx=10, ny=10)
    grid.process_grid(SetMaxDepth(max_depth=100.0))
    np.testing.assert_array_almost_equal(grid.topo().mean(), 100.0)
    grid.process_grid(SetMaxDepth(), max_depth=1000.0)
    np.testing.assert_array_almost_equal(grid.topo().mean(), 100.0)
