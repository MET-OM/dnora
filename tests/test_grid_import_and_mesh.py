from dnora.grd import Grid
from dnora.grd.read import ConstantTopo
import numpy as np
def test_constant_depth():
    grid = Grid(lon=(1,2), lat=(0,3))
    grid.set_spacing(nx=10, ny=10)
    grid.import_topo(ConstantTopo(depth=10))
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(),10)
    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(),10)
    np.testing.assert_almost_equal(grid.topo().max(),10)
    np.testing.assert_almost_equal(grid.topo().min(),10)