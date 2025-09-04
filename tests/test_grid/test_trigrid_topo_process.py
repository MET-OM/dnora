from dnora.grid import TriGrid
from dnora.read.generic import ConstantData
from dnora.grid.process import SetMinDepth, SetMaxDepth
import numpy as np


def test_import_empty():
    grid = TriGrid(
        lon=[0, 1, 2, 3], lat=[6, 7, 8, 9], ntriang=range(5), corner=range(3)
    )
    topo_reader = ConstantData(vars={"topo": 999.0})
    grid.import_topo(topo_reader)
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 999)
    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 999)
    assert np.all(grid.sea_mask())
    assert np.all(np.logical_not(grid.land_mask()))


def test_import_land():
    grid = TriGrid(
        lon=[0, 1, 2, 3], lat=[6, 7, 8, 9], ntriang=range(5), corner=range(3)
    )
    topo_reader = ConstantData(vars={"topo": 0.0})
    grid.import_topo(topo_reader)
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 0)
    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 0)
    assert np.all(grid.land_mask())
    assert np.all(np.logical_not(grid.sea_mask()))


def test_min_depth():
    grid = TriGrid(
        lon=[0, 1, 2, 3], lat=[6, 7, 8, 9], ntriang=range(5), corner=range(3)
    )
    topo_reader = ConstantData(vars={"topo": 999.0})
    grid.import_topo(topo_reader)
    grid.process_grid(SetMinDepth(min_depth=1000), raw=True)
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 1000)

    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 1000)
    assert np.all(grid.sea_mask())
    assert np.all(np.logical_not(grid.land_mask()))


def test_max_depth():
    grid = TriGrid(
        lon=[0, 1, 2, 3], lat=[6, 7, 8, 9], ntriang=range(5), corner=range(3)
    )
    topo_reader = ConstantData()
    grid.import_topo(topo_reader, topo=999.0)
    grid.process_grid(SetMaxDepth(max_depth=100), raw=True)
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 100)

    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 100)
    assert np.all(grid.sea_mask())
    assert np.all(np.logical_not(grid.land_mask()))
