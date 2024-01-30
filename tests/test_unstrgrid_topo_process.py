from dnora.grid import TriGrid
from dnora.grid.read import ConstantTopo
from dnora.grid.process import SetConstantDepth, SetMinDepth, SetMaxDepth
import numpy as np


def test_import_empty():
    grid = TriGrid(lon=[0, 1, 2, 3], lat=[6, 7, 8, 9])
    topo_reader = ConstantTopo()
    grid.import_topo(topo_reader)
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 999)
    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 999)
    assert np.all(grid.sea_mask())
    assert np.all(np.logical_not(grid.land_mask()))
    assert grid.is_empty("topo")


def test_import_land():
    grid = TriGrid(lon=[0, 1, 2, 3], lat=[6, 7, 8, 9])
    topo_reader = ConstantTopo(depth=0.0)
    grid.import_topo(topo_reader)
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 0)
    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 0)
    assert np.all(grid.land_mask())
    assert np.all(np.logical_not(grid.sea_mask()))
    assert grid.is_empty("topo")


def test_min_depth():
    grid = TriGrid(lon=[0, 1, 2, 3], lat=[6, 7, 8, 9])
    topo_reader = ConstantTopo()
    grid.import_topo(topo_reader)
    grid.raw().process_grid(SetMinDepth(depth=1000))
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 1000)

    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 1000)
    assert np.all(grid.sea_mask())
    assert np.all(np.logical_not(grid.land_mask()))
    assert not grid.is_empty("topo")


def test_max_depth():
    grid = TriGrid(lon=[0, 1, 2, 3], lat=[6, 7, 8, 9])
    topo_reader = ConstantTopo()
    grid.import_topo(topo_reader)
    grid.raw().process_grid(SetMaxDepth(depth=100))
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 100)

    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 100)
    assert np.all(grid.sea_mask())
    assert np.all(np.logical_not(grid.land_mask()))
    assert not grid.is_empty("topo")
