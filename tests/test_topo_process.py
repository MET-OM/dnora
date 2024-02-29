from dnora.grid import Grid
from dnora.readers.generic_readers import ConstantGriddedData
from dnora.grid.process import SetConstantDepth, SetMinDepth, SetMaxDepth
import numpy as np


def test_import_empty():
    grid = Grid(lon=(1, 2), lat=(-1, 3))
    grid.set_spacing(nx=5, ny=5)
    topo_reader = ConstantGriddedData(topo=999.0)
    grid.import_topo(topo_reader)
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 999)
    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 999)
    assert np.all(grid.sea_mask())
    assert np.all(np.logical_not(grid.land_mask()))
    assert grid.is_empty("topo")
    grid.x()
    grid.y()


def test_import_constant():
    grid = Grid(lon=(1, 2), lat=(0, 3))
    grid.set_spacing(nx=10, ny=10)
    topo_reader = ConstantGriddedData(topo=100.0)
    grid.import_topo(topo_reader)
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 100)
    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 100)
    assert np.all(grid.sea_mask())
    assert np.all(np.logical_not(grid.land_mask()))
    assert not grid.is_empty("topo")


def test_import_land():
    grid = Grid(lon=(1, 2), lat=(0, 3))
    grid.set_spacing(nx=10, ny=10)
    topo_reader = ConstantGriddedData(topo=0.0)
    grid.import_topo(topo_reader)
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 0)
    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 0)
    assert np.all(grid.land_mask())
    assert np.all(np.logical_not(grid.sea_mask()))
    assert grid.is_empty("topo")


def test_process_constant():
    grid = Grid(lon=(1, 2), lat=(0, 3))
    grid.set_spacing(nx=10, ny=10)
    topo_reader = ConstantGriddedData(topo=999.0)
    grid.import_topo(topo_reader)
    grid.process_grid(SetConstantDepth(depth=10), raw=True)
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 10)
    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 10)


def test_process_min_depth():
    grid = Grid(lon=(1, 2), lat=(0, 3))
    grid.set_spacing(nx=10, ny=10)
    topo_reader = ConstantGriddedData(topo=999.0)
    grid.import_topo(topo_reader)

    grid.process_grid(SetMinDepth(depth=1000), raw=True)
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 1000)
    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 1000)

    grid.process_grid(SetMinDepth(depth=10), raw=True)
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 1000)
    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 1000)


def test_process_max_depth():
    grid = Grid(lon=(1, 2), lat=(0, 3))
    grid.set_spacing(nx=10, ny=10)
    topo_reader = ConstantGriddedData(topo=999.0)
    grid.import_topo(topo_reader)

    grid.process_grid(SetMaxDepth(depth=10), raw=True)
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 10)
    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 10)

    grid.process_grid(SetMaxDepth(depth=1000), raw=True)
    np.testing.assert_array_almost_equal(grid.raw().topo().mean(), 10)
    grid.mesh_grid()
    np.testing.assert_array_almost_equal(grid.topo().mean(), 10)
