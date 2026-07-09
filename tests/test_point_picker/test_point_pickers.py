from dnora.grid import Grid

from dnora import pick
from geo_skeletons import PointSkeleton, GriddedSkeleton
import numpy as np
import pytest


@pytest.fixture(scope="session")
def grid():
    grid = Grid(lon=(10, 20), lat=(60, 65))
    grid.set_spacing(nx=20, ny=20)
    return grid


@pytest.fixture(scope="session")
def all_points(grid):
    all_points = PointSkeleton.from_skeleton(grid)
    return all_points


def test_trivial_picker(grid, all_points):
    inds = pick.Trivial()(grid, all_points)
    np.testing.assert_array_almost_equal(inds, all_points.inds())


def test_area_picker(grid, all_points):
    small_grid = Grid(lon=(15, 16), lat=(62, 63))
    inds = pick.Area()(small_grid, all_points, expansion_factor=1.0)

    assert np.all(all_points.lon(inds=list(inds)) <= small_grid.lon()[-1])
    assert np.all(all_points.lon(inds=list(inds)) >= small_grid.lon()[0])
    assert np.all(all_points.lat(inds=list(inds)) <= small_grid.lat()[-1])
    assert np.all(all_points.lat(inds=list(inds)) >= small_grid.lat()[0])


def test_area_picker_expand(all_points):
    small_grid = Grid(lon=(15, 16), lat=(62, 63))
    inds = pick.Area()(small_grid, all_points, expansion_factor=2.0)

    assert np.any(all_points.lon(inds=list(inds)) > small_grid.lon()[-1])
    assert np.any(all_points.lon(inds=list(inds)) < small_grid.lon()[0])
    assert np.any(all_points.lat(inds=list(inds)) > small_grid.lat()[-1])
    assert np.any(all_points.lat(inds=list(inds)) < small_grid.lat()[0])
    assert np.all(all_points.lon(inds=list(inds)) <= 17)
    assert np.all(all_points.lon(inds=list(inds)) >= 14)
    assert np.all(all_points.lat(inds=list(inds)) <= 64)
    assert np.all(all_points.lat(inds=list(inds)) >= 61)


def test_selected_points(grid, all_points):
    selected_points = PointSkeleton.from_skeleton(
        all_points, mask=all_points.lon() < 10.01
    )
    inds = pick.NearestGridPoint()(grid, all_points, selected_points)
    assert max(all_points.lon(inds=list(inds))) == 10
    assert min(all_points.lon(inds=list(inds))) == 10
    lat0, lat1 = selected_points.edges("lat")
    assert min(all_points.lat(inds=list(inds))) == lat0
    assert max(all_points.lat(inds=list(inds))) == lat1


def test_area_180_wrap():
    all_points = PointSkeleton(lon=(-175, 179), lat=(20,20))
    grid = GriddedSkeleton(lon=(-179,-160), lat=(19, 21))
    inds = pick.Area()(grid=grid, all_points=all_points, expansion_factor=1)
    assert len(inds) == 1
    assert inds[0] == 0
    inds = pick.Area()(grid=grid, all_points=all_points, expansion_factor=2)
    np.testing.assert_array_almost_equal(inds,[0,1])

def test_nearest_180_wrap():
    # Nearest point wraps around
    all_points = PointSkeleton(lon=(-179, 175), lat=(20,20))
    selected_points = PointSkeleton(lon=(175.1,179.9), lat=(20, 20))
    inds = pick.NearestGridPoint()(grid=None, all_points=all_points, selected_points=selected_points)
    np.testing.assert_array_almost_equal(inds,[0,1])
    
    # Nearest point is not wrapped around
    all_points = PointSkeleton(lon=(-149, 175), lat=(20,20))
    selected_points = PointSkeleton(lon=(175.1,179.9), lat=(20, 20))
    inds = pick.NearestGridPoint()(grid=None, all_points=all_points, selected_points=selected_points)
    np.testing.assert_array_almost_equal(inds,[1])