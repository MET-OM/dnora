import numpy as np
from abc import ABC, abstractmethod
from typing import Union
# Import objects
from ..grd.grd_mod import Grid, UnstrGrid

# Import aux_funcsiliry functions
from .. import msg
from ..aux_funcs import min_distance, expand_area

class PointPicker(ABC):
    """PointPickers take in longitude and latitude values, and returns indeces
    of the chosen points."""
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, grid: Union[Grid, UnstrGrid], all_points: UnstrGrid, expansion_factor: float) -> np.ndarray:
        return inds

class TrivialPicker(PointPicker):
    """Choose all the points in the list."""
    def __init__(self):
        pass

    def __call__(self, grid: Union[Grid, UnstrGrid],
                    all_points: UnstrGrid,
                    expansion_factor: float) -> np.ndarray:
        return all_points.inds()

class NearestGridPoint(PointPicker):
    """Choose the nearest grid point to each boundary point in the grid.
    Set a maximum allowed distance using `max_dist` (in km) at instantiation time.
    """
    def __init__(self, max_dist=None):
        self.max_dist = max_dist
        pass

    def __call__(self, grid: Union[Grid, UnstrGrid],
                    all_points: UnstrGrid,
                    expansion_factor: float) -> np.ndarray:

        lon, lat = grid.boundary_points('lon')

        # Go through all points where we want output and find the nearest available point
        inds = []
        for n in range(len(lat)):
            dx, ind = min_distance(lon[n], lat[n], all_points.lon(), all_points.lat())
            ms = f"Point {n}: lat: {lat[n]:10.7f}, lon: {lon[n]:10.7f} <<< ({all_points.lat()[ind]: .7f}, {all_points.lon()[ind]: .7f}). Distance: {dx:.1f} km"
            if self.max_dist is None or dx <= self.max_dist:
                msg.plain(ms)
                inds.append(ind)
            else:
                msg.plain('DISCARDED, too far: '+ms)

        inds = np.array(inds)
        return inds

class Area(PointPicker):
    """Choose all the points within a certain area around the grid."""
    def __call__(self, grid: Union[Grid, UnstrGrid],
                    all_points: UnstrGrid,
                    expansion_factor: float) -> np.ndarray:

        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")

        # Define area to search in
        if grid.is_cartesian():
            number, zone = grid.utm()
            all_points.set_utm(number, zone)
            x, y = expand_area(grid.edges('x'), grid.edges('y'), expansion_factor)
            x_all, y_all = all_points.xy()
        else:
            x, y = expand_area(grid.edges('lon'), grid.edges('lat'), expansion_factor)
            x_all, y_all = all_points.lonlat()

        maskx = np.logical_and(x_all > x[0], x_all < x[1])
        masky = np.logical_and(y_all > y[0], y_all < y[1])
        mask = np.logical_and(maskx, masky)

        inds = np.where(mask)[0]

        msg.info(f"Found {len(inds)} points inside {x[0]:10.7f}-{x[1]:10.7f}, {y[0]:10.7f}-{y[1]:10.7f}.")

        return inds
