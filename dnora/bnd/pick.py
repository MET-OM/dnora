import numpy as np
from abc import ABC, abstractmethod

# Import objects
from ..grd.grd_mod import Grid

# Import aux_funcsiliry functions
from .. import msg
from ..aux_funcs import min_distance, expand_area

class PointPicker(ABC):
    """PointPickers take in longitude and latitude values, and returns indeces
    of the chosen points."""
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, grid: Grid, bnd_lon, bnd_lat):
        return

class TrivialPicker(PointPicker):
    """Choose all the points in the list."""
    def __init__(self):
        pass

    def __call__(self, grid: Grid, bnd_lon, bnd_lat, expansion_factor):
        inds = np.array(range(len(bnd_lon)))
        return inds

class NearestGridPoint(PointPicker):
    """Choose the nearest grid point to each boundary point in the grid.
    Set a maximum allowed distance using `max_dist` (in km) at instantiation time.
    """
    def __init__(self, max_dist=None):
        self.max_dist = max_dist
        pass

    def __call__(self, grid, bnd_lon, bnd_lat, expansion_factor):
        lon, lat = grid.boundary_points()

        # Go through all points where we want output and find the nearest available point
        inds = []
        for n in range(len(lat)):
            dx, ind = min_distance(lon[n], lat[n], bnd_lon, bnd_lat)
            ms = f"Point {n}: lat: {lat[n]:10.7f}, lon: {lon[n]:10.7f} <<< ({bnd_lat[ind]: .7f}, {bnd_lon[ind]: .7f}). Distance: {dx:.1f} km"
            if self.max_dist is None or dx <= self.max_dist:
                msg.plain(ms)
                inds.append(ind)
            else:
                msg.plain('DISCARDED, too far: '+ms)

        inds = np.array(inds)
        return inds

class Area(PointPicker):
    """Choose all the points within a certain area around the grid."""
    def __call__(self, grid: Grid, bnd_lon, bnd_lat, expansion_factor):
        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")

        # Define area to search in
        lon, lat = expand_area(grid.edges('lon'), grid.edges('lat'), expansion_factor)

        masklon = np.logical_and(bnd_lon > lon[0], bnd_lon < lon[1])
        masklat = np.logical_and(bnd_lat > lat[0], bnd_lat < lat[1])
        mask=np.logical_and(masklon, masklat)

        inds = np.where(mask)[0]

        msg.info(f"Found {len(inds)} points inside {lon[0]:10.7f}-{lon[1]:10.7f}, {lat[0]:10.7f}-{lat[1]:10.7f}.")

        return inds
