import numpy as np
from .. import msg
from ..aux import min_distance, expand_area
from abc import ABC, abstractmethod

from ..grd.grd_mod import Grid # Grid object

class PointPicker(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, grid: Grid, bnd_lon, bnd_lat):
        return

class TrivialPicker(PointPicker):
    def __init__(self):
        pass

    def __call__(self, grid: Grid, bnd_lon, bnd_lat):
        inds = np.array(range(len(bnd_lon)))
        return inds

class NearestGridPoint(PointPicker):
    def __init__(self):
        pass

    def __call__(self, grid, bnd_lon, bnd_lat):
        bnd_points = grid.boundary_points()
        lon = bnd_points[:,0]
        lat = bnd_points[:,1]

        # Go through all points where we want output and find the nearest available point
        inds = []
        for n in range(len(lat)):
            dx, ind = min_distance(lon[n], lat[n], bnd_lon, bnd_lat)
            msg.plain(f"Point {n}: lat: {lat[n]:10.7f}, lon: {lon[n]:10.7f} <<< ({bnd_lat[ind]: .7f}, {bnd_lon[ind]: .7f}). Distance: {dx:.1f} km")
            inds.append(ind)

        inds = np.array(inds)
        return inds

class Area(PointPicker):
    def __init__(self, expansion_factor=1.5):
        self.expansion_factor = expansion_factor
        return

    def __call__(self, grid: Grid, bnd_lon, bnd_lat):
        msg.info(f"Using expansion_factor = {self.expansion_factor:.2f}")

        # Define area to search in
        lon_min, lon_max, lat_min, lat_max = expand_area(grid.lon()[0], grid.lon()[-1], grid.lat()[0], grid.lat()[-1], self.expansion_factor)

        masklon = np.logical_and(bnd_lon > lon_min, bnd_lon < lon_max)
        masklat = np.logical_and(bnd_lat > lat_min, bnd_lat < lat_max)
        mask=np.logical_and(masklon, masklat)

        inds = np.where(mask)[0]

        msg.info(f"Found {len(inds)} points inside {lon_min:10.7f}-{lon_max:10.7f}, {lat_min:10.7f}-{lat_max:10.7f}.")

        return inds
