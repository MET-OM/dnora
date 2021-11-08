import numpy as np
from .. import msg
from ..aux import min_distance

from ..bnd_mod import PointPicker # Abstract class


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
    def __init__(self, expansion_factor = 1.5):
        self.expansion_factor = expansion_factor
        return

    def __call__(self, grid, bnd_lon, bnd_lat):
        msg.info(f"Using expansion_factor = {self.expansion_factor:.2f}")
        # Define area to search in
        expand_lon = (grid.lon()[-1] - grid.lon()[0])*(self.expansion_factor-1)*0.5
        expand_lat = (grid.lat()[-1] - grid.lat()[0])*(self.expansion_factor-1)*0.5

        # Get all the spectra in this area
        lon0=grid.lon()[0] - expand_lon
        lon1=grid.lon()[-1] + expand_lon

        lat0=grid.lat()[0] - expand_lat
        lat1=grid.lat()[-1] + expand_lat

        masklon = np.logical_and(bnd_lon < lon1, bnd_lon > lon0)
        masklat = np.logical_and(bnd_lat > lat0, bnd_lat < lat1)
        mask=np.logical_and(masklon, masklat)

        inds = np.where(mask)[0]

        msg.info(f"Found {len(inds)} points inside {lon0:10.7f}-{lon1:10.7f}, {lat0:10.7f}-{lat1:10.7f}.")

        return inds

