import numpy as np
import pandas as pd
import xarray as xr
import utm
from copy import copy
from skeleton import Skeleton

class GriddedSkeleton(Skeleton):
    def __init__(self, x=None, y=None, lon=None, lat=None, time=None, name='GriddedData'):
        self.data = super()._create_structure(x, y, lon, lat, time)
        self.data.attrs['name'] = name

    def _init_ds(self, x: np.ndarray, y: np.ndarray, time=None) -> xr.Dataset:
        coords_dict = {self.y_str: y, self.x_str: x}
        if time is not None:
            coords_dict['time'] = time
        return xr.Dataset(coords=coords_dict, attrs={'name': self.name})

    def _ds_coords_dict(self):
        """Return coordinate dictionary for creating xarray Dataset"""
        coords_dict = {self.y_str: super().native_y(), self.x_str: super().native_x()}
        if super().time() is not None:
            coords_dict['time'] = super().time()
        return coords_dict

    def size(self) -> tuple[int, int]:
        """Returns the size (nx, ny) of the grid."""
        return (super().ny(), super().nx())

    def lonlat(self, mask: np.array=None, order_by: str='lat') -> tuple[np.ndarray, np.ndarray]:
        """Returns a tuple of longitude and latitude of all points.

        mask is a boolean array (default True for all points)
        order_by = 'lat' (default) or 'lon'
        """
        if mask is None:
            mask = np.full((super().nx(), super().ny()), True)
        mask = mask.ravel()

        lon, lat = self.native_xy(mask, order_by)

        # Transforms lon/lat to x/y if necessary
        lon, lat = super()._lonlat(lon, lat)

        return lon[mask], lat[mask]

    def xy(self, mask: np.array=None, order_by: str='y') -> tuple[np.ndarray, np.ndarray]:
        """Returns a tuple of x and y of all points.

        mask is a boolean array (default True for all points)
        order_by = 'y' (default) or 'x'
        """
        if mask is None:
            mask = np.full((super().nx(), super().ny()), True)
        mask = mask.ravel()

        x, y = self.native_xy(mask, order_by)

        # Transforms lon/lat to x/y if necessary
        x, y = super()._xy(x, y)

        return x[mask], y[mask]


    def native_xy(self, mask: np.array=None, order_by: str='y') -> tuple[np.ndarray, np.ndarray]:
        """Returns a tuple of native x and y of all points.

        mask is a boolean array (default True for all points)
        order_by = 'y' (default) or 'x'
        """
        if order_by == 'y' or order_by == 'lat':
            x, y = np.meshgrid(super().native_x(), super().native_y())
        elif order_by == 'x' or order_by == 'lon':
            y, x = np.meshgrid(super().native_y(), super().native_x())
        else:
            raise ValueError("order_by should be 'y' (/'lat') or 'x' (/'lon')")

        if mask is None:
            mask = np.full((super().nx(), super().ny()), True)
        mask = mask.ravel()

        return x.ravel()[mask], y.ravel()[mask]
