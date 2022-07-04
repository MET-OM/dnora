import numpy as np
import pandas as pd
import xarray as xr
import utm
from copy import copy
from skeleton import Skeleton
from dataset_manager import DatasetManager

class GriddedSkeleton(Skeleton):
    """Gives a gridded structure to the Skeleton.

    In practise this means that:

    1) Grid coordinates are defined as x,y / lon,lat.
    2) Methods x(), y() / lon(), lat() will return the vectors defining the grid.
    3) Methods xy() / lonlat() will return a list of all points of the grid
    (i.e. raveled meshgrid).
    """

    # def __init__(self, x=None, y=None, lon=None, lat=None, name='GriddedData'):
    #     self.data = super()._create_structure(x, y, lon, lat)
    #     self.data.attrs['name'] = name


    def _initial_coords(self) -> list[str]:
        return ['y', 'x']

    def _initial_vars(self) -> dict:
        return {}

    def lonlat(self, mask: np.array=None, order_by: str='lat', strict=False, **kwargs) -> tuple[np.ndarray, np.ndarray]:
        """Returns a tuple of longitude and latitude of all points.
        If strict=True, then None is returned if grid is cartesian.

        mask is a boolean array (default True for all points)
        order_by = 'lat' (default) or 'lon'
        """
        if mask is None:
            mask = np.full(super().size('spatial', **kwargs), True)
        mask = mask.ravel()
        lon, lat = self.native_xy(mask, order_by,**kwargs)

        # Transforms lon/lat to x/y if necessary
        lon, lat = super()._lonlat(lon, lat, strict=strict)
        if lon is None:
            return None, None

        return lon[mask], lat[mask]

    def xy(self, mask: np.array=None, order_by: str='y', strict=False, **kwargs) -> tuple[np.ndarray, np.ndarray]:
        """Returns a tuple of x and y of all points.
        If strict=True, then None is returned if grid is sperical.

        mask is a boolean array (default True for all points)
        order_by = 'y' (default) or 'x'
        """
        if mask is None:
            mask = np.full(super().size('spatial'), True)
        mask = mask.ravel()

        x, y = self.native_xy(mask, order_by, **kwargs)

        # Transforms lon/lat to x/y if necessary
        x, y = super()._xy(x, y, strict=strict)

        if x is None:
            return None, None

        return x[mask], y[mask]


    def native_xy(self, mask: np.array=None, order_by: str='y', **kwargs) -> tuple[np.ndarray, np.ndarray]:
        """Returns a tuple of native x and y of all points.

        mask is a boolean array (default True for all points)
        order_by = 'y' (default) or 'x'
        """
        if order_by == 'y' or order_by == 'lat':
            x, y = np.meshgrid(super().native_x(**kwargs), super().native_y(**kwargs))
        elif order_by == 'x' or order_by == 'lon':
            y, x = np.meshgrid(super().native_y(**kwargs), super().native_x(**kwargs))
        else:
            raise ValueError("order_by should be 'y' (/'lat') or 'x' (/'lon')")

        if mask is None:
            mask = np.full(super().size('spatial', **kwargs), True)
        mask = mask.ravel()

        return x.ravel()[mask], y.ravel()[mask]
