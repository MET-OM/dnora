import numpy as np
import pandas as pd
import xarray as xr
import utm
from copy import copy
from skeleton import Skeleton
from dataset_manager import DatasetManager

class PointSkeleton(Skeleton):
    """Gives a unstructured structure to the Skeleton.

    In practise this means that:

    1) Grid coordinates are defined with and index (inds),
    2) x,y / lon/lat values are data variables of the index.
    3) Methods x(), y() / lon(), lat() will returns all points.
    4) Methods xy() / lonlat() are identical to e.g. (x(), y()).
    """

    # def __init__(self, x=None, y=None, lon=None, lat=None, time=None, name='PointyData'):
    #     self.data = super()._create_structure(x, y, lon, lat, time)
    #     self.data.attrs['name'] = name

    def _initial_coords(self) -> list[str]:
        return ['inds']

    def _initial_vars(self) -> dict:
        return {'x': 'inds', 'y': 'inds'}

    def lonlat(self, mask: np.array=None, strict=False, **kwargs) -> tuple[np.ndarray, np.ndarray]:
        """Returns a tuple of longitude and latitude of all points.
        If strict=True, then None is returned if grid is cartesian.

        Identical to (.lon(), .lat()) (with no mask)
        mask is a boolean array (default True for all points)
        """
        if mask is None:
            mask = np.full(super().size('spatial', **kwargs), True)

        lon, lat = super().lon(strict=strict, **kwargs)[mask], super().lat(strict=strict, **kwargs)[mask]

        if lon is None:
            return None, None

        return lon[mask], lat[mask]

    def xy(self, mask: np.array=None, strict=False, **kwargs) -> tuple[np.ndarray, np.ndarray]:
        """Returns a tuple of x and y of all points.
        If strict=True, then None is returned if grid is sperical.

        Identical to (.x(), .y()) (with no mask)
        mask is a boolean array (default True for all points)
        """
        if mask is None:
            mask = np.full(super().size('spatial', **kwargs), True)

        x, y = super().x(strict=strict, **kwargs)[mask], super().y(strict=strict, **kwargs)[mask]

        if x is None:
            return None, None

        return x[mask], y[mask]

    def native_xy(self, mask: np.array=None, **kwargs) -> tuple[np.ndarray, np.ndarray]:
        """Returns a tuple of native x and y of all points.
        Identical to (.native_x(), .native_y()) (with no mask)

        mask is a boolean array (default True for all points)
        """
        if mask is None:
            mask = np.full(super().size('spatial', **kwargs), True)

        return super().native_x(**kwargs)[mask], super().native_y(**kwargs)[mask]
