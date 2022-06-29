import numpy as np
import pandas as pd
import xarray as xr
import utm
from copy import copy
from skeleton import Skeleton

class PointSkeleton(Skeleton):
    def __init__(self, x=None, y=None, lon=None, lat=None, time=None, name='PointyData'):
        self.data = super()._create_structure(x, y, lon, lat, time)
        self.data.attrs['name'] = name

    def _init_ds(self, x: np.ndarray, y: np.ndarray, time=None, **kwargs) -> xr.Dataset:
        coords_dict = {'inds': np.arange(len(x))}

        if time is not None:
            coords_dict['time'] = (['inds'], time)

        for key, value in kwargs.items():
            coords_dict[key] = value

        vars_dict = {self.x_str: (['inds'], x), self.y_str: (['inds'], y)}

        return xr.Dataset(coords=coords_dict, data_vars=vars_dict, attrs={'name': self.name})

    def _ds_coords_dict(self):
        """Return coordinate dictionary for creating xarray Dataset"""
        coords_dict = {'inds': self.inds()}
        if super().time() is not None:
            coords_dict['time'] = super().time()
        return coords_dict

    def _ds_vars_dict(self):
        """Return variable dictionary for creating xarray Dataset"""
        vars_dict = {self.x_str: (['inds'], super().native_x()), self.y_str: (['inds'], super().native_y())}
        return vars_dict

    def inds(self) -> np.ndarray:
        if hasattr(self.data, 'inds') :
            return self.data.inds.values
        else:
            return None

    def size(self) -> tuple[int]:
        """Returns the size of the object.

        Spatial, temporal and possible added dimensions."""
        list = [super().nx()]
        if self.nt() is not None:
            list.append(super().nt())

        for coord in self._additional_coords():
            if self._additional_coord_val(coord) is not None:
                list.append(len(self._additional_coord_val(coord)))
        return tuple(list)

    def lonlat(self, mask: np.array=None) -> tuple[np.ndarray, np.ndarray]:
        """Returns a tuple of longitude and latitude of all points.
        Identical to (.lon(), .lat()) (with no mask)

        mask is a boolean array (default True for all points)
        """
        if mask is None:
            mask = np.full((self.nx(),), True)

        return super().lon()[mask], super().lat()[mask]

    def xy(self, mask: np.array=None) -> tuple[np.ndarray, np.ndarray]:
        """Returns a tuple of x and y of all points.
        Identical to (.x(), .y()) (with no mask)

        mask is a boolean array (default True for all points)
        """
        if mask is None:
            mask = np.full((super().nx(),), True)

        return super().x()[mask], super().y()[mask]

    def native_xy(self, mask: np.array=None, **kwargs) -> tuple[np.ndarray, np.ndarray]:
        """Returns a tuple of native x and y of all points.
        Identical to (.native_x(), .native_y()) (with no mask)

        mask is a boolean array (default True for all points)
        """
        if mask is None:
            mask = np.full((super().nx(),), True)

        return super().native_x()[mask], super().native_y()[mask]
