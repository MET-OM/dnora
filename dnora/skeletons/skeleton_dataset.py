import numpy as np
import pandas as pd
import xarray as xr
import utm
from copy import copy

class SkeletonDataset:
    """Contains methods related to the creataon and handling of the Xarray
    Dataset that will be used in any object that inherits from Skeleton."""

    _spatial_coord_list = ['x', 'y', 'lon', 'lat', 'inds']

    def _create_structure(self, x=None, y=None, lon=None, lat=None, **kwargs):
        """Create the first Dataset with either x,y (Cartesian) or lon, lat (Spherical)
        depending on which variables are provided."""
        def check_consistency():
            ds_coords = list(ds.coords)
            # Check spatial coordinates
            xy_set = 'x' in ds_coords and 'y' in ds_coords
            lonlat_set = 'lon' in ds_coords and 'lat' in ds_coords
            inds_set = 'inds' in ds_coords
            if not (xy_set or lonlat_set or inds_set):
                raise ValueError("A proper spatial grid is not set: Requires 'x' and 'y', 'lon' and 'lat' or 'inds'!")
            if sum([xy_set, lonlat_set, inds_set]) > 1:
                raise ValueError("A well defined spatial grid is not set: Requires 'x' and 'y', 'lon' and 'lat' or 'inds'!")

            # Check that all added coordinates are porvided
            for coord in self._added_coords():
                if coord not in ds_coords:
                    raise ValueError(f"Coordinate '{coord}' has been added (by a decorator?) but it was not provided when the Dataset ({ds_coords}) was created!")

            # Check that all provided coordinates have been added
            for coord in set(ds_coords)-set(self._spatial_coords(strict=False)):
                if coord not in self._added_coords():
                    raise Warning(f"Coordinate '{coord}' has been provided, but has not been added ({self._added_coords()})! Missing a decorator?")


        native_x, native_y, xvec, yvec = will_grid_be_spherical_or_cartesian(x, y, lon, lat)
        self.x_str = native_x
        self.y_str = native_y

        ds = self._init_ds(x=xvec,y=yvec, **kwargs)
        check_consistency()

        return ds

    def _init_ds(self) -> xr.Dataset:
        """Return a Dataset with only the grid coordinates and time.
        """
        raise NotImplementedError("This method depends on the specific gridding and needs to be defined in the subclass (e.g. GriddedSkeleton or PointSkeleton)")

    def ds(self):
        """Resturns the Dataset (None if doesn't exist)."""
        if not hasattr(self, 'data'):
            raise Warnign('No Dataset found. Returning None.')
            return None
        return self.data

    def _set(self, data: np.ndarray, data_name: str, only_grid_coords=False) -> None:
        self._merge_in_ds(self._compile_to_ds(data, data_name, only_grid_coords))

    def _get(self, data_name: str, default_data=None):
        """Gets data from Dataset"""
        ds = self.ds()
        if ds is None:
            return None

        data = ds.get(data_name, default_data)

        if isinstance(data, xr.DataArray):
            data = data.values.copy()

        return data

    def _merge_in_ds(self, ds_list: list[xr.Dataset]):
        """Merge in Datasets with some data into the existing Dataset of the
        Skeleton.
        """
        if not isinstance(ds_list, list):
            ds_list = [ds_list]
        for ds in ds_list:
            self.data = ds.merge(self.data, compat='override')


    def _compile_to_ds(self, data: np.ndarray, data_name:str, only_grid_coords):
        """This is used to compile a Dataset containing the given data using the
        coordinates of the Skeleton.
        """
        def check_coord_consistency():
            for i, item in enumerate(coords_dict.items()):
                if i > len(data.shape)-1:
                    raise Exception(f'{item[0]} coordinate is {len(item[1])} long, but that dimension doesnt exist in the data!!!')
                if len(item[1]) != data.shape[i]:
                    raise Exception(f'{item[0]} coordinate is {len(item[1])} long, but size of data in that dimension (dim {i}) is {data.shape[i]}!!!')

            if i < len(data.shape)-1:
                raise Warning(f'The data had {len(data.shape)} dimensions but only {i} dimensions have been defined. Missing a decorator?')


        if only_grid_coords:
            coords_dict = self._grid_coords_dict()
        else:
            coords_dict = self._coords_dict()
        check_coord_consistency()

        # Data variables
        vars_dict= {}
        vars_dict[data_name] = (coords_dict.keys(),data)

        return xr.Dataset(data_vars=vars_dict, coords=coords_dict)

    def _vars(self) -> list[str]:
        """Returns a list of the variables in the Dataset."""
        if hasattr(self, 'data'):
            return list(self.data.keys())
        return []

    def _vars_dict(self) -> list[str]:
        """Returns a dict of the variables in the Dataset."""
        return self._keys_to_dict(self._vars())

    def _all_coords(self) -> list[str]:
        """Returns a list of the coordinates in the Dataset."""
        if hasattr(self, 'data'):
            return list(self.data.coords)
        return []

    def _spatial_coords(self, strict: bool):
        """Returns a list of all spatial coords explicitly known to the
        Skeleton.

        If strict = True, then only coordinates that are found in the
        Dataset are returned.
        """

        if not strict:
            return self._spatial_coord_list

        if len(self._all_coords())>0:
            return list(set(self._all_coords()).intersection(set(self._spatial_coord_list)))

        return []

    def _added_grid_coords(self):
        """Returns a list of the grid coordinates that have been set in addition
        to the 1D-2D spatial coordinates (i.e. the skeleton coords)."""
        if hasattr(self, '_grid_coord_list'):
            added_grid_coords = self._grid_coord_list
        else:
            added_grid_coords = []
        return added_grid_coords

    def _added_gridpoint_coords(self):
        """Returns a list of the grid coordinates that have been set in addition
        to the 1D-2D spatial coordinates (i.e. the skeleton coords)."""
        if hasattr(self, '_gridpoint_coord_list'):
            added_gridpoint_coords = self._gridpoint_coord_list
        else:
            added_gridpoint_coords = []
        return added_gridpoint_coords

    def _added_coords(self):
        """Returns a list of the coordinates that have been set in addition to
        the 1D-2D spatial coordinates (i.e. the skeleton coords)."""
        return self._added_grid_coords() + self._added_gridpoint_coords()

    def  _grid_coords(self):
        """Return a list of coordinates that are in use in the Dataset to
        describe the 'outer' structure of the data (typically spatial and
        temporal dimensions).
        """

        return self._spatial_coords(strict=True) + self._added_grid_coords()

    def _gridpoint_coords(self):
        """Return a list of coordinates that are in use in the Dataset to
        describe the 'inner' structure of the data, i.e. dimensions of singe
        'grid point' (typically frequency and directions).
        """
        return self._added_gridpoint_coords()

    def _keys_to_dict(self, coords: list[str]) -> dict:
        """Takes a list of coordinates and returns a dictionary."""
        coords_dict = {}
        for coord in coords:
            coords_dict[coord] = self._get(coord)
        return coords_dict

    def _coords_dict(self):
        """Return variable dictionary of the Dataset.
        """
        return self._keys_to_dict(self._all_coords())

    def _grid_coords_dict(self):
        """Return variable dictionary of the grid coordinates of the Dataset.
        """
        return self._keys_to_dict(self._grid_coords())

    def _gridpoint_coords_dict(self):
        """Return variable dictionary of the gridpoint coordinates of the Dataset.
        """
        return self._keys_to_dict(self._gridpoint_coords())

    def _added_coords_dict(self):
        """Return variable dictionary of the added coordinates of the Dataset.
        """
        return self._keys_to_dict(self._added_coords())

    def _coords_to_size(self, coords: list[str]) -> tuple[int]:
        list = []

        for coord, val in self.ds().dims.items():
            if coord in coords:
                list.append(val)
        return tuple(list)

    def size(self) -> tuple[int]:
        """Returns the size of the Dataset."""
        return self._coords_to_size(self._all_coords())

    def spatial_size(self) -> tuple[int]:
        """Returns the size of the skeleton grid."""
        return self._coords_to_size(self._skeleton_coords(strict=True))

    def grid_size(self) -> tuple[int, int]:
        """Returns the size of the grid."""

        return self._coords_to_size(self._grid_coords())

    def gridpoint_size(self) -> tuple[int, int]:
        """Returns the inner size of one grid point."""

        return self._coords_to_size(self._gridpoint_coords())

def will_grid_be_spherical_or_cartesian(x, y, lon, lat):
    """Determines if the grid will be spherical or cartesian based on which
    inputs are given and which are None.

    Returns the ringth vector and string to identify the native values.
    """
    xy = False
    lonlat = False

    if x is not None and y is not None:
        xy = True
        native_x = 'x'
        native_y = 'y'
        xvec = x
        yvec = y

    if lon is not None and lat is not None:
        lonlat = True
        native_x = 'lon'
        native_y = 'lat'
        xvec = lon
        yvec = lat

    if xy and lonlat:
        raise ValueError("Can't set both lon/lat and x/y!")

    if not xy and not lonlat:
        raise ValueError('Have to set either lon/lat or x/y!')

    return native_x, native_y, np.array(xvec), np.array(yvec)
