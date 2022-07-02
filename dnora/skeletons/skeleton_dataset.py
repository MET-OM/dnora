import numpy as np
import pandas as pd
import xarray as xr
import utm
from copy import copy

class SkeletonDataset:
    """Contains methods related to the creataon and handling of the Xarray
    Dataset that will be used in any object that inherits from Skeleton."""

    _spatial_coord_list = ['x', 'y', 'lon', 'lat', 'inds']

    def _create_structure(self, x=None, y=None, lon=None, lat=None, **kwargs) -> None:
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
            for coord in self._added_coords('all'):
                if coord not in ds_coords:
                    raise ValueError(f"Coordinate '{coord}' has been added (by a decorator?) but it was not provided when the Dataset ({ds_coords}) was created!")

            # Check that all provided coordinates have been added
            for coord in set(ds_coords)-set(self._spatial_coord_list):
                if coord not in self._added_coords('all'):
                    raise Warning(f"Coordinate '{coord}' has been provided, but has not been added ({self._added_coords('all')})! Missing a decorator?")


        native_x, native_y, xvec, yvec = will_grid_be_spherical_or_cartesian(x, y, lon, lat)
        self.x_str = native_x
        self.y_str = native_y

        ds = self._init_ds(x=xvec,y=yvec, **kwargs)
        check_consistency()
        self.data = ds
        self._reset_masks()


    def _init_ds(self, x: np.ndarray, y: np.ndarray, **kwargs) -> xr.Dataset:
        """Creates a Dataset containing the spatial variables as coordinates.
        x, y are assumed to be in native_format.

        Any additional keyword arguments are also added to the coordinate list.
        """

        def determine_coords(x, y):
            coords = self._initial_coords() # Get list of coords from subclass

            coord_dict = {}
            if 'y' in coords:
                coord_dict[self.y_str] = y
            if 'x' in coords:
                coord_dict[self.x_str] = x
            if 'inds' in coords:
                coord_dict['inds'] = np.arange(len(x))

            # Add in other possible coordinates that are set at initialization
            for key, value in kwargs.items():
                coord_dict[key] = value

            return coord_dict

        def determine_vars(x, y, set_coords):
            vars = self._initial_vars() # Get dict of variables from subclass
            var_dict = {}

            if 'y' in vars.keys():
                if vars['y'] not in set_coords:
                    raise ValueError(f"Trying to make variable 'y' depend on {vars['y']}, but {vars['y']} is not set as a coordinate!")
                var_dict[self.y_str] = ([vars['y']], y)
            if 'x' in vars.keys():
                if vars['x'] not in set_coords:
                    raise ValueError(f"Trying to make variable 'x' depend on {vars['x']}, but {vars['x']} is not set as a coordinate!")
                var_dict[self.x_str] = ([vars['x']], x)

            return var_dict

        coord_dict = determine_coords(x, y)
        var_dict = determine_vars(x, y, coord_dict.keys())

        return xr.Dataset(coords=coord_dict, data_vars=var_dict)

    def _reset_masks(self) -> None:
        """Resets the mask to default values."""
        masks = getattr(self, '_mask_dict', {})

        for name, mask_tuple in masks.items():
            eval(f'self._update_{name}_mask()')

    def ds(self):
        """Resturns the Dataset (None if doesn't exist)."""
        if not hasattr(self, 'data'):
            raise Warning('No Dataset found. Returning None.')
            return None
        return self.data

    def _set(self, data: np.ndarray, data_name: str, coords: str='all') -> None:
        self._merge_in_ds(self._compile_to_ds(data, data_name, coords))

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


    def _compile_to_ds(self, data: np.ndarray, data_name:str, type: str):
        """This is used to compile a Dataset containing the given data using the
        coordinates of the Skeleton.

        'type' determines over which coordinates to set the mask:

        'all': all coordinates in the Dataset
        'spatial': Dataset coordinates from the Skeleton (x, y, lon, lat, inds)
        'grid': coordinates for the grid (e.g. z, time)
        'gridpoint': coordinates for a grid point (e.g. frequency, direcion or time)
        """
        def check_coord_consistency():
            for i, item in enumerate(coords_dict.items()):
                if i > len(data.shape)-1:
                    raise Exception(f'{item[0]} coordinate is {len(item[1])} long, but that dimension doesnt exist in the data!!!')
                if len(item[1]) != data.shape[i]:
                    breakpoint()
                    raise Exception(f'{item[0]} coordinate is {len(item[1])} long, but size of data in that dimension (dim {i}) is {data.shape[i]}!!!')

            if i < len(data.shape)-1:
                raise Warning(f'The data had {len(data.shape)} dimensions but only {i} dimensions have been defined. Missing a decorator?')


        coords_dict = self._coords_dict(type)

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

    def _added_coords(self, type: str='all') -> list[str]:
        """Returns list of coordinates that have been added to the fixed
        Skeleton coords.

        'all': All added coordinates
        'grid': coordinates for the grid (e.g. z, time)
        'gridpoint': coordinates for a grid point (e.g. frequency, direcion or time)
        """

        if type == 'all':
            return self._added_coords('grid') + self._added_coords('gridpoint')
        elif type == 'grid':
            return getattr(self, '_grid_coord_list', [])
        elif type == 'gridpoint':
            return getattr(self, '_gridpoint_coord_list', [])
        else:
            print("Type needs to be 'all', 'grid' or 'gridpoint'.")
            return None


    def coords(self, type: str='all') -> list[str]:
        """Returns a list of the coordinates from the Dataset.

        'all': all coordinates in the Dataset
        'spatial': Dataset coordinates from the Skeleton (x, y, lon, lat, inds)
        'grid': coordinates for the grid (e.g. z, time)
        'gridpoint': coordinates for a grid point (e.g. frequency, direcion or time)
        """
        if not hasattr(self, 'data'):
            return []

        all_coords = list(self.data.coords)
        spatial_coords = self._spatial_coord_list

        if type == 'all':
            return all_coords
        elif type == 'spatial':
            return list(set(all_coords).intersection(set(spatial_coords)))
        elif type == 'grid':
            return self.coords('spatial') + self._added_coords('grid')
        elif type == 'gridpoint':
            return self._added_coords('gridpoint')
        else:
            print("Type needs to be 'full', 'spatial', 'grid' or 'gridpoint'.")
            return None

    def _keys_to_dict(self, coords: list[str]) -> dict:
        """Takes a list of coordinates and returns a dictionary."""
        coords_dict = {}
        for coord in coords:
            coords_dict[coord] = self._get(coord)
        return coords_dict

    def _coords_dict(self, type: str='all') -> dict:
        """Return variable dictionary of the Dataset.

        'all': all coordinates in the Dataset
        'spatial': Dataset coordinates from the Skeleton (x, y, lon, lat, inds)
        'grid': coordinates for the grid (e.g. z, time)
        'gridpoint': coordinates for a grid point (e.g. frequency, direcion or time)
        """
        return self._keys_to_dict(self.coords(type))

    def _coords_to_size(self, coords: list[str]) -> tuple[int]:
        list = []

        for coord, val in self.ds().dims.items():
            if coord in coords:
                list.append(val)
        return tuple(list)

    def size(self, type: str='all') -> tuple[int]:
        """Returns the size of the Dataset.

        'all': size of entire Dataset
        'spatial': size over coordinates from the Skeleton (x, y, lon, lat, inds)
        'grid': size over coordinates for the grid (e.g. z, time)
        'gridpoint': size over coordinates for a grid point (e.g. frequency, direcion or time)
        """
        return self._coords_to_size(self.coords(type))


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

    if isinstance(xvec, float) or isinstance(xvec, int):
        xvec = [xvec]

    if isinstance(yvec, float) or isinstance(yvec, int):
        yvec = [yvec]

    return native_x, native_y, np.array(xvec), np.array(yvec)
