import numpy as np
import pandas as pd
import xarray as xr
import utm
from copy import copy
from .coordinate_manager import CoordinateManager
from ..aux_funcs import min_distance
from .. import msg
def move_time_dim_to_front(coord_list) -> list[str]:
    if 'time' not in coord_list:
        return coord_list
    coord_list.insert(0, coord_list.pop(coord_list.index('time')))
    return coord_list

class DatasetManager:
    """Contains methods related to the creataon and handling of the Xarray
    Dataset that will be used in any object that inherits from Skeleton."""

    def __init__(self, coordinate_manager: CoordinateManager) -> None:
        self.coord_manager = coordinate_manager

    def create_structure(self, x: np.ndarray, y: np.ndarray,
                        x_str: str, y_str: str,
                        **kwargs) -> xr.Dataset:
        """Create a Dataset containing only the relevant coordinates.

        x_str, y_str = 'x', 'y' means x, y are cartesiant coordinates
        x_str, y_str = 'lon', 'lat' means x, y are spherical coordinates

        **kwargs contains any additional coordinates (e.g. time)
        """

        def check_consistency() -> None:
            """Checks that the provided coordinates are consistent with the
            coordinates that the Skeleton is defined over."""
            ds_coords = list(coord_dict.keys())
            # Check spatial coordinates
            xy_set = 'x' in ds_coords and 'y' in ds_coords
            lonlat_set = 'lon' in ds_coords and 'lat' in ds_coords
            inds_set = 'inds' in ds_coords
            if inds_set:
                ind_len = len(coord_dict['inds'])
                for key, value in var_dict.items():
                    if len(value[1]) != ind_len:
                        raise ValueError(f"Variable {key} is {len(value[1])} long but the index variable is {ind_len} long!")
            if not (xy_set or lonlat_set or inds_set):
                raise ValueError("A proper spatial grid is not set: Requires 'x' and 'y', 'lon' and 'lat' or 'inds'!")
            if sum([xy_set, lonlat_set, inds_set]) > 1:
                raise ValueError("A well defined spatial grid is not set: Requires 'x' and 'y', 'lon' and 'lat' or 'inds'!")

            # Check that all added coordinates are provided
            for coord in self.coord_manager.added_coords('all'):
                if coord not in ds_coords:
                    raise ValueError(f"Coordinate '{coord}' has been added (by a decorator?) but it was not provided when the Dataset ({ds_coords}) was created!")

            # Check that all provided coordinates have been added
            for coord in set(ds_coords)-set(self.coord_manager.spatial_coords):
                if coord not in self.coord_manager.added_coords('all'):
                    raise Warning(f"Coordinate '{coord}' has been provided, but has not been added ({self.coord_manager.added_coords('all')})! Missing a decorator?")

        def determine_coords() -> dict:
            """Creates dictonary of the coordinates"""
            coord_dict = {}
            if 'y' in self.coord_manager.initial_coords():
                coord_dict[y_str] = y
            if 'x' in self.coord_manager.initial_coords():
                coord_dict[x_str] = x
            if 'inds' in self.coord_manager.initial_coords():
                coord_dict['inds'] = np.arange(len(x))

            # Add in other possible coordinates that are set at initialization
            for key, value in kwargs.items():
                coord_dict[key] = value

            coord_dict = {c: coord_dict[c] for c in move_time_dim_to_front(list(coord_dict))}

            return coord_dict

        def determine_vars() -> dict:
            """Creates dictionary of variables"""
            var_dict = {}
            initial_vars = self.coord_manager.initial_vars()
            if 'y' in initial_vars.keys():
                if initial_vars['y'] not in coord_dict.keys():
                    raise ValueError(f"Trying to make variable 'y' depend on {initial_vars['y']}, but {initial_vars['y']} is not set as a coordinate!")
                var_dict[y_str] = ([initial_vars['y']], y)
            if 'x' in initial_vars.keys():
                if initial_vars['x'] not in coord_dict.keys():
                    raise ValueError(f"Trying to make variable 'x' depend on {initial_vars['x']}, but {initial_vars['x']} is not set as a coordinate!")
                var_dict[x_str] = ([initial_vars['x']], x)

            return var_dict


        indexed = 'inds' in self.coord_manager.initial_coords()
        x, y = clean_coordinate_vectors(x, y, is_cartesian=(x_str=='x'), indexed=indexed)
        coord_dict = determine_coords()
        var_dict = determine_vars()

        check_consistency()
        self.set_new_ds(xr.Dataset(coords=coord_dict, data_vars=var_dict))

    def set_new_ds(self, ds: xr.Dataset) -> None:
        self.data = ds

    def ds(self):
        """Resturns the Dataset (None if doesn't exist)."""
        if not hasattr(self, 'data'):
            return None
        return self.data

    def set(self, data: np.ndarray, data_name: str, coord_type: str='all') -> None:
        """Adds in new data to the Dataset.

        coord_type = 'all', 'spatial', 'grid' or 'gridpoint'
        """
        self._merge_in_ds(self.compile_to_ds(data, data_name, coord_type))

    def get(self, data_name: str, default_data=None, method='exact', **kwargs) -> xr.DataArray:
        """Gets data from Dataset.

        **kwargs can be used for slicing data.

        method = 'exact' slices
        method = 'nearest' gets nearest lon/lat point instead
        """
        ds = self.ds()
        if ds is None:
            return None

        data = ds.get(data_name, default_data)
        if isinstance(data, xr.DataArray):
            if method == 'nearest':
                if kwargs.get('lon') is None or kwargs.get('lat') is None:
                    raise Exception("Define both lon and lat when using 'nearest'")
                else:
                    __, ind = min_distance(kwargs['lon'], kwargs['lat'], ds.get('lon'), ds.get('lat'))
                    if 'inds' in ds.coords:
                        kwargs['inds'] = ind
                    else:
                        kwargs['lon'] = ds.get('lon')[ind]
                        kwargs['lat'] = ds.get('lat')[ind]

            data = self._slice_data(data, **kwargs)

        return data

    def set_attrs(self, attributes: dict, da_name: str=None) -> None:
        """Sets attributes to DataArray da_name.

        If da_name is not given, sets global attributes
        """
        if da_name is None:
            self.data.attrs = attributes
        else:
            self.data.get(da_name).attrs = attributes

    def _slice_data(self, data, **kwargs) -> xr.DataArray:
        for key, value in kwargs.items():
            if key in list(data.coords):
                data = eval(f'data.sel({key}={value})')
        return data

    def _merge_in_ds(self, ds_list: list[xr.Dataset]) -> None:
        """Merge in Datasets with some data into the existing Dataset of the
        Skeleton.
        """
        if not isinstance(ds_list, list):
            ds_list = [ds_list]
        for ds in ds_list:
            self.set_new_ds(ds.merge(self.ds(), compat='override'))


    def compile_to_ds(self, data: np.ndarray, data_name:str, coord_type: str) -> xr.Dataset:
        """This is used to compile a Dataset containing the given data using the
        coordinates of the Skeleton.

        coord_type determines over which coordinates to set the mask:

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
                        raise Exception(f'{item[0]} coordinate is {len(item[1])} long, but size of data in that dimension (dim {i}) is {data.shape[i]}!!!')

            if i < len(data.shape)-1:
                raise Exception(f'The data had {len(data.shape)} dimensions but only {i} dimensions have been defined. Missing a decorator?')

        coords_dict = self.coords_dict(coord_type)
        check_coord_consistency()

        # Data variables
        vars_dict= {}
        vars_dict[data_name] = (coords_dict.keys(),data)

        ds = xr.Dataset(data_vars=vars_dict, coords=coords_dict)
        return ds

    def vars(self) -> list[str]:
        """Returns a list of the variables in the Dataset."""
        if hasattr(self, 'data'):
            return list(self.data.keys())
        return []

    def vars_dict(self) -> list[str]:
        """Returns a dict of the variables in the Dataset."""
        return self.keys_to_dict(self.vars())

    def coords(self, type: str='all') -> list[str]:
        """Returns a list of the coordinates from the Dataset.

        'all': all coordinates in the Dataset
        'spatial': Dataset coordinates from the Skeleton (x, y, lon, lat, inds)
        'grid': coordinates for the grid (e.g. z, time)
        'gridpoint': coordinates for a grid point (e.g. frequency, direcion or time)
        """
        def list_intersection(list1, list2):
            """Uning intersections of sets doesn't necessarily preserve order"""
            list3 = []
            for val in list1:
                if val in list2:
                    list3.append(val)
            return list3

        if type not in ['all', 'spatial', 'grid', 'gridpoint']:
            raise ValueError("Type needs to be 'all', 'spatial', 'grid' or 'gridpoint'.")

        if not hasattr(self, 'data'):
            return []

        all_coords = list(self.ds().coords)
        spatial_coords = self.coord_manager.spatial_coords

        if type == 'all':
            return all_coords
        if type == 'spatial':
            return list_intersection(all_coords, spatial_coords)
        if type == 'grid':
            return move_time_dim_to_front(self.coords('spatial') + self.coord_manager.added_coords('grid'))
        if type == 'gridpoint':
            return self.coord_manager.added_coords('gridpoint')


    def keys_to_dict(self, coords: list[str]) -> dict:
        """Takes a list of coordinates and returns a dictionary."""
        coords_dict = {}
        for coord in coords:
            coords_dict[coord] = self.get(coord)
        return coords_dict

    def coords_dict(self, type: str='all') -> dict:
        """Return variable dictionary of the Dataset.

        'all': all coordinates in the Dataset
        'spatial': Dataset coordinates from the Skeleton (x, y, lon, lat, inds)
        'grid': coordinates for the grid (e.g. z, time)
        'gridpoint': coordinates for a grid point (e.g. frequency, direcion or time)
        """
        return self.keys_to_dict(self.coords(type))

    def coords_to_size(self, coords: list[str], **kwargs) -> tuple[int]:
        list = []
        data = self._slice_data(self.ds(), **kwargs)
        for coord in coords:
                list.append(len(data.get(coord)))
        # for coord, val in data.dims.items():
        #     if coord in coords:
        #         list.append(val)
        return tuple(list)


def clean_coordinate_vectors(x, y, is_cartesian, indexed):
    """Cleans up the coordinate vectors to make sure they are numpy arrays and
    have the right dimensions in case of single points etc.
    """
    def clean_lons(lon):
        mask = lon<-180
        lon[mask] = lon[mask] + 360
        mask = lon>180
        lon[mask] = lon[mask] - 360
        return lon

    # def utm_lat_mask(lat):
    #     return np.logical_and(lat < 84, lat > -80)

    x = np.array(x)
    y = np.array(y)

    if not x.shape:
        x = np.array([x])

    if not y.shape:
        y = np.array([y])

    if not is_cartesian:
        # force lon to be -180, 180
        x = clean_lons(x)

    if not indexed:
        if len(np.unique(x)) == 1 and len(x) == 2: # e.g. lon=(4.0, 4.0) should behave like lon=4.0
            x = np.unique(x)

        if len(np.unique(y)) == 1 and len(y) == 2:
            y = np.unique(y)

    return x, y
