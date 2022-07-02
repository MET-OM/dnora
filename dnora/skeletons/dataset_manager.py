import numpy as np
import pandas as pd
import xarray as xr
import utm
from copy import copy
from coordinate_manager import CoordinateManager
class DatasetManager:
    """Contains methods related to the creataon and handling of the Xarray
    Dataset that will be used in any object that inherits from Skeleton."""

    def __init__(self, coordinate_manager: CoordinateManager) -> None:
        self.coord_manager = coordinate_manager

    def create_structure(self, x: np.ndarray, y: np.ndarray,
                        x_str: str, y_str: str,
                        **kwargs) -> None:
        """Create the first Dataset. x and y are assummed to be in the native
        format (specified by x_str = 'x'/'lon' and y_str='y'/'lat')"""

        def check_consistency() -> None:
            ds_coords = list(ds.coords)
            # Check spatial coordinates
            xy_set = 'x' in ds_coords and 'y' in ds_coordsx_str
            lonlat_set = 'lon' in ds_coords and 'lat' in ds_coords
            inds_set = 'inds' in ds_coords
            if not (xy_set or lonlat_set or inds_set):
                raise ValueError("A proper spatial grid is not set: Requires 'x' and 'y', 'lon' and 'lat' or 'inds'!")
            if sum([xy_set, lonlat_set, inds_set]) > 1:
                raise ValueError("A well defined spatial grid is not set: Requires 'x' and 'y', 'lon' and 'lat' or 'inds'!")

            # Check that all added coordinates are porvided
            for coord in self.coord_manager.added_coords('all'):
                if coord not in ds_coords:
                    raise ValueError(f"Coordinate '{coord}' has been added (by a decorator?) but it was not provided when the Dataset ({ds_coords}) was created!")

            # Check that all provided coordinates have been added
            for coord in set(ds_coords)-set(self.coord_manager.spatial_coords):
                if coord not in self.coord_manager.added_coords('all'):
                    raise Warning(f"Coordinate '{coord}' has been provided, but has not been added ({self.coord_manager.added_coords('all')})! Missing a decorator?")

        def determine_coords():
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

            return coord_dict

        def determine_vars():
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

        coord_dict = determine_coords()
        var_dict = determine_vars()
        ds = xr.Dataset(coords=coord_dict, data_vars=var_dict)
        #ds = self.init_ds(x=x,y=y,)
        check_consistency()
        self.data = ds
        self.reset_masks()
        self.reset_datavars()


    def reset_masks(self) -> None:
        """Resets the mask to default values."""
        for name in self.coord_manager.added_masks():
            # update-method sets empty mask when no is provided
            self.update_mask(name)

    def reset_datavars(self) -> None:
        """Resets the data variables to default values."""
        for name in self.coord_manager.added_vars():
            # update-method sets empty mask when no is provided
            self.update_datavar(name)

    def update_mask(self, name: str, updated_mask=None) -> None:
        masks = self.coord_manager.added_masks()
        if name not in masks.keys():
            raise ValueError(f'A mask named {name} has not been defines ({list(masks.keys())})')
        coords, get_empty = masks[name]

        if updated_mask is None:
            updated_mask = get_empty(self)
        self.set(data=updated_mask, data_name=f'{name}_mask', coords=coords)

    def update_datavar(self, name: str, updated_var=None) -> None:
        vars = self.coord_manager.added_vars()
        if name not in vars.keys():
            raise ValueError(f'A data variable named {name} has not been defines ({list(vars.keys())})')
        coords, get_empty = vars[name]
        if updated_var is None:
            updated_var = get_empty(self)
        self.set(data=updated_var, data_name=name, coords=coords)

    def ds(self):
        """Resturns the Dataset (None if doesn't exist)."""
        if not hasattr(self, 'data'):
            raise Warning('No Dataset found. Returning None.')
            return None
        return self.data

    def set(self, data: np.ndarray, data_name: str, coords: str='all') -> None:
        self.merge_in_ds(self.compile_to_ds(data, data_name, coords))

    def get(self, data_name: str, default_data=None, empty=False):
        """Gets data from Dataset"""
        if empty:
            ## For 'sea_mask' check both 'sea' in the mask dict
            data_tuple = self.coord_manager.added_vars().get(data_name) or self.coord_manager.added_masks().get(data_name[0:-5])
            if data_tuple is None:
                return None
            __, get_empty = data_tuple

            return get_empty(self)

        ds = self.ds()
        if ds is None:
            return None

        data = ds.get(data_name, default_data)

        if isinstance(data, xr.DataArray):
            data = data.values.copy()

        return data

    def is_empty(self, data_name):
        """Checks if a Dataset variable is empty."""
        data = self.get(data_name)
        empty_data = self.get(data_name, empty=True)
        if data is None:
            return False
        return np.allclose(data.astype(float), empty_data.astype(float))

    def merge_in_ds(self, ds_list: list[xr.Dataset]):
        """Merge in Datasets with some data into the existing Dataset of the
        Skeleton.
        """
        if not isinstance(ds_list, list):
            ds_list = [ds_list]
        for ds in ds_list:
            if list(ds.coords)[0] == 'lon':
                breakpoint()
            self.data = ds.merge(self.data, compat='override')



    def compile_to_ds(self, data: np.ndarray, data_name:str, type: str):
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


        coords_dict = self.coords_dict(type)
        check_coord_consistency()

        if list(coords_dict.keys())[0] == 'lon':
            breakpoint()
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

        all_coords = list(self.data.coords)
        spatial_coords = self.coord_manager.spatial_coords

        if type == 'all':
            return all_coords
        elif type == 'spatial':
            #return list(set(all_coords).intersection(set(spatial_coords)))
            return list_intersection(all_coords, spatial_coords)
        elif type == 'grid':
            return self.coords('spatial') + self.coord_manager.added_coords('grid')
        elif type == 'gridpoint':
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

    def coords_to_size(self, coords: list[str]) -> tuple[int]:
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
        return self.coords_to_size(self.coords(type))
