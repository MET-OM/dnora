from functools import partial
from .. import aux_funcs
from datetime import datetime
import pandas as pd
import numpy as np
from .coordinate_manager import CoordinateManager

def coord_decorator(coord_name, grid_coord, c, stash_get=False):
    """stash_get = True means that the coordinate data can be accessed
    by method ._{coord_name}() instead of .{coord_name}()

    This allows for alternative definitions of the get-method elsewere."""
    def get_coord(self, data_array=False, **kwargs):
        if not self._structure_initialized():
            return None
        data = self.ds_manager.get(coord_name, **kwargs)
        if data_array:
            return data
        return data.values.copy()

    if not hasattr(c, '_coord_manager'):
        c._coord_manager =  CoordinateManager()

    c._coord_manager.add_coord(coord_name, grid_coord)
    if stash_get:
        exec(f'c._{coord_name} = get_coord')
    else:
        exec(f'c.{coord_name} = get_coord')
    return c

def add_time(grid_coord: bool=False, coord_name: str='time'):
    def wrapper(c):
        def unique_times(times, strf: str):
            return np.unique(np.array(pd.to_datetime(times).strftime(strf).to_list()))


        def hours(self, datetime=True):
            """Determins a Pandas data range of all the days in the time span."""
            if not self._structure_initialized():
                return None
            times = self.ds_manager.get(coord_name).values.copy()
            if datetime:
                return pd.to_datetime(unique_times(times, '%Y-%m-%d %H'))
            else:
                return list(unique_times(times, '%Y-%m-%d %H'))

        def days(self, datetime=True):
            """Determins a Pandas data range of all the days in the time span."""
            if not self._structure_initialized():
                return None
            times = self.ds_manager.get(coord_name).values.copy()
            if datetime:
                return pd.to_datetime(unique_times(times, '%Y-%m-%d'))
            else:
                return list(unique_times(times, '%Y-%m-%d'))

        def months(self, datetime=True):
            """Determins a Pandas data range of all the months in the time span."""
            if not self._structure_initialized():
                return None
            times = self.ds_manager.get(coord_name).values.copy()
            if datetime:
                return pd.to_datetime(unique_times(times, '%Y-%m'))
            else:
                return list(unique_times(times, '%Y-%m'))

        def years(self, datetime=True):
            """Determins a Pandas data range of all the months in the time span."""
            if not self._structure_initialized():
                return None
            times = self.ds_manager.get(coord_name).values.copy()
            if datetime:
                return pd.to_datetime(unique_times(times, '%Y'))
            else:
                return list(unique_times(times, '%Y'))

        def get_time(self, data_array=False, **kwargs):
            if not self._structure_initialized():
                return (None, None)
            data = self.ds_manager.get(coord_name, **kwargs)
            if data_array:
                return data
            return pd.to_datetime(data.values.copy())

        if not hasattr(c, '_coord_manager'):
            c._coord_manager =  CoordinateManager()

        c._coord_manager.add_coord(coord_name, grid_coord)
        exec(f'c.{coord_name} = get_time')

        c.days = days
        c.months = months
        c.years = years
        return c

    return wrapper

def add_frequency(grid_coord: bool=False, coord_name: str='freq'):
    def wrapper(c):
        def get_freq(self, angular=False):
            if not self._structure_initialized():
                return None
            freq = self.ds_manager.get(coord_name).values.copy()
            if angular:
                freq = 2*np.pi*freq
            return freq

        def df(self, angular=False):
            if not self._structure_initialized():
                return None
            freq = get_freq(self, angular=angular).values.copy()
            return (freq[-1]-freq[0])/(len(freq)-1)

        if not hasattr(c, '_coord_manager'):
            c._coord_manager =  CoordinateManager()
        c._coord_manager.add_coord(coord_name, grid_coord)
        exec(f'c.{coord_name} = get_freq')
        c.df = df



        return c

    return wrapper

def add_direction(grid_coord: bool=False, coord_name: str='dirs'):
    def wrapper(c):
        def get_dirs(self, radians=False):
            if not self._structure_initialized():
                return None
            dirs = self.ds_manager.get(coord_name).values.copy()
            if radians:
                dirs = dirs*np.pi/180
            return dirs

        def ddir(self, radians=False):
            if not self._structure_initialized():
                return None
            dirs = get_dirs(self, radians=False).values.copy()
            dmax = 2*np.pi if radians else 360
            return dmax/len(dirs)

        if not hasattr(c, '_coord_manager'):
            c._coord_manager =  CoordinateManager()
        c._coord_manager.add_coord(coord_name, grid_coord)
        exec(f'c.{coord_name} = get_dirs')
        c.dd = ddir
        return c

    return wrapper


def add_dummy(grid_coord: bool=False, coord_name: str='dummy'):
    return partial(coord_decorator, coord_name, grid_coord)
