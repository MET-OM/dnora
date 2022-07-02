from functools import partial
from dnora import aux_funcs
from datetime import datetime
import pandas as pd
from coordinate_manager import CoordinateManager
# def add_coord_name_to_list(c, coord_name, grid_coord):
#     if grid_coord:
#         if not hasattr(c, '_grid_coord_list'):
#             c._grid_coord_list = []
#         c._grid_coord_list.append(coord_name)
#     else:
#         if not hasattr(c, '_gridpoint_coord_list'):
#             c._gridpoint_coord_list = []
#         c._gridpoint_coord_list.append(coord_name)
#     return c

def coord_decorator(coord_name, grid_coord, c):
    def get_coord(self):
        return self.ds_manager.get(coord_name)

    if not hasattr(c, '_coord_manager'):
        c._coord_manager =  CoordinateManager()

    c._coord_manager.add_coord(coord_name, grid_coord)
    #c = add_coord_name_to_list(o, coord_name, grid_coord)
    exec(f'c.{coord_name} = get_coord')
    return c

def add_time(grid_coord: bool=False, coord_name: str='time'):
    def wrapper(c):
        # def days(self):
        #     """Determins a Pandas data range of all the days in the time span."""
        #     times = pd.to_datetime(self.ds_manager.get(coord_name))
        #     return aux_funcs.day_list(start_time=times[0], end_time=times[1])
        #
        # def months(self):
        #     """Determins a Pandas data range of all the months in the time span."""
        #     times = pd.to_datetime(self.ds_manager.(coord_name))
        #     return aux_funcs.month_list(start_time=times[0], end_time=times[1])
        #
        # def years(self):
        #     """Determins a Pandas data range of all the months in the time span."""
        #     times = pd.to_datetime(self.ds_manager.(coord_name))
        #     return aux_funcs.year_list(start_time=times[0], end_time=times[1])

        def get_coord(self):
            return pd.to_datetime(self.ds_manager.get(coord_name))

        if not hasattr(c, '_coord_manager'):
            c._coord_manager =  CoordinateManager()

        c._coord_manager.add_coord(coord_name, grid_coord)
        exec(f'c.{coord_name} = get_coord')

        # c.days = days
        # c.months = months
        # c.years = years
        return c

    return wrapper

def add_frequency(grid_coord: bool=False, coord_name: str='freq'):
    return partial(coord_decorator, coord_name, grid_coord)

def add_direction(grid_coord: bool=False, coord_name: str='dirs'):
    return partial(coord_decorator, coord_name, grid_coord)
