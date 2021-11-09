from abc import ABC, abstractmethod
import xarray as xr
import numpy as np
from copy import copy
import sys
import matplotlib.pyplot as plt
from .. import msg
from ..aux import distance_2points, day_list

class ForcingReader(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, start_time, end_time, inds):
        pass

    def __str__(self):
        return (f"{self.start_time} - {self.end_time}")


class ForcingWriter(ABC):
    @abstractmethod
    def __call__(self, forcing_out):
        pass


class Forcing:
    def __init__(self, grid, name='AnonymousForcing'):
        self.grid = copy(grid)
        self.name = name

    def import_forcing(self, start_time: str, end_time: str, forcing_fetcher: ForcingReader, expansion_factor=1.2):
        self.start_time = copy(start_time)
        self.end_time = copy(end_time)

        msg.header(
            f"{type(forcing_fetcher).__name__}: Loading wind forcing...")
        self.data = forcing_fetcher(
            self.grid, start_time, end_time, expansion_factor)

        return

    def days(self):
        """Determins a Pandas data range of all the days in the time span."""
        days = day_list(start_time=self.start_time, end_time=self.end_time)
        return days

    def time(self):
        return copy(self.data.time.values)

    def u(self):
        return copy(self.data.u.values)

    def v(self):
        return copy(self.data.v.values)

    def nx(self):
        return (self.data.x_wind_10m.shape[0])

    def ny(self):
        return (self.data.x_wind_10m.shape[1])

    def slice_data(self, start_time: str = '', end_time: str = ''):
        if not start_time:
            # This is not a string, but slicing works also with this input
            start_time = self.time()[0]

        if not end_time:
            # This is not a string, but slicing works also with this input
            end_time = self.time()[-1]

        sliced_data = self.data.sel(time=slice(start_time, end_time))

        return sliced_data

    def times_in_day(self, day):
        """Determines time stamps of one given day."""
        t0 = day.strftime('%Y-%m-%d') + "T00:00:00"
        t1 = day.strftime('%Y-%m-%d') + "T23:59:59"

        times = self.slice_data(start_time=t0, end_time=t1).time.values
        return times
