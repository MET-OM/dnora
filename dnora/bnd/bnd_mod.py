from abc import ABC, abstractmethod
import xarray as xr
import numpy as np
from copy import copy
from typing import List
import sys
import matplotlib.pyplot as plt
from .. import msg
from ..aux import distance_2points, day_list

#from .bnd_abc import BoundaryReader, PointPicker, SpectralProcessor
#from .bnd import pick_Trivial, process_Multiply

from ..grd.grd_mod import Grid # Grid object





class BoundaryReader(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def get_coordinates(self, start_time):
        pass

    @abstractmethod
    def __call__(self, start_time, end_time, inds):
        pass

    def get_time_limits_day(self, ind):
        """Determines star and end time for the day. First and last day doesn't start at 00:00 or end at 23:59"""

        days = day_list(start_time = self.start_time, end_time = self.end_time)

        if ind == 0:
            t0 = self.start_time
            t1 = days[0].strftime('%Y-%m-%d') + "T23:59:59"
        elif ind == (len(days)-1):
            t0 = days[-1].strftime('%Y-%m-%d') + "T00:00:00"
            t1 = self.end_time
        else:
            t0 = days[ind].strftime('%Y-%m-%d') + "T00:00:00"
            t1 = days[ind].strftime('%Y-%m-%d') + "T23:59:59"
        return t0, t1

    def __str__(self):
        return (f"{self.start_time} - {self.end_time}")

class PointPicker(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, grid: Grid, bnd_lon, bnd_lat):
        return

class TrivialPicker(PointPicker):
    def __init__(self):
        pass

    def __call__(self, grid: Grid, bnd_lon, bnd_lat):
        inds = np.array(range(len(bnd_lon)))
        return inds



class SpectralProcessor(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, spec, freqs, dirs):
        pass

    @abstractmethod
    def __str__(self):
        """Describes how the spectral values as processed"""
        pass


class Multiply(SpectralProcessor):
    def __init__(self, calib_spec = 1):
        self.calib_spec = calib_spec
        return

    def __call__(self, spec, dirs, freq):
        new_spec = copy(spec)*self.calib_spec
        return new_spec, dirs, freq

    def __str__(self):
        return(f"Multiplying spectral values with {self.calib_spec}")

class BoundaryWriter(ABC):
    bnd_in: list
    bnd_points: np.array
    message: str
    @abstractmethod
    def __call__(self, bnd_out):
        pass



class Boundary:
    def __init__(self, grid: Grid, name: str = "AnonymousBoundary"):
        self.grid = copy(grid)
        self.name = copy(name)
        return

    def import_boundary(self, start_time: str, end_time: str, boundary_fetcher: BoundaryReader,  point_picker: PointPicker = TrivialPicker()):
        self.start_time = copy(start_time)
        self.end_time = copy(end_time)

        msg.header(f"{type(boundary_fetcher).__name__}: Reading coordinats of spectra...")
        lon_all, lat_all = boundary_fetcher.get_coordinates(self.start_time)


        msg.header(f"Choosing spectra with {type(point_picker).__name__}")
        inds = point_picker(self.grid, lon_all, lat_all)

        msg.header(f"{type(boundary_fetcher).__name__}: Loading boundary spectra...")
        time, freq, dirs, spec, lon, lat, source = boundary_fetcher(self.start_time, end_time, inds)

        self.data = self.compile_to_xr(time, freq, dirs, spec, lon, lat, source)
        self.mask = [True]*len(self.x())

        return


    # def process_spectra(self, spectral_processors: SpectralProcessor = Multiply(calib_spec = 1)):
    #
    #     if not isinstance(spectral_processors, list):
    #         spectral_processors = [spectral_processors]
    #
    #     for n in range (len(spectral_processors)):
    #         spectral_processor = spectral_processors[n]
    #
    #         msg.process(f"Processing spectra with {type(spectral_processor).__name__}")
    #         new_spec, new_mask, new_freq, new_dirs = spectral_processor(self.spec(), self.freq(), self.dirs(), self.time(), self.x(), self.lon(), self.lat(), self.mask)
    #
    #         self.data.spec.values = new_spec
    #         self.mask = new_mask
    #
    #         self.data = self.data.assign_coords(dirs=new_dirs)
    #         self.data = self.data.assign_coords(freq=new_freq)
    #
    #     return

    def process_spectra(self, spectral_processors: List[SpectralProcessor] = [Multiply(calib_spec = 1)]):

        if not isinstance(spectral_processors, list):
            spectral_processors = [spectral_processors]

        for n in range (len(spectral_processors)):
            spectral_processor = spectral_processors[n]

            msg.process(f"Processing spectra with {type(spectral_processor).__name__}")
            print(spectral_processor)
            for x in range(len(self.x())):
                for t in range(len(self.time())):
                    new_spec, new_dirs, new_freq = spectral_processor(self.spec()[t,x,:,:], self.dirs(), self.freq())

                    self.data.spec.values[t,x,:,:] = new_spec

            self.data = self.data.assign_coords(dirs=new_dirs)
            self.data = self.data.assign_coords(freq=new_freq)

        return


    def compile_to_xr(self, time, freq, dirs, spec, lon, lat, source):
        x = np.array(range(spec.shape[1]))
        data = xr.Dataset(
            data_vars=dict(
                spec=(["time", "x", "freq", "dirs"], spec),
            ),
            coords=dict(
                freq=freq,
                dirs=dirs,
                x=x,
                lon=(["x"], lon),
                lat=(["x"], lat),
                time=time,
            ),
            attrs=dict(source=source,
                name=self.name
            ),
            )
        return data

    def slice_data(self, start_time: str = '', end_time: str = '', x: List[int] = []):
        if isinstance(x, int):
            x = [x]
        elif not x:
            x=self.x()

        if not start_time:
            # This is not a string, but slicing works also with this input
            start_time = self.time()[0]

        if not end_time:
            # This is not a string, but slicing works also with this input
            end_time = self.time()[-1]

        sliced_data = self.data.sel(time=slice(start_time, end_time), x = x)

        return sliced_data

    def spec(self, start_time: str = '', end_time: str = '', x = []):
        spec = self.slice_data(start_time, end_time, x).spec.values

        return spec


    def time(self):
        return copy(self.data.time.values)

    def freq(self):
        return copy(self.data.freq.values)

    def dirs(self):
        return copy(self.data.dirs.values)

    def lon(self):
        return copy(self.data.lon.values)

    def lat(self):
        return copy(self.data.lat.values)

    def x(self):
        return copy(self.data.x.values)

    def days(self):
        """Determins a Pandas data range of all the days in the time span."""
        days = day_list(start_time = self.start_time, end_time = self.end_time)
        return days

    def times_in_day(self, day):
        """Determines time stamps of one given day."""
        t0 = day.strftime('%Y-%m-%d') + "T00:00:00"
        t1 = day.strftime('%Y-%m-%d') + "T23:59:59"

        times = self.slice_data(start_time = t0, end_time = t1, x = 0).time.values
        return times
