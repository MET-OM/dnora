from abc import ABC, abstractmethod
from copy import copy
import numpy as np
import xarray as xr
from typing import List
import pandas as pd
# Import objects
from ..grd.grd_mod import Grid

# Import abstract classes and needed instances of them
from ..bnd.pick import PointPicker, TrivialPicker
from .read import SpectralReader
from .. import msg
class Spectra:
    def __init__(self, grid: Grid, name: str="AnonymousSpectra"):
        self.grid = copy(grid)
        self._name = copy(name)
        self._convention = None
        self._history = []

    def import_spectra(self, start_time: str, end_time: str, spectral_reader: SpectralReader,  point_picker: PointPicker = TrivialPicker()) -> None:
        """Imports omnidirectional spectra from a certain source.

        Spectra are import between start_time and end_time from the source
        defined in the boundary_reader. Which spectra to choose spatially
        are determined by the point_picker.
        """

        self.start_time = copy(start_time)
        self.end_time = copy(end_time)
        self._history.append(copy(spectral_reader))

        msg.header(spectral_reader, "Reading coordinates of spectra...")
        lon_all, lat_all = spectral_reader.get_coordinates(self.start_time)

        msg.header(point_picker, "Choosing spectra...")
        inds = point_picker(self.grid, lon_all, lat_all)

        msg.header(spectral_reader, "Loading omnidirectional spectra...")
        time, freq, spec, mdir, spr, lon, lat, source = spectral_reader(self.start_time, self.end_time, inds)

        self.data = self.compile_to_xr(time, freq, spec, mdir, spr, lon, lat, source)
        self.mask = [True]*len(self.x())

        # E.g. are the spectra oceanic convention etc.
        self._convention = spectral_reader.convention()

    def compile_to_xr(self, time, freq, spec, mdir, spr, lon, lat, source):
        """Data from .import_spectra() is stored as an xarray Dataset."""
        x = np.array(range(spec.shape[1]))
        data = xr.Dataset(
            data_vars=dict(
                spec=(["time", "x", "freq"], spec),
                mdir=(["time", "x", "freq"], mdir),
                spr=(["time", "x", "freq"], spr),
            ),
            coords=dict(
                freq=freq,
                x=x,
                lon=(["x"], lon),
                lat=(["x"], lat),
                time=time,
            ),
            attrs=dict(source=source,
                name=self.name()
            ),
            )
        return data

    def slice_data(self, start_time: str='', end_time: str='', x: List[int]=None):
        """Slice data in space (x) and time. Returns an xarray dataset."""

        if x is None:
            x=self.x()

        if not start_time:
            # This is not a string, but slicing works also with this input
            start_time = self.time()[0]

        if not end_time:
            # This is not a string, but slicing works also with this input
            end_time = self.time()[-1]

        sliced_data = self.data.sel(time=slice(start_time, end_time), x = x)

        return sliced_data

    def spec(self, start_time: str='', end_time: str='', x: List[int]=None):
        """Slice spectra in space (x) and time. Returns an numpy array."""

        spec = self.slice_data(start_time, end_time, x).spec.values
        return spec

    def mdir(self, start_time: str='', end_time: str='', x: List[int]=None):
        """Slice mean direction in space (x) and time. Returns an numpy array."""

        return self.slice_data(start_time, end_time, x).mdir.values

    def spr(self, start_time: str='', end_time: str='', x: List[int]=None):
        """Slice mean direction in space (x) and time. Returns an numpy array."""

        return self.slice_data(start_time, end_time, x).spr.values

    def time(self):
        return pd.to_datetime(self.data.time.values)

    def dt(self) -> float:
        """ Returns time step of boundary spectra in hours."""
        return self.time().to_series().diff().dt.total_seconds().values[-1]/3600

    def freq(self):
        if hasattr(self, 'data'):
            return copy(self.data.freq.values)
        else:
            return np.array([])

    def lon(self):
        if hasattr(self, 'data'):
            return copy(self.data.lon.values)
        else:
            return np.array([])

    def lat(self):
        if hasattr(self, 'data'):
            return copy(self.data.lat.values)
        else:
            return np.array([])

    def x(self):
        if hasattr(self, 'data'):
            return copy(self.data.x.values)
        else:
            return np.array([])

    def days(self):
        """Determins a Pandas data range of all the days in the time span."""
        days = day_list(start_time = self.start_time, end_time = self.end_time)
        return days

    def size(self):
        return (len(self.time()), len(self.x()))


    def name(self) -> str:
        """Return the name of the grid (set at initialization)."""
        return copy(self._name)

    def convention(self) -> str:
        """Returns the convention (WW3/Ocean/Met/Math) of the spectra"""
        if hasattr(self, '_convention'):
            return copy(self._convention)
        else:
            return None

    def times_in_day(self, day):
        """Determines time stamps of one given day."""
        t0 = day.strftime('%Y-%m-%d') + "T00:00:00"
        t1 = day.strftime('%Y-%m-%d') + "T23:59:59"

        times = self.slice_data(start_time=t0, end_time=t1, x=[0]).time.values
        return times

    def __str__(self) -> str:
        """Prints status of spectra."""

        msg.header(self, f"Status of spectra {self.name()}")
        msg.plain(f"Contains data ({len(self.x())} points) for {self.start_time} - {self.end_time}")
        msg.plain(f"Data covers: lon: {min(self.lon())} - {max(self.lon())}, lat: {min(self.lat())} - {max(self.lat())}")
        if len(self._history) > 0:
            msg.blank()
            msg.plain("Object has the following history:")
            for obj in self._history:
                msg.process(f"{obj.__class__.__bases__[0].__name__}: {type(obj).__name__}")
        #msg.print_line()
        #msg.plain("The Boundary is for the following Grid:")
        #print(self.grid)

        msg.print_line()

        return ''
