import xarray as xr
import numpy as np
from copy import copy
import pandas as pd
import sys
import re
import glob, os
from calendar import monthrange
# Import abstract classes and needed instances of them
from .read import ForcingReader, DnoraNc
from .write import ForcingWriter

# Import default values and aux_funcsiliry functions
from .. import msg
from .. import aux_funcs
class Forcing:
    def __init__(self, grid, name='AnonymousForcing'):
        self.grid = copy(grid)
        self._name = copy(name)
        self._history = []
        return

    def import_forcing(self, start_time: str, end_time: str,
                    forcing_reader: ForcingReader,
                    expansion_factor: float=1.2,
                    write_cache: bool=False,
                    read_cache: bool=False,
                    cache_name: str='#Grid_#Lon0_#Lon1_#Lat0_#Lat1'):
        """Imports forcing data from a certain source.

        Data are import between start_time and end_time from the source
        defined in the forcing_reader. Data are read around an area defined
        by the Grid object passed at initialization of this object.
        """

        self.start_time = copy(start_time)
        self.end_time = copy(end_time)
        self._history.append(copy(forcing_reader))

        if write_cache or read_cache:
            cache_folder, cache_name, cache_empty = aux_funcs.setup_cache('wnd', forcing_reader.name(), cache_name, self.grid)

        if read_cache and not cache_empty:
            msg.info('Reading wind forcing data from cache!!!')
            original_forcing_reader = copy(forcing_reader)
            forcing_reader = DnoraNc(files=glob.glob(f'{cache_folder}/{cache_name}_wnd.nc'))

        msg.header(forcing_reader, "Loading wind forcing...")
        self.data = forcing_reader(
            self.grid, start_time, end_time, expansion_factor)

        ### Patch data if read from cache and all data not found
        if read_cache and not cache_empty:
            patch_start, patch_end = aux_funcs.determine_patch_periods(self.time(), start_time, end_time)
            if patch_start:
                msg.info(
                    f'Timesteps stored in {cache_name}_wnd.nc does not cover the entire range. set read_cache to false and rerun.')
                exit(-1)

        if write_cache:
            msg.info('Caching data:')
            cache_file = f"{cache_name}_wnd.nc"
            outfile = f'{cache_folder}/{cache_file}'
            if os.path.exists(outfile):
                os.remove(outfile)
            self.data.to_netcdf(outfile)
            msg.to_file(outfile)

        return

    def days(self):
        """Determins a Pandas data range of all the days in the time span."""

        days = aux_funcs.day_list(start_time=self.start_time, end_time=self.end_time)
        return days

    def months(self):
        """Determins a Pandas data range of all the months in the time span."""
        if len(self.time()) == 0:
            return []
        t0 = self.time()[0].strftime('%Y-%m') + '-01'
        t1 = self.time()[-1].strftime('%Y-%m') + '-01'

        return pd.date_range(start=t0, end=t1, freq='MS')

    def name(self) -> str:
        """Return the name of the grid (set at initialization)."""

        return copy(self._name)

    def time(self):
        if hasattr(self, 'data'):
            return pd.to_datetime(self.data.time.values)
        else:
            return None

    def dt(self) -> float:
        """ Returns time step of forcing data in hours."""
        if hasattr(self, 'data'):
            return self.time().to_series().diff().dt.total_seconds().values[-1]/3600
        else:
            return 0

    def u(self):
        if hasattr(self, 'data'):
            return copy(self.data.u.values)
        else:
            return np.array([[[]]])

    def v(self):
        if hasattr(self, 'data'):
            return copy(self.data.v.values)
        else:
            return np.array([[[]]])

    def magnitude(self):
        if np.min(self.u()).shape == 0:
            return np.array([[[]]])
        if np.min(self.v()).shape == 0:
            return self.u()

        return (self.u()**2 + self.v()**2)**0.5

    def nx(self):
        if min(self.size()) > 0:
            return self.size()[2]
        else:
            return 0

    def ny(self):
        if min(self.size()) > 0:
            return (self.size()[1])
        else:
            return 0

    def nt(self):
        if min(self.size()) > 0:
            return (self.size()[0])
        else:
            return 0

    def lon(self):
        """Returns a longitude vector of the grid."""

        if hasattr(self, 'data'):
            lon = copy(self.data.lon.values)
        else:
            lon = np.array([])
        return lon

    def lat(self):
        """Returns a latitude vector of the grid."""

        if hasattr(self, 'data'):
            lat = copy(self.data.lat.values)
        else:
            lat = np.array([])
        return lat

    def size(self) -> tuple:
        """Returns the size (nt, ny, nx) of the grid."""

        return self.u().shape

    def _point_list(self, mask):
        """Provides a list on longitudes and latitudes with a given mask.

        Used to e.g. generate list of boundary points or land points.
        """

        meshlon, meshlat=np.meshgrid(self.lon(),self.lat())
        lonlat_flat = np.column_stack((meshlon.ravel(),meshlat.ravel()))
        mask_flat = mask.ravel()

        return lonlat_flat[mask_flat]


    def slice_data(self, start_time: str='', end_time: str=''):
        """Slice data in time. Returns an xarray dataset."""

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

    def __str__(self) -> str:
        """Prints status of forcing."""

        msg.header(self, f"Status of forcing {self.name()}")
        if self.time() is not None:
            msg.plain(f"Contains data for {self.time()[0]} - {self.time()[-1]}")
            msg.plain(f"\t dt={self.dt()} hours, i.e. ({self.nt()} time steps)")
            msg.plain(f"Data covers: lon: {min(self.lon())} - {max(self.lon())}, lat: {min(self.lat())} - {max(self.lat())}")
            msg.plain(f"\t {self.ny()}x{self.nx()} grid points, dlon/dlat={np.mean(np.diff(self.lon()))}/{np.mean(np.diff(self.lat()))}")
        if len(self._history) > 0:
            msg.blank()
            msg.plain("Object has the following history:")
            for obj in self._history:
                msg.process(f"{obj.__class__.__bases__[0].__name__}: {type(obj).__name__}")
        #msg.print_line()
        #msg.plain("The Forcing is for the following Grid:")
        #print(self.grid)

        msg.print_line()

        return ''
