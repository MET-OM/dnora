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
from ..skeletons.gridded_skeleton import GriddedSkeleton
from ..skeletons.coordinate_factory import add_time
from ..skeletons.datavar_factory import add_datavar

@add_datavar(name='v', default_value=0.)
@add_datavar(name='u', default_value=0.)
@add_time(grid_coord=True)
class Forcing(GriddedSkeleton):
    def __init__(self, grid, name='AnonymousForcing'):
        self.grid = copy(grid)
        self._name = copy(name)
        self._history = []

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
            forcing_reader = DnoraNc(files=glob.glob(f'{cache_folder}/{cache_name}*'))

        msg.header(forcing_reader, "Loading wind forcing...")
        time, u, v, lon, lat, x, y, attributes = forcing_reader(
            self.grid, start_time, end_time, expansion_factor)

        self._init_structure(x, y, lon, lat, time=time)
        self.ds_manager.set(u, 'u', coord_type='all')
        self.ds_manager.set(v, 'v', coord_type='all')

        self.ds_manager.set_attrs(attributes)
        ### Patch data if read from cache and all data not found
        if read_cache and not cache_empty:
            patch_start, patch_end = aux_funcs.determine_patch_periods(self.time(), start_time, end_time)
            if patch_start:
                msg.info('Not all data found in cache. Patching from original source...')
                wnd_list = [self.data]
                for t0, t1 in zip(patch_start, patch_end):
                    wnd_ds = original_forcing_reader(self.grid, t0, t1, expansion_factor)
                    wnd_list.append(wnd_ds)

                self.data = xr.concat(wnd_list, dim="time").sortby('time')

        if write_cache:
            msg.info('Caching data:')
            for month in self.months():
                cache_file = f"{cache_name}_{month.strftime('%Y-%m')}.nc"
                t0 = f"{month.strftime('%Y-%m-01')}"
                d1 = monthrange(int(month.strftime('%Y')), int(month.strftime('%m')))[1]
                t1 = f"{month.strftime(f'%Y-%m-{d1}')}"
                outfile = f'{cache_folder}/{cache_file}'
                if os.path.exists(outfile):
                    os.remove(outfile)
                self.data.sel(time=slice(t0, t1)).to_netcdf(outfile)
                msg.to_file(outfile)

        return

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
