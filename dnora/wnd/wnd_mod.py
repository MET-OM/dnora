import xarray as xr
import numpy as np
from copy import copy
import pandas as pd
import sys
import re
from calendar import monthrange
# Import abstract classes and needed instances of them
from .read import ForcingReader, DnoraNc
from .write import ForcingWriter
from ..grd.grd_mod import Grid
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
        self._grid = grid
        self._name = name
        self._history = []

    def import_forcing(self, start_time: str, end_time: str,
                    forcing_reader: ForcingReader,
                    **kwargs):
        """Imports forcing data from a certain source.

        Data are import between start_time and end_time from the source
        defined in the forcing_reader. Data are read around an area defined
        by the Grid object passed at initialization of this object.
        """

        self._history.append(copy(forcing_reader))

        msg.header(forcing_reader, "Loading wind forcing...")
        time, u, v, lon, lat, x, y, attributes = forcing_reader(
            self.grid(), start_time, end_time, **kwargs)

        self._init_structure(x, y, lon, lat, time=time)
        self.ds_manager.set(u, 'u', coord_type='all')
        self.ds_manager.set(v, 'v', coord_type='all')

        self.ds_manager.set_attrs(attributes)

        return

    def grid(self) -> Grid:
        if hasattr(self, '_grid'):
            return self._grid
        return None

    def magnitude(self):
        return (self.u()**2 + self.v()**2)**0.5

    def __str__(self) -> str:
        """Prints status of forcing."""

        msg.header(self, f"Status of forcing {self.name}")
        if not self.time() == (None, None):
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
        #print(self.grid())

        msg.print_line()

        return ''
