from abc import ABC, abstractmethod
import xarray as xr
import numpy as np
from copy import copy
import pandas as pd
import sys
import re
import matplotlib.pyplot as plt
from .. import msg
from ..aux import distance_2points, day_list, create_filename_obj, create_filename_time, add_file_extension
from ..defaults import dflt_frc

class ForcingReader(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, start_time, end_time, inds):
        pass

    def __str__(self):
        return (f"{self.start_time} - {self.end_time}")


class Forcing:
    def __init__(self, grid, name='AnonymousForcing'):
        self.grid = copy(grid)
        self._name = copy(name)
        return

    def import_forcing(self, start_time: str, end_time: str, forcing_reader: ForcingReader, expansion_factor: float = 1.2):
        self.start_time = copy(start_time)
        self.end_time = copy(end_time)

        msg.header(
            f"{type(forcing_reader).__name__}: Loading wind forcing...")
        self.data = forcing_reader(
            self.grid, start_time, end_time, expansion_factor)

        return

    def export_forcing(self, forcing_writer):
        output_file, output_folder = forcing_writer(self)

        # This is set as info in case an input file needs to be generated
        self._written_as = output_file
        self._written_to = output_folder
        return


    def days(self):
        """Determins a Pandas data range of all the days in the time span."""
        days = day_list(start_time=self.start_time, end_time=self.end_time)
        return days

    def name(self):
        """Return the name of the grid (set at initialization)."""
        return copy(self._name)


    def filename(self, filestring: str=dflt_frc['fs']['General'], datestring: str=dflt_frc['ds']['General'], extension: str='', defaults: str=''):
        # E.g. defaults='SWAN' uses all SWAN defaults
        if defaults:
            filestring = dflt_frc['fs'][defaults]
            datestring = dflt_frc['ds'][defaults]
            extension = dflt_frc['ext'][defaults]

        # Substitute placeholders for objects ($Grid etc.)
        filename = create_filename_obj(filestring=filestring, objects=[self, self.grid])
        # Substitute placeholders for times ($T0 etc.)
        filename = create_filename_time(filestring=filename, times=[self.start_time, self.end_time], datestring=datestring)

        # Possible clean up
        filename = re.sub(f"__", '_', filename)
        filename = re.sub(f"_$", '', filename)

        filename = add_file_extension(filename, extension=extension)

        return filename

    def written_as(self, filestring: str=dflt_frc['fs']['General'], datestring: str=dflt_frc['ds']['General'], extension: str='', defaults: str=''):
        # E.g. defaults='SWAN' uses all SWAN defaults
        if defaults:
            filestring = dflt_frc['fs'][defaults]
            datestring = dflt_frc['ds'][defaults]
            extension = dflt_frc['ext'][defaults]

        if hasattr(self, '_written_as'):
            filename = self._written_as
        else:
            filename = self.filename(filestring=filestring, datestring=datestring, extension=extension)

        return filename

    def written_to(self, folder: str=dflt_frc['fldr']['General']):
        if hasattr(self, '_written_to'):
            return self._written_to
        else:
            return folder

    def is_written(self):
        return hasattr(self, '_written_as')

    def time(self):
        return copy(pd.to_datetime(self.data.time.values))

    # def start_time(self):
    #     return self.time()[0]
    #
    # def end_time(self):
    #     return self.time()[-1]

    def u(self):
        return copy(self.data.u.values)

    def v(self):
        return copy(self.data.v.values)

    def nx(self):
        return (self.data.u.shape[1])

    def ny(self):
        return (self.data.u.shape[2])

    def nt(self):
        return (self.data.u.shape[0])

    def lon(self):
        """Returns a longitude vector of the grid."""
        if hasattr(self.data, 'lon'):
            lon = copy(self.data.lon.values)
        else:
            lon = np.array([])
        return lon

    def lat(self):
        """Returns a latitude vector of the grid."""
        if hasattr(self.data, 'lat'):
            lat = copy(self.data.lat.values)
        else:
            lat = np.array([])
        return lat

    def size(self) -> tuple:
        """Returns the size (nx, ny) of the grid."""
        return self.data.u.shape

    def _point_list(self, mask):
        """Provides a list on longitudes and latitudes with a given mask.

        Used to e.g. generate list of boundary points or land points.
        """
        meshlon, meshlat=np.meshgrid(self.lon(),self.lat())
        lonlat_flat = np.column_stack((meshlon.ravel(),meshlat.ravel()))
        mask_flat = mask.ravel()

        return lonlat_flat[mask_flat]


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

class ForcingWriter(ABC):

    @abstractmethod
    def __call__(self, forcing_out: Forcing) -> None:
        pass

    def create_filename(self, forcing_out: Forcing, forcing_in_filename: bool=True, grid_in_filename: bool=True, time_in_filename: bool=True) -> str:
        """Creates a filename based on the boolean swithes set in __init__ and the meta data in the objects"""

        forcing_fn = ''
        grid_fn = ''
        time_fn = ''

        if forcing_in_filename:
            forcing_fn = f"_{forcing_out.name()}"

        if grid_in_filename:
            grid_fn = f"_{forcing_out.grid.name()}"

        if time_in_filename:
            time_fn = f"_{str(forcing_out.time()[0])[0:10]}_{str(forcing_out.time()[-1])[0:10]}"

        filename = forcing_fn + grid_fn + time_fn

        return filename
