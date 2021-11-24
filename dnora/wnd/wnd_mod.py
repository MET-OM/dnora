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

from .read import ForcingReader
from .write import ForcingWriter

class Forcing:
    def __init__(self, grid, name='AnonymousForcing'):
        self.grid = copy(grid)
        self._name = copy(name)
        return

    def import_forcing(self, start_time: str, end_time: str, forcing_reader: ForcingReader, expansion_factor: float=1.2):
        """Imports forcing data from a certain source.

        Data are import between start_time and end_time from the source
        defined in the forcing_reader. Data are read around an area defined
        by the Grid object passed at initialization of this object.
        """

        self.start_time = copy(start_time)
        self.end_time = copy(end_time)

        msg.header(
            f"{type(forcing_reader).__name__}: Loading wind forcing...")
        self.data = forcing_reader(
            self.grid, start_time, end_time, expansion_factor)

        return

    def export_forcing(self, forcing_writer) -> None:
        """Exports the forcing data to a file.

        The forcing_writer defines the file format.
        """

        output_file, output_folder = forcing_writer(self)

        # This is set as info in case an input file needs to be generated
        self._written_as = output_file
        self._written_to = output_folder

        return


    def days(self):
        """Determins a Pandas data range of all the days in the time span."""

        days = day_list(start_time=self.start_time, end_time=self.end_time)
        return days

    def name(self) -> str:
        """Return the name of the grid (set at initialization)."""

        return copy(self._name)


    def filename(self, filestring: str=dflt_frc['fs']['General'], datestring: str=dflt_frc['ds']['General'], extension: str='', defaults: str=''):
        """Creates a filename for the object.

        The filename can be based on e.g. the name of the Grid or Boundary
        object itself, or the start and end times.

        This is typically called by a ForcingWriter object when using
        the .export_forcing() method.
        """

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
        """Provide the filename the object has been exported to.

        If it has not been exported, a filename is created based on the
        metadata of the object / filestring provided in the function call.

        This is typically called when an input file for the model run needs
        to be created.
        """

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
        """Provide the folder the object has been exported to.

        If it has not been exported, a folder is created based on the
        metadata of the object / filestring provided in the function call.

        This is typically called when an input file for the model run needs
        to be created.
        """

        if hasattr(self, '_written_to'):
            return self._written_to
        else:
            return folder

    def is_written(self):
        """True / False statement to check if the object has ever been
        exported with .export_forcing()."""

        return hasattr(self, '_written_as')

    def time(self):
        return copy(pd.to_datetime(self.data.time.values))

    def u(self):
        return copy(self.data.u.values)

    def v(self):
        return copy(self.data.v.values)

    def nx(self):
        return (self.data.u.shape[2])

    def ny(self):
        return (self.data.u.shape[1])

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
