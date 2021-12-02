from __future__ import annotations # For TYPE_CHECKING

from abc import ABC, abstractmethod
import xarray as xr
import numpy as np
from copy import copy
import pandas as pd
from typing import List
import sys
import re

# Import objects
from ..grd.grd_mod import Grid

# Import abstract classes and needed instances of them
from .process import BoundaryProcessor, Multiply
from .pick import PointPicker, TrivialPicker
from .read import BoundaryReader
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .write import BoundaryWriter # Abstract class

# Import default values and auxiliry functions
from .. import msg
from ..aux import day_list, create_filename_obj, create_filename_time, create_filename_lonlat, clean_filename, check_if_folder
from .process import processor_for_convention_change
from ..defaults import dflt_bnd, list_of_placeholders


class Boundary:
    def __init__(self, grid: Grid, name: str="AnonymousBoundary"):
        self.grid = copy(grid)
        self._name = copy(name)
        self._convention = None
        return

    def import_boundary(self, start_time: str, end_time: str, boundary_reader: BoundaryReader,  point_picker: PointPicker = TrivialPicker()):
        """Imports boundary spectra from a certain source.

        Spectra are import between start_time and end_time from the source
        defined in the boundary_reader. Which spectra to choose spatically
        are determined by the point_picker.
        """

        self.start_time = copy(start_time)
        self.end_time = copy(end_time)

        msg.header(boundary_reader, "Reading coordinates of spectra...")
        lon_all, lat_all = boundary_reader.get_coordinates(self.start_time)

        msg.header(point_picker, "Choosing spectra...")
        inds = point_picker(self.grid, lon_all, lat_all)

        msg.header(boundary_reader, "Loading boundary spectra...")
        time, freq, dirs, spec, lon, lat, source = boundary_reader(self.start_time, end_time, inds)

        self.data = self.compile_to_xr(time, freq, dirs, spec, lon, lat, source)
        self.mask = [True]*len(self.x())

        # E.g. are the spectra oceanic convention etc.
        self._convention = boundary_reader.convention()

        return


    def process_boundary(self, boundary_processors: List[BoundaryProcessor]=[Multiply(calib_spec = 1)]):
        """Process all the individual spectra of the boundary object.

        E.g. change convention form WW3 to Oceanic, interpolate spectra to
        new directional grid, or multiply everything with a constant.
        """

        if not isinstance(boundary_processors, list):
            boundary_processors = [boundary_processors]

        for n in range (len(boundary_processors)):
            boundary_processor = boundary_processors[n]

            msg.process(f"Processing spectra with {type(boundary_processor).__name__}")
            print(boundary_processor)

            new_spec, new_dirs, new_freq = boundary_processor(self.spec(), self.dirs(), self.freq())

            self.data.spec.values = new_spec
            self.data = self.data.assign_coords(dirs=new_dirs)
            self.data = self.data.assign_coords(freq=new_freq)

            # Set new convention if the processor changed it
            new_convention = boundary_processor.convention()
            if new_convention is not None:
                msg.info(f"Setting new convention to {new_convention}")
                self._convention = new_convention
        return

    def export_boundary(self, boundary_writer: BoundaryWriter, out_format: str=None, filestring: str=None, datestring: str=None, folder: str=None) -> None:
        """Exports the boundary spectra to a file.

        The bounday_writer defines the file format.
        """

        # For setting the file name
        if out_format is None:
            out_format = boundary_writer._preferred_format()

        if filestring is None:
            filestring=dflt_bnd['fs'][out_format]

        if datestring is None:
            datestring=dflt_bnd['ds'][out_format]

        filename = self.filename(filestring=filestring, datestring=datestring)

        if folder is not None:
            folder = self.filename(filestring=folder)
        else:
            folder = self.filename(filestring=dflt_bnd['fldr'][out_format])

        existed = check_if_folder(folder=folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {folder}")

        msg.header(boundary_writer, f"Writing boundary spectra from {self.name()}")

        # Make sure convention is right for the reader
        wanted_convention = boundary_writer.convention()
        self.change_convention(wanted_convention=wanted_convention)

        output_files, output_folder = boundary_writer(self, filename=filename, folder=folder)

        return output_files, output_folder

    def filename(self, filestring: str=dflt_bnd['fs']['General'], datestring: str=dflt_bnd['ds']['General'], defaults: str=''):
        """Creates a filename for the object.

        The filename can be based on e.g. the name of the Grid or Boundary
        object itself, or the start and end times.

        This is typically called by a BoundaryWriter object when using
        the .export_boundary() method.
        """

        # E.g. defaults='SWAN' uses all SWAN defaults
        if defaults:
            filestring = dflt_bnd['fs'][defaults]
            datestring = dflt_bnd['ds'][defaults]

        # Substitute placeholders for objects ($Grid etc.)
        filename = create_filename_obj(filestring=filestring, objects=[self, self.grid])
        # Substitute placeholders for times ($T0 etc.)
        filename = create_filename_time(filestring=filename, times=[self.start_time, self.end_time], datestring=datestring)
        filename = clean_filename(filename, list_of_placeholders)

        # Substitute $Lon and $Lat if a single output point is specified
        # if n is not None:
        #     filename = create_filename_lonlat(filename, lon=self.lon()[n], lat=self.lat()[n])
        # else:
        #     # Trying to remove possible mentions to longitude and latitude
        #     filename = re.sub(f"E\$Lon", '', filename)
        #     filename = re.sub(f"N\$Lat", '', filename)
        #     filename = re.sub(f"\$Lon", '', filename)
        #     filename = re.sub(f"\$Lat", '', filename)
        #
        # # Possible clean up
        # filename = re.sub(f"__", '_', filename)
        # filename = re.sub(f"_$", '', filename)


        return filename

    def folder(self, folderstring: str=dflt_bnd['fldr']['General'], datestring: str=dflt_bnd['ds']['General']) -> str:
        # Substitute placeholders for $Grid
        folder = create_filename_obj(filestring=folderstring, objects=[self.grid, self])
        folder = create_filename_time(filestring=folder, times=[self.start_time, self.end_time], datestring=datestring)
        folder = clean_filename(folder, list_of_placeholders)

        return folder

    def change_convention(self, wanted_convention: str='') -> None:
        """Changes the convention of the spectra.

        The conventions to choose from are predetermined:

        'Ocean':    Oceanic convention
                    Directional vector monotonically increasing.
                    Direction to. North = 0, East = 90.

        'Met':      Meteorological convention
                    Directional vector monotonically increasing.
                    Direction from. North = 0, East = 90.

        'Math':     Mathematical convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        'WW3':      WAVEWATCH III output convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 0, East = 90.
        """

        boundary_processor = processor_for_convention_change(current_convention = self.convention(), wanted_convention = wanted_convention)
        if boundary_processor is not None:
            self.process_boundary(boundary_processor)
        return

    def compile_to_xr(self, time, freq, dirs, spec, lon, lat, source):
        """Data from .import_boundary() is stored as an xarray Dataset."""

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

    def time(self):
        return copy(pd.to_datetime(self.data.time.values))

    # def dt(self) -> float:
    #     """ Returns time step of boundary spectra in hours."""
    #     return self.time().to_series().diff().dt.total_seconds().values[-1]/3600

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

    def name(self) -> str:
        """Return the name of the grid (set at initialization)."""
        return copy(self._name)

    def convention(self) -> str:
        """Returns the convention (WW3/Ocean/Met/Math) of the spectra"""
        return copy(self._convention)

    def times_in_day(self, day):
        """Determines time stamps of one given day."""
        t0 = day.strftime('%Y-%m-%d') + "T00:00:00"
        t1 = day.strftime('%Y-%m-%d') + "T23:59:59"

        times = self.slice_data(start_time=t0, end_time=t1, x=[0]).time.values
        return times
