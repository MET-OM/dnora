from __future__ import annotations # For TYPE_CHECKING
from abc import ABC, abstractmethod
import xarray as xr
import numpy as np
from copy import copy
import pandas as pd
from typing import List
import sys
import matplotlib.pyplot as plt
import re

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .write import BoundaryWriter # Abstract class

# dnora imports
from .. import msg
from ..aux import distance_2points, day_list, create_filename_obj, create_filename_time, create_filename_lonlat, add_file_extension

from .process import processor_for_convention_change, Multiply, SpectralProcessor
from .pick import PointPicker, TrivialPicker
from .read import BoundaryReader
from .write import BoundaryWriter
from ..grd.grd_mod import Grid # Grid object

from ..defaults import dflt_bnd


class Boundary:
    def __init__(self, grid: Grid, name: str = "AnonymousBoundary"):
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

        msg.header(f"{type(boundary_reader).__name__}: Reading coordinats of spectra...")
        lon_all, lat_all = boundary_reader.get_coordinates(self.start_time)

        msg.header(f"Choosing spectra with {type(point_picker).__name__}")
        inds = point_picker(self.grid, lon_all, lat_all)

        msg.header(f"{type(boundary_reader).__name__}: Loading boundary spectra...")
        time, freq, dirs, spec, lon, lat, source = boundary_reader(self.start_time, end_time, inds)

        self.data = self.compile_to_xr(time, freq, dirs, spec, lon, lat, source)
        self.mask = [True]*len(self.x())

        # E.g. are the spectra oceanic convention etc.
        self._convention = boundary_reader.get_convention()

        return


    def process_spectra(self, spectral_processors: List[SpectralProcessor] = [Multiply(calib_spec = 1)]):
        """Process all the individual spectra of the boundary object.

        E.g. change convention form WW3 to Oceanic, interpolate spectra to
        new directional grid, or multiply everything with a constant.
        """

        if not isinstance(spectral_processors, list):
            spectral_processors = [spectral_processors]

        for n in range (len(spectral_processors)):
            spectral_processor = spectral_processors[n]

            msg.process(f"Processing spectra with {type(spectral_processor).__name__}")
            print(spectral_processor)

            new_spec, new_dirs, new_freq = spectral_processor(self.spec(), self.dirs(), self.freq())

            self.data.spec.values = new_spec
            self.data = self.data.assign_coords(dirs=new_dirs)
            self.data = self.data.assign_coords(freq=new_freq)

            # Set new convention if the processor changed it
            new_convention = spectral_processor.get_convention()
            if new_convention is not None:
                msg.info(f"Setting new convention to {new_convention}")
                self._convention = new_convention
        return

    def export_boundary(self, boundary_writer: BoundaryWriter) -> None:
        """Exports the boundary spectra to a file.

        The bounday_writer defines the file format.
        """

        output_file, output_folder = boundary_writer(self)

        # This is set as info in case an input file needs to be generated
        self._written_as = output_file
        self._written_to = output_folder
        return

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

        spectral_processor = processor_for_convention_change(current_convention = self.convention(), wanted_convention = wanted_convention)
        if spectral_processor is not None:
            self.process_spectra(spectral_processor)
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

    def filename(self, filestring: str=dflt_bnd['fs']['General'], datestring: str=dflt_bnd['ds']['General'], n: int=None, extension: str='', defaults: str=''):
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
            extension = dflt_bnd['ext'][defaults]

        # Substitute placeholders for objects ($Grid etc.)
        filename = create_filename_obj(filestring=filestring, objects=[self, self.grid])
        # Substitute placeholders for times ($T0 etc.)
        filename = create_filename_time(filestring=filename, times=[self.start_time, self.end_time], datestring=datestring)

        # Substitute $Lon and $Lat if a single output point is specified
        if n is not None:
            filename = create_filename_lonlat(filename, lon=self.lon()[n], lat=self.lat()[n])
        else:
            # Trying to remove possible mentions to longitude and latitude
            filename = re.sub(f"E\$Lon", '', filename)
            filename = re.sub(f"N\$Lat", '', filename)
            filename = re.sub(f"\$Lon", '', filename)
            filename = re.sub(f"\$Lat", '', filename)

        # Possible clean up
        filename = re.sub(f"__", '_', filename)
        filename = re.sub(f"_$", '', filename)

        filename = add_file_extension(filename, extension=extension)

        return filename

    def written_as(self, filestring: str=dflt_bnd['fs']['General'], datestring: str=dflt_bnd['ds']['General'], extension: str='', defaults: str='') -> str:
        """Provide the filename the object has been exported to.

        If it has not been exported, a filename is created based on the
        metadata of the object / filestring provided in the function call.

        This is typically called when an input file for the model run needs
        to be created.
        """

        # E.g. defaults='SWAN' uses all SWAN defaults
        if defaults:
            filestring = dflt_bnd['fs'][defaults]
            datestring = dflt_bnd['ds'][defaults]
            extension = dflt_bnd['ext'][defaults]

        if hasattr(self, '_written_as'):
            filename = self._written_as
        else:
            filename =  self.filename(filestring=filestring, datestring=datestring, extension=extenstion)

        return filename

    def written_to(self, folder: str=dflt_bnd['fldr']['General']) -> str:
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

    def is_written(self) -> bool:
        """True / False statement to check if the object has ever been
        exported with .export_boundary()."""

        return hasattr(self, '_written_as')

    def time(self):
        return copy(pd.to_datetime(self.data.time.values))

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
