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
from ..skeletons.point_skeleton import PointSkeleton
from ..skeletons.coordinate_factory import add_time, add_frequency
from ..skeletons.mask_factory import add_mask
from ..skeletons.datavar_factory import add_datavar

from .process import SpectralProcessor

@add_mask(name='bad', coords='all', default_value=0)
@add_datavar(name='spec', coords='all', default_value=0.)
@add_datavar(name='mdir', coords='all', default_value=0.)
@add_datavar(name='spr', coords='all', default_value=0.)
@add_frequency(grid_coord=False)
@add_time(grid_coord=True)
class Spectra(PointSkeleton):
    def __init__(self, grid: Grid, name: str="AnonymousSpectra"):
        self.grid = copy(grid)
        self._name = copy(name)
        self._convention = None
        self._history = []

    def import_spectra(self, start_time: str, end_time: str, spectral_reader: SpectralReader,  point_picker: PointPicker = TrivialPicker()) -> None:
        """Imports omnidirectional spectra from a certain source.

        Spectra are import between start_time and end_time from the source
        defined in the spectral_reader. Which spectra to choose spatially
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
        time, freq, spec, mdir, spr, lon, lat, x, y, attributes = spectral_reader(self.start_time, self.end_time, inds)

        self._init_structure(x, y, lon, lat, time=time, freq=freq)

        self.ds_manager.set(spec, 'spec', coord_type='all')
        self.ds_manager.set(mdir, 'mdir', coord_type='all')
        self.ds_manager.set(spr, 'spr', coord_type='all')

        self.ds_manager.set_attrs(attributes)

        # self.data = self.compile_to_xr(time, freq, spec, mdir, spr, lon, lat, source)
        # self.mask = [True]*len(self.x())

        # E.g. are the spectra oceanic convention etc.
        self._convention = spectral_reader.convention()

    def process_spectra(self, spectral_processors: List[SpectralProcessor]=None):
        """Process all the individual spectra of the spectra object.

        E.g. change convention form WW3 to Oceanic, interpolate spectra to
        new frequency grid, or multiply everything with a constant.
        """

        if spectral_processors is None:
            msg.info("No SpectralProcessor provided. Doing Nothing.")
            return

        if not isinstance(spectral_processors, list):
            spectral_processors = [spectral_processors]

        convention_warning = False

        for processor in spectral_processors:

            msg.process(f"Processing spectra with {type(processor).__name__}")
            self._history.append(copy(processor))
            old_convention = processor._convention_in()
            if old_convention is not None:
                if old_convention != self.convention():
                    msg.warning(f"Spectral convention ({self.convention()}) doesn't match that expected by the processor ({old_convention})!")
                    convention_warning=True


            new_spec, new_dirs, new_freq = processor(self.spec(), self.mdir(), self.freq())
            self._init_structure(x=self.x(strict=True), y=self.y(strict=True),
                            lon=self.lon(strict=True), lat=self.lat(strict=True),
                            time=self.time(), freq=new_freq)
            self.ds_manager.set(new_spec, 'spec', coord_type='all')

            # self.data.spec.values = new_spec
            # self.data = self.data.assign_coords(dirs=new_dirs)
            # self.data = self.data.assign_coords(freq=new_freq)

            # Set new convention if the processor changed it
            new_convention = processor._convention_out()
            if new_convention is not None:
                self._convention = new_convention
                if convention_warning:
                    msg.warning(f"Convention variable set to {new_convention}, but this might be wrong...")
                else:
                    msg.info(f"Changing convention from {old_convention} >>> {new_convention}")

            print(processor)
            msg.blank()
        return


    def convention(self):
        """Returns the convention (OCEAN/MET/MATH) of the spectra"""
        if not hasattr(self, '_convention'):
            return None
        return copy(self._convention)


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
