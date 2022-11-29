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
from ..grd.grd_mod import Grid, UnstrGrid
from .conventions import SpectralConvention
from .process import boundary_processor_for_convention_change
# Import abstract classes and needed instances of them
from .process import BoundaryProcessor, Multiply
from .pick import PointPicker, TrivialPicker
from .read import BoundaryReader
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .write import BoundaryWriter # Abstract class

# Import default values and aux_funcsiliry functions
from .. import msg
from .. import aux_funcs
from ..skeletons.point_skeleton import PointSkeleton
from ..skeletons.coordinate_factory import add_time, add_frequency, add_direction
from ..skeletons.mask_factory import add_mask
from ..skeletons.datavar_factory import add_datavar

@add_mask(name='bad', coords='all', default_value=0)
@add_datavar(name='spec', coords='all', default_value=0.)
@add_direction(grid_coord=False)
@add_frequency(grid_coord=False)
@add_time(grid_coord=True)
class Boundary(PointSkeleton):
    def __init__(self, grid: Grid, name: str="AnonymousBoundary"):
        self._grid = grid
        self._name = name
        self._convention = None
        self._history = []

    def import_boundary(self, start_time: str, end_time: str,
                        boundary_reader: BoundaryReader,
                        point_picker: PointPicker,
                        expansion_factor: float=1.5):
        """Imports boundary spectra from a certain source.

        Spectra are import between start_time and end_time from the source
        defined in the boundary_reader. Which spectra to choose spatically
        are determined by the point_picker.
        """

        self._history.append(copy(boundary_reader))

        msg.header(boundary_reader, "Reading coordinates of spectra...")
        lon_all, lat_all, x_all, y_all = boundary_reader.get_coordinates(self.grid(), start_time)
        all_points = UnstrGrid(lon=lon_all, lat=lat_all, x=x_all, y=y_all)

        msg.header(point_picker, "Choosing boundary spectra...")
        inds = point_picker(self.grid(), all_points, expansion_factor)

        if len(inds) < 1:
            msg.warning("PointPicker didn't find any points. Aborting import of boundary.")
            return

        # Main reading happens here
        msg.header(boundary_reader, "Loading boundary spectra...")

        time, freq, dirs, spec, lon, lat, x, y, metadata = boundary_reader(self.grid(), start_time, end_time, inds)

        self._init_structure(x, y, lon, lat, time=time, freq=freq, dirs=dirs)

        self.ds_manager.set(spec, 'spec', coord_type='all')
        self.set_metadata(metadata)
        # E.g. are the spectra oceanic convention etc.
        self._convention = boundary_reader.convention()

        self.set_metadata({'spectral_convention': self.convention().value}, append=True)

        if boundary_reader.post_processing() is not None:
            self.process_boundary(boundary_reader.post_processing())
        return

    def process_boundary(self, boundary_processors: List[BoundaryProcessor]=None):
        """Process all the individual spectra of the boundary object.

        E.g. change convention form WW3 to Oceanic, interpolate spectra to
        new directional grid, or multiply everything with a constant.
        """

        if boundary_processors is None:
            msg.info("No BoundaryProcessor provided. Doing Nothing.")
            return

        if not isinstance(boundary_processors, list):
            boundary_processors = [boundary_processors]

        convention_warning = False

        for processor in boundary_processors:

            msg.process(f"Processing spectra with {type(processor).__name__}")
            print(processor)
            self._history.append(copy(processor))
            old_convention = processor._convention_in()
            if old_convention is not None:
                if old_convention != self.convention():
                    msg.warning(f"Boundary convention ({self.convention()}) doesn't match that expected by the processor ({old_convention})!")
                    convention_warning=True


            new_spec, new_dirs, new_freq, new_inds = processor(self.spec(), self.dirs(), self.freq(), self.inds())
            new_inds = list(new_inds)

            if new_inds:
                metadata = self.metadata()

                self._init_structure(x=self.x(strict=True, inds=new_inds), y=self.y(strict=True, inds=new_inds),
                                lon=self.lon(strict=True, inds=new_inds), lat=self.lat(strict=True, inds=new_inds),
                                time=self.time(), freq=new_freq, dirs=new_dirs)
                self.ds_manager.set(new_spec, 'spec', coord_type='all')
                self.set_metadata(metadata) # Global attributes

                # Set new convention if the processor changed it
                new_convention = processor._convention_out()
                if new_convention is not None:
                    self._set_convention(new_convention, process=False)
                    if convention_warning:
                        msg.warning(f"Convention variable set to {new_convention}, but this might be wrong...")
                    else:
                        msg.info(f"Changing convention from {old_convention} >>> {new_convention}")
            else:
                self.ds_manager.set_new_ds(None)
                msg.warning(f"No boundary spectra left after processing. Removing all data.")

            msg.blank()
        return

    def _set_convention(self, convention: SpectralConvention, process: bool=True) -> None:
        """Sets a new spectral directional convention. To not touch spectra, use process=False."""
        if isinstance(convention, str):
            convention = SpectralConvention[convention.upper()]

        boundary_processor = boundary_processor_for_convention_change(
                            current_convention = self.convention(),
                            wanted_convention = convention)

        if boundary_processor is None:
            msg.info(f"Convention ({self.convention()}) already equals wanted convention ({convention}).")
            return

        if process:
            self.process_boundary(boundary_processor)
        else:
            self._convention = convention
            self.set_metadata({'spectral_convention': self.convention().value}, append=True)
            print(f'Spectral convention is now: {self.convention()}')

    def convention(self):
        """Returns the convention (WW3/OCEAN/MET/MATH/MATHVEC) of the spectra"""
        if not hasattr(self, '_convention'):
            return None
        return copy(self._convention)

    def grid(self) -> Grid:
        if hasattr(self, '_grid'):
            return self._grid
        return None

    def __str__(self) -> str:
        """Prints status of boundary."""

        msg.header(self, f"Status of boundary {self.name}")
        if self.x() is not None:
            msg.plain(f"Contains data ({len(self.x())} points) for {self.start_time()} - {self.end_time()}")
            msg.plain(f"Data covers: lon: {min(self.lon())} - {max(self.lon())}, lat: {min(self.lat())} - {max(self.lat())}")
        if len(self._history) > 0:
            msg.blank()
            msg.plain("Object has the following history:")
            for obj in self._history:
                msg.process(f"{obj.__class__.__bases__[0].__name__}: {type(obj).__name__}")

        msg.print_line()

        return ''
