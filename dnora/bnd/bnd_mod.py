from __future__ import annotations # For TYPE_CHECKING

from abc import ABC, abstractmethod
import xarray as xr
import numpy as np
from copy import copy
import pandas as pd
from typing import List
import sys
import re
import glob, os
# Import objects
from ..grd.grd_mod import Grid
from .conventions import SpectralConvention, convention_from_string
from .process import boundary_processor_for_convention_change
# Import abstract classes and needed instances of them
from .process import BoundaryProcessor, Multiply
from .pick import PointPicker, TrivialPicker
from .read import BoundaryReader, DnoraNc
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .write import BoundaryWriter # Abstract class

# Import default values and aux_funcsiliry functions
from .. import msg
from .. import aux_funcs
from ..cacher import Cacher
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
        self.grid = copy(grid)
        self._name = copy(name)
        self._convention = None
        self._history = []

    def import_boundary(self, start_time: str, end_time: str,
                        boundary_reader: BoundaryReader,
                        point_picker: PointPicker = TrivialPicker(),
                        write_cache: bool=False,
                        read_cache: bool=False,
                        cache_name: str=None):
        """Imports boundary spectra from a certain source.

        Spectra are import between start_time and end_time from the source
        defined in the boundary_reader. Which spectra to choose spatically
        are determined by the point_picker.
        """
        ### Some boundary readers may need to define an area and grid when
        ### reading the database. This is set to make that possible.
        ### The regular grid doesn't necessarily match the boundary points
        ### exactly.
        if not np.all(np.logical_not(self.grid.boundary_mask())): # Boundary mask empty?
            boundary_point_grid = Grid(lon=self.grid.edges('lon'),
                                        lat=self.grid.edges('lat'),
                                        name='boundary_points')
            if len(boundary_point_grid.lon()) > 1 and len(boundary_point_grid.lat()) > 1:
                boundary_point_grid.set_spacing(nx=self.grid.boundary_nx(),
                                                ny=self.grid.boundary_ny())
            boundary_reader.set_restricted_area(boundary_point_grid)

        # Prepare for working with cahced data if we have to
        if write_cache or read_cache:
            cacher = Cacher(self, boundary_reader.name(), cache_name)

        # Read whatever we have in the chached data to start with
        # Setting the reader to read standard DNORA netcdf-files
        if read_cache and not cacher.empty():
            msg.info('Reading boundary data from cache!!!')
            original_boundary_reader = copy(boundary_reader)
            boundary_reader = DnoraNc(files=glob.glob(f'{cacher.filepath(extension=False)}*'), convention = boundary_reader.convention())

        self._history.append(copy(boundary_reader))

        msg.header(boundary_reader, "Reading coordinates of spectra...")
        lon_all, lat_all, x_all, y_all = boundary_reader.get_coordinates(start_time)

        msg.header(point_picker, "Choosing spectra...")
        inds = point_picker(self.grid, lon_all, lat_all)

        if len(inds) < 1:
            msg.warning("PointPicker didn't find any points. Aborting import of boundary.")
            return

        # Main reading happens here
        msg.header(boundary_reader, "Loading boundary spectra...")
        time, freq, dirs, spec, lon, lat, x, y, attributes = boundary_reader(start_time, end_time, inds)
        self._init_structure(x, y, lon, lat, time=time, freq=freq, dirs=dirs)
        self.ds_manager.set(spec, 'spec', coord_type='all')
        self.ds_manager.set_attrs(attributes)

        # Patch data if read from cache and all data not found
        if read_cache and not cacher.empty():
            patch_start, patch_end = cacher.determine_patch_periods(start_time, end_time)
            if patch_start:
                msg.info('Not all data found in cache. Patching from original source...')

                for t0, t1 in zip(patch_start, patch_end):
                    boundary_temp = Boundary(self.grid)
                    boundary_temp.import_boundary(start_time=t0, end_time=t1,
                                    boundary_reader=original_boundary_reader,
                                    point_picker=point_picker)
                    self._absorb_object(boundary_temp, 'time')

        # Dump monthly netcdf-files that will now be in standard DNORA format
        if write_cache:
            msg.info('Caching data:')
            cacher.write_cache()

        # E.g. are the spectra oceanic convention etc.
        self._convention = boundary_reader.convention()
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
            self._history.append(copy(processor))
            old_convention = processor._convention_in()
            if old_convention is not None:
                if old_convention != self.convention():
                    msg.warning(f"Boundary convention ({self.convention()}) doesn't match that expected by the processor ({old_convention})!")
                    convention_warning=True


            new_spec, new_dirs, new_freq = processor(self.spec(), self.dirs(), self.freq())
            self._init_structure(x=self.x(strict=True), y=self.y(strict=True),
                            lon=self.lon(strict=True), lat=self.lat(strict=True),
                            time=self.time(), freq=new_freq, dirs=new_dirs)
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

    def _set_convention(self, convention: SpectralConvention) -> None:
        boundary_processor = boundary_processor_for_convention_change(
                            current_convention = self.convention(),
                            wanted_convention = convention)

        if boundary_processor is None:
            msg.info(f"Convention ({self.convention()}) already equals wanted convention ({convention}).")
        else:
            self.process_boundary(boundary_processor)

    def convention(self):
        """Returns the convention (WW3/OCEAN/MET/MATH/MATHVEC) of the spectra"""
        if not hasattr(self, '_convention'):
            return None
        return copy(self._convention)


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
