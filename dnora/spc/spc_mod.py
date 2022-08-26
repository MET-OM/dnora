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

@add_mask(name='bad', coords='all', default_value=0)
@add_datavar(name='spec', coords='all', default_value=0.)
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
        defined in the boundary_reader. Which spectra to choose spatially
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
        time, freq, spec, mdir, spr, lon, lat, x, y, source = spectral_reader(self.start_time, self.end_time, inds)

        self._init_structure(x, y, lon, lat, time=time, freq=freq)

        self.ds_manager.set(spec, 'spec', coord_type='all')
        self.ds_manager.set(mdir, 'mdir', coord_type='all')
        self.ds_manager.set(spr, 'spr', coord_type='all')

        # self.data = self.compile_to_xr(time, freq, spec, mdir, spr, lon, lat, source)
        # self.mask = [True]*len(self.x())

        # E.g. are the spectra oceanic convention etc.
        self._convention = spectral_reader.convention()
