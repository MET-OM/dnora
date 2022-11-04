import xarray as xr
from typing import List
import numpy as np
import matplotlib.pyplot as plt
from typing import List
from copy import copy
import pandas as pd
from .. import msg
from .read import WaveSeriesReader
from ..grd import Grid
from .wave_parameters import WaveParameter
from ..skeletons.point_skeleton import PointSkeleton
from ..skeletons.datavar_factory import add_datavar
from ..skeletons.coordinate_factory import add_time, add_frequency, add_direction
from ..bnd.pick import PointPicker, TrivialPicker
#from ..skeletons.mask_factory import add_mask
#from ..skeletons.datavar_factory import add_datavar
from copy import deepcopy
#@add_mask(name='bad', coords='all', default_value=0)
#@add_datavar(name='spec', coords='all', default_value=0.)
@add_time(grid_coord=True)
class WaveSeries(PointSkeleton):
    #def __call__(self, parameter: str) -> np.ndarray:
    #    return self.data[parameter].values
    def __init__(self, grid: Grid, name: str="LonelyWaveSeries"):
        self._grid = grid
        self._name = name
        self._convention = None
        self._history = []
        self._coord_manager = deepcopy(WaveSeries._coord_manager) # We are dynamically adding data variables to the instance

    def import_waveseries(self, start_time: str, end_time: str,
                        waveseries_reader: WaveSeriesReader,
                        point_picker: PointPicker,
                        expansion_factor: float=1.5) -> None:
        self.start_time = copy(start_time)
        self.end_time = copy(end_time)
        self._history.append(copy(waveseries_reader))

        msg.header(waveseries_reader, "Reading coordinates of WaveSeries...")
        lon_all, lat_all = waveseries_reader.get_coordinates(self.grid(), self.start_time)

        msg.header(point_picker, "Choosing waves series points...")
        inds = point_picker(self.grid(), lon_all, lat_all, expansion_factor)

        msg.header(waveseries_reader, "Loading wave series data...")
        time, data_dict, lon, lat, x, y, attributes = waveseries_reader(self.grid(), start_time, end_time, inds)

        self._init_structure(x, y, lon, lat, time=time)

        for wp, data in data_dict.items():
            self.ds_manager.set(data, wp.name(), coord_type='all')
            self.ds_manager.set_attrs({'name': wp.name(), 'unit': wp.unit(), 'standard_name': wp.standard_name()}, wp.name())
            self = add_datavar(wp.name(), aftermath=True)(self) # Creates .hs() etc. methods
        self.ds_manager.set_attrs(attributes) # Global attributes

    def grid(self) -> Grid:
        if hasattr(self, '_grid'):
            return self._grid
        return None
