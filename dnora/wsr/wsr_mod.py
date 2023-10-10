from geo_skeletons import PointSkeleton
from geo_skeletons.decorators import add_mask, add_time

from copy import deepcopy
import numpy as np
import pandas as pd


# @add_mask(name='bad', coords='all', default_value=0)
# @add_datavar(name='spec', coords='all', default_value=0.)
@add_mask(name="buoy", coords="grid", default_value=1)
@add_time(grid_coord=False)
class WaveSeries(PointSkeleton):
    def __init__(self, x, y, lon, lat, **kwargs):
        # if np.all([a is None for a in [x, y, lon, lat]]):
        #     x, y = 0, 0
        super().__init__(x, y, lon, lat, **kwargs)
        self._coord_manager = deepcopy(
            WaveSeries._coord_manager
        )  # We are dynamically adding data variables to the instance
