
from skeletons.point_skeleton import PointSkeleton
from skeletons.mask_factory import add_mask
from skeletons.coordinate_factory import add_time
from copy import deepcopy
import numpy as np
import pandas as pd
#@add_mask(name='bad', coords='all', default_value=0)
#@add_datavar(name='spec', coords='all', default_value=0.)
@add_mask(name='buoy', coords='grid', default_value=1)
@add_time(grid_coord=False)
class WaveSeries(PointSkeleton):
    def __init__(self, x=None, y=None, lon=None, lat=None, time= pd.date_range('1990-01-01 00:00', '1990-01-01 00:00', freq='H'), name: str='LonelyWaveSeries', **kwargs):
        if np.all([a is None for a in [x,y,lon,lat]]):
            x, y = 0, 0
        super().__init__(x=x, y=y, lon=lon, lat=lat, name=name, time=time, **kwargs)
        self._coord_manager = deepcopy(WaveSeries._coord_manager) # We are dynamically adding data variables to the instance



