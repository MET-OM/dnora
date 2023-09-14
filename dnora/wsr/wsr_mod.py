
from skeletons.point_skeleton import PointSkeleton
from skeletons.mask_factory import add_mask
from skeletons.coordinate_factory import add_time
from copy import deepcopy

#@add_mask(name='bad', coords='all', default_value=0)
#@add_datavar(name='spec', coords='all', default_value=0.)
@add_mask(name='buoy', coords='grid', default_value=1)
@add_time(grid_coord=False)
class WaveSeries(PointSkeleton):
    def __init__(self, x, y, lon, lat, name, **kwargs):
        super().__init__(x, y, lon, lat, name, **kwargs)
        self._coord_manager = deepcopy(WaveSeries._coord_manager) # We are dynamically adding data variables to the instance



