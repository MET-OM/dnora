from geo_skeletons import PointSkeleton
from geo_skeletons.decorators import add_mask, add_time

from copy import deepcopy


@add_mask(name="buoy", coord_group="grid", default_value=1)
@add_time(grid_coord=False)
class WaveSeries(PointSkeleton):
    pass
