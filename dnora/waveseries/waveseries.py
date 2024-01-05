from geo_skeletons import PointSkeleton
from geo_skeletons.decorators import add_mask, add_time

from copy import deepcopy


@add_mask(name="buoy", coords="grid", default_value=1)
@add_time(grid_coord=False)
class WaveSeries(PointSkeleton):
    def __init__(self, x=None, y=None, lon=None, lat=None, **kwargs):
        super().__init__(x, y, lon, lat, **kwargs)
        self._coord_manager = deepcopy(
            WaveSeries._coord_manager
        )  # We are dynamically adding data variables to the instance
