from geo_skeletons import PointSkeleton
from geo_skeletons.decorators import add_mask, add_time, add_datavar
import geo_parameters as gp

from copy import deepcopy

@add_datavar(gp.wave.Tm02)
@add_datavar(gp.wave.Tm_10)
@add_datavar(gp.wave.Tm01)
@add_datavar(gp.wave.Dirm)
@add_datavar(gp.wave.Dirp)
@add_datavar(gp.wave.Tp)
@add_datavar(gp.wave.Hs)
@add_mask(name="buoy", coord_group="grid", default_value=1)
@add_time(grid_coord=False)
class WaveSeries(PointSkeleton):
    pass
