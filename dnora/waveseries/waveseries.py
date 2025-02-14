from geo_skeletons import PointSkeleton
from geo_skeletons.decorators import add_mask, add_time, add_datavar, add_magnitude
import geo_parameters as gp

from copy import deepcopy
import numpy as np


@add_magnitude(
    gp.ocean.Current, x="x_current", y="y_current", direction=gp.ocean.CurrentDir
)
@add_datavar(gp.ocean.YCurrent, default_value=0.1)
@add_datavar(gp.ocean.XCurrent, default_value=0.1)
@add_magnitude(gp.wind.Wind, x="x_wind", y="y_wind", direction=gp.wind.WindDir)
@add_datavar(gp.wind.YWind, default_value=0.1)
@add_datavar(gp.wind.XWind, default_value=0.1)
@add_datavar(gp.ocean.WaterDepth)
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
