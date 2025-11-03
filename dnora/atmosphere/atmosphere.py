from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_time, add_datavar, add_magnitude
import geo_parameters as gp
import numpy as np


@add_datavar(name=gp.atm.MeanSeaLevelPressure("mslp"), default_value=np.nan)
@add_datavar(name=gp.atm.RelativeHumidity("r"), default_value=np.nan)
@add_datavar(name=gp.atm.AirTemperature("t2m"), default_value=np.nan)
@add_time(grid_coord=True)
class Atmosphere(GriddedSkeleton):
    pass
