from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_time, add_datavar, add_magnitude
import numpy as np
import geo_parameters as gp


@add_datavar(name=gp.ocean.SeaSurfaceSalinity("sss"), default_value=np.nan)
@add_datavar(name=gp.ocean.SeaSurfaceTemperature("sst"), default_value=np.nan)
@add_time(grid_coord=True)
class Ocean(GriddedSkeleton):
    pass
