from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_datavar, add_time

import geo_parameters as gp


@add_datavar(gp.ocean.IceFraction("sic"), default_value=0.0)
@add_datavar(gp.ocean.IceThickness("sit"), default_value=0.0)
@add_time(grid_coord=True)
class Ice(GriddedSkeleton):
    pass
