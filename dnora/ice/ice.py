from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_datavar, add_time

import geo_parameters as gp


@add_datavar(name="thickness", default_value=0.0)
@add_datavar(name="concentration", default_value=0.0)
@add_time(grid_coord=True)
class Ice(GriddedSkeleton):
    meta_dict = {
        "thickness": gp.ocean.IceThickness,
        "concentration": gp.ocean.IceFraction,
    }
    pass
