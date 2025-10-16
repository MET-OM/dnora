from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_time, add_datavar, add_magnitude


@add_datavar(name="sss", default_value=0.0)
@add_datavar(name="sst", default_value=0.0)
@add_time(grid_coord=True)
class Ocean(GriddedSkeleton):
    pass
