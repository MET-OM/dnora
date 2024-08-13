from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_datavar, add_time, add_magnitude
import geo_parameters as gp


@add_magnitude(
    gp.ocean.Current("mag"), x="u", y="v", direction=gp.ocean.CurrentDir("dir")
)
@add_datavar(name=gp.ocean.YCurrent("v"), default_value=0.0)
@add_datavar(name=gp.ocean.XCurrent("u"), default_value=0.0)
@add_time(grid_coord=True)
class Current(GriddedSkeleton):
    pass
