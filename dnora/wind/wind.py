from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_time, add_datavar, add_magnitude

import geo_parameters as gp


@add_magnitude(gp.wind.Wind("mag"), x="u", y="v", direction=gp.wind.WindDir("dir"))
@add_datavar(name=gp.wind.YWind("v"), default_value=0.0)
@add_datavar(name=gp.wind.XWind("u"), default_value=0.0)
@add_time(grid_coord=True)
class Wind(GriddedSkeleton):
    pass
