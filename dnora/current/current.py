from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_datavar, add_time

from dnora.aux_funcs import speed_dir_from_u_v

import geo_parameters as gp


@add_datavar(name=gp.ocean.YCurrent("v"), default_value=0.0)
@add_datavar(name=gp.ocean.XCurrent("u"), default_value=0.0)
@add_time(grid_coord=True)
class Current(GriddedSkeleton):
    pass
    # def magnitude(self):
    #     ws, __ = speed_dir_from_u_v(self.u(), self.v())
    #     return ws

    # def direction(self):
    #     __, wdir = speed_dir_from_u_v(self.u(), self.v())
    #     return wdir
