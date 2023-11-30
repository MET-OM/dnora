from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_time, add_datavar
from ..aux_funcs import speed_dir_from_u_v


@add_datavar(name="v", default_value=0.0)
@add_datavar(name="u", default_value=0.0)
@add_time(grid_coord=True)
class Wind(GriddedSkeleton):
    def magnitude(self):
        ws, __ = speed_dir_from_u_v(self.u(), self.v())
        return ws

    def direction(self):
        __, wdir = speed_dir_from_u_v(self.u(), self.v())
        return wdir
