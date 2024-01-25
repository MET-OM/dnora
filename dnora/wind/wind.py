from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_time, add_datavar
from dnora.aux_funcs import speed_dir_from_u_v

from dnora.metaparameter.parameters import XWind, YWind


@add_datavar(name="v", default_value=0.0)
@add_datavar(name="u", default_value=0.0)
@add_time(grid_coord=True)
class Wind(GriddedSkeleton):
    meta_dict = {"u": XWind, "v": YWind}

    def magnitude(self):
        ws, __ = speed_dir_from_u_v(self.u(), self.v())
        return ws

    def direction(self):
        __, wdir = speed_dir_from_u_v(self.u(), self.v())
        return wdir
