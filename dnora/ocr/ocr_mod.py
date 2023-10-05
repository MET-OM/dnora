from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_datavar, add_time

import pandas as pd
import numpy as np
from ..aux_funcs import speed_dir_from_u_v


@add_datavar(name="v", default_value=0.0)
@add_datavar(name="u", default_value=0.0)
@add_time(grid_coord=True)
class OceanCurrent(GriddedSkeleton):
    def __init__(
        self,
        x=None,
        y=None,
        lon=None,
        lat=None,
        time=pd.date_range("1990-01-01 00:00", "1990-01-01 01:00", freq="H"),
        name="LonelyOceanCurrent",
        **kwargs
    ):
        if np.all([a is None for a in [x, y, lon, lat]]):
            x, y = 0, 0
        super().__init__(x=x, y=y, lon=lon, lat=lat, name=name, time=time, **kwargs)

    def magnitude(self):
        ws, __ = speed_dir_from_u_v(self.u(), self.v())
        return ws

    def direction(self):
        __, wdir = speed_dir_from_u_v(self.u(), self.v())
        return wdir
