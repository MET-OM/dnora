from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_datavar, add_time

import pandas as pd
import numpy as np


@add_datavar(name="waterlevel", default_value=0.0)
@add_time(grid_coord=True)
class WaterLevel(GriddedSkeleton):
    def __init__(
        self,
        x=None,
        y=None,
        lon=None,
        lat=None,
        time=pd.date_range("1990-01-01 00:00", "1990-01-01 01:00", freq="H"),
        name="LonelyWaterLevel",
        **kwargs
    ):
        if np.all([a is None for a in [x, y, lon, lat]]):
            x, y = 0, 0
        super().__init__(x=x, y=y, lon=lon, lat=lat, name=name, time=time, **kwargs)
