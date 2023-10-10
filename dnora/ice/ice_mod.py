from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_datavar, add_time
import numpy as np
import pandas as pd


@add_datavar(name="thickness", default_value=0.0)
@add_datavar(name="concentration", default_value=0.0)
@add_time(grid_coord=True)
class IceForcing(GriddedSkeleton):
    def __init__(
        self,
        x=None,
        y=None,
        lon=None,
        lat=None,
        time=pd.date_range("1990-01-01 00:00", "1990-01-01 01:00", freq="H"),
        name="LonelyIceForcing",
        **kwargs
    ):
        if np.all([a is None for a in [x, y, lon, lat]]):
            x, y = 0, 0
        super().__init__(x=x, y=y, lon=lon, lat=lat, name=name, time=time, **kwargs)
