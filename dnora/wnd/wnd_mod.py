from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_time, add_datavar
import numpy as np
import pandas as pd
from ..aux_funcs import speed_dir_from_u_v


@add_datavar(name="v", default_value=0.0)
@add_datavar(name="u", default_value=0.0)
@add_time(grid_coord=True)
class Forcing(GriddedSkeleton):
    def __init__(
        self,
        x=None,
        y=None,
        lon=None,
        lat=None,
        time=pd.date_range("1990-01-01 00:00", "1990-01-01 01:00", freq="H"),
        name="LonelyForcing",
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

    # def __str__(self) -> str:
    #     """Prints status of forcing."""

    #     msg.header(self, f"Status of forcing {self.name}")
    #     if not self.time() == (None, None):
    #         msg.plain(f"Contains data for {self.time()[0]} - {self.time()[-1]}")
    #         msg.plain(f"\t dt={self.dt()} hours, i.e. ({self.nt()} time steps)")
    #         msg.plain(f"Data covers: lon: {min(self.lon())} - {max(self.lon())}, lat: {min(self.lat())} - {max(self.lat())}")
    #         msg.plain(f"\t {self.ny()}x{self.nx()} grid points, dlon/dlat={np.mean(np.diff(self.lon()))}/{np.mean(np.diff(self.lat()))}")
    #     if len(self._history) > 0:
    #         msg.blank()
    #         msg.plain("Object has the following history:")
    #         for obj in self._history:
    #             msg.process(f"{obj.__class__.__bases__[0].__name__}: {type(obj).__name__}")
    #     #msg.print_line()
    #     #msg.plain("The Forcing is for the following Grid:")
    #     #print(self.grid())

    #     msg.print_line()

    #     return ''
