from abc import ABC, abstractmethod
import xarray as xr

# Import objects
from ..grid.grid import Grid
from .. import aux_funcs
from .. import msg
import pandas as pd
import numpy as np
from geo_skeletons import PointSkeleton
from ..data_sources import DataSource


class CurrentReader(ABC):
    """Reads forcing data from some source and provide it to the object.

    The area is defined from the Grid object that is passed.
    """

    @abstractmethod
    def __call__(
        self, grid: Grid, start_time: str, end_time: str, source: str, **kwargs
    ):
        """Reads in the forcing witih grid and between start_time and end_time.

        The variables needed to be returned are:

        time:   Time stamps as numpy.datetime64 array
        u:      west-to-east velocity [lon, lat, time] as numpy array
        v:      south-to-north velocity [lon, lat, time] as numpy array
        lon:    Longitude vector as numpy array (None if Cartesian)
        lat:    Latitude vector as numpy array (None if Cartesian)
        x:      Longitude vector as numpy array (None if Spherical)
        y:      Latitude vector as numpy array (None if Spherical)
        metadata: dict{key, value} will be set as attributes of the xr.Dataset
        """

        return time, u, v, lon, lat, x, y, metadata

    def name(self) -> str:
        return type(self).__name__

    def default_data_source(self) -> DataSource:
        return DataSource.UNDEFINED


class ConstantOceanCurrent(CurrentReader):
    def __init__(self, u: float = 1, v: float = 2, metadata: dict = None):
        self.u = u
        self.v = v
        self.metadata = metadata

    def __call__(self, grid, start_time, end_time, source: DataSource, **kwargs):
        time = pd.date_range(start=start_time, end=end_time, freq="H").values

        lon = grid.lon(strict=True)
        lat = grid.lat(strict=True)
        x = grid.x(strict=True)
        y = grid.y(strict=True)

        u = np.full((len(time), grid.ny(), grid.nx()), self.u)
        v = np.full((len(time), grid.ny(), grid.nx()), self.v)
        metadata = {"metadata": "this is a constant forcing"}

        return time, u, v, lon, lat, x, y, metadata


class DnoraNc(CurrentReader):
    def __init__(self, files: str) -> None:
        self.files = files

    def __call__(
        self,
        grid,
        start_time,
        end_time,
        source: DataSource,
        expansion_factor: float = 1.2,
        **kwargs,
    ):
        def _crop(ds):
            if lon is not None:
                return ds.sel(
                    time=slice(start_time, end_time),
                    lon=slice(lon[0], lon[1]),
                    lat=slice(lat[0], lat[1]),
                )
            else:
                return ds.sel(
                    time=slice(start_time, end_time),
                    x=slice(x[0], x[1]),
                    y=slice(y[0], y[1]),
                )

        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        ds0 = xr.open_dataset(self.files[0])
        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds0)

        if lon is not None:
            lon, lat = aux_funcs.expand_area(
                grid.edges("lon"), grid.edges("lat"), expansion_factor
            )
        else:
            x, y = aux_funcs.expand_area(
                grid.edges("x"), grid.edges("y"), expansion_factor
            )
        msg.info(
            f"Getting wind forcing from cached netcdf (e.g. {self.files[0]}) from {start_time} to {end_time}"
        )

        # These files might get deleted, so we don't want to use dask for a lazy load
        ds = xr.open_mfdataset(self.files, preprocess=_crop, data_vars="minimal")
        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds)

        return ds.time.values, ds.u.values, ds.v.values, lon, lat, x, y, ds.attrs
