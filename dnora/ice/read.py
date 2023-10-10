from abc import ABC, abstractmethod
import xarray as xr
import pandas as pd
import numpy as np

from ..grd.grd_mod import Grid
from .. import aux_funcs
from .. import msg


class IceForcingReader(ABC):
    """Reads ice forcing data from some source and provide it to the object.

    The area is defined from the Grid object that is passed.
    """

    @abstractmethod
    def __call__(
        self, grid: Grid, start_time: str, end_time: str, expansion_factor: float
    ):
        return time, concentration, thickness, lon, lat, x, y, metadata

    def name(self) -> str:
        return type(self).__name__


class DnoraNc(IceForcingReader):
    def __init__(self, files: str) -> None:
        self.files = files

    def __call__(
        self, grid, start_time, end_time, expansion_factor: float = 1.2, **kwargs
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
            f"Getting ice forcing from cached netcdf (e.g. {self.files[0]}) from {start_time} to {end_time}"
        )

        # These files might get deleted, so we don't want to use dask for a lazy load
        ds = xr.open_mfdataset(self.files, preprocess=_crop, data_vars="minimal")
        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds)

        return (
            ds.time.values,
            ds.ice_concentration.values,
            ds.ice_thickness.values,
            lon,
            lat,
            x,
            y,
            ds.attrs,
        )


class ConstantIceForcing(IceForcingReader):
    def __init__(
        self, concentration: float = 0.4, thickness: float = 0.5, metadata: dict = None
    ):
        self.concentration = concentration
        self.thickness = thickness
        self.metadata = metadata

    def __call__(self, grid, start_time, end_time, **kwargs):
        time = pd.date_range(start=start_time, end=end_time, freq="H").values

        lon = grid.lon(strict=True)
        lat = grid.lat(strict=True)
        x = grid.x(strict=True)
        y = grid.y(strict=True)

        concentration = np.full((len(time), grid.ny(), grid.nx()), self.concentration)
        thickness = np.full((len(time), grid.ny(), grid.nx()), self.thickness)
        metadata = {"metadata": "this is a constant ice forcing"}

        return time, concentration, thickness, lon, lat, x, y, metadata
