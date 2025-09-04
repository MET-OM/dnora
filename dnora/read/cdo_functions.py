import xarray as xr
import numpy as np
import pandas as pd
from dnora.type_manager.dnora_types import DnoraDataType
from dnora.utils.distance import lon_in_km
from subprocess import call
from typing import Callable
from cdo import Cdo
from geo_skeletons import GriddedSkeleton

from cdo import Cdo


def ds_cdo_read(
    start_time,
    end_time,
    url,
    lon,
    lat,
    resolution_in_km: float,
    data_vars: list[str],
    data_type: DnoraDataType,
    name: str,
) -> xr.Dataset:
    # create grid_definition_file
    def write_grid_definition(Grid: GriddedSkeleton):
        with open("cdo_grid.txt", "w") as f:
            f.write(
                f"""gridtype = lonlat
    xsize    = {grid.nx()}
    ysize    = {grid.ny()}
    xfirst   = {grid.lon()[0]}
    xinc     = {grid.dlon()}
    yfirst   = {grid.lat()[0]}
    yinc     = {grid.dlat()}
    """
            )

    grid = GriddedSkeleton(lon=lon, lat=lat)
    grid.set_spacing(dm=resolution_in_km * 1000)

    write_grid_definition(grid)
    cdo = Cdo()
    breakpoint()
    cdo.remapnn("cdo_grid.txt", input=url, output="cdo_out.nc")
    ds = xr.open_dataset("cdo_out.nc")
