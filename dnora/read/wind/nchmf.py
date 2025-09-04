import numpy as np
import xarray as xr
from subprocess import call
import os
import pandas as pd
from functools import partial
from dnora.read.file_structure import FileStructure
import re

# Import objects
from dnora.grid import Grid

# Import aux_funcsiliry functions
from dnora import msg
from dnora import utils

from dnora.type_manager.data_sources import DataSource
from dnora.type_manager.dnora_types import DnoraDataType
from dnora.read.abstract_readers import DataReader
from dnora.read.ds_read_functions import read_ds_list, setup_temp_dir

from dnora.read.product_readers import ProductReader
from dnora.read.product_configuration import ProductConfiguration
from dnora.read.ds_read_functions import ftp_read
from dnora.read.file_structure import FileStructure
from dnora.read.depreciation_decorator import deprecated_class_call


@deprecated_class_call("NCHMF", "nchmf", "wind")
class ECMWF(ProductReader):
    """Connects to MET Norways ftp server and downloads the ECMWF operational atmospheric data that covers Vietnam"""

    product_configuration = ProductConfiguration(
        filename="vietnam_atmos_%Y%m%d_00.nc",
        default_folders={},
        ds_creator_function=ftp_read,
        data_vars=["x_wind_10m", "y_wind_10m"],
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure(stride=24, hours_per_file=160)


def ds_xarray_read(
    start_time: pd.Timestamp,
    end_time: pd.Timestamp,
    url: str,
    lon: np.ndarray,
    lat: np.ndarray,
    data_vars: list[str],
):
    with xr.open_dataset(url) as f:
        ds = f.sel(
            lon=slice(lon[0], lon[1]),
            lat=slice(lat[0], lat[1]),
        )[data_vars]
        ds["time"] = start_time + pd.to_timedelta(ds.time.values)
        ds = ds.sel(time=slice(start_time, end_time))

    return ds


@deprecated_class_call("NCHMF", "nchmf", "wind")
class Oper(DataReader):
    """ """

    _default_filename = "wind%Y%m%d%H.nc"

    def default_data_source(self) -> DataSource:
        return DataSource.LOCAL

    def __init__(
        self,
        stride: int = 12,
        hours_per_file: int = 240,
        offset: int = 0,
        last_file: str = "",
        lead_time: int = 0,
    ):
        """The data is currently in daily files. Do not change the default
        setting unless you have a good reason to do so.
        """

        self.file_structure = FileStructure(
            stride=stride,
            hours_per_file=hours_per_file,
            last_file=last_file,
            lead_time=lead_time,
            offset=offset,
        )
        return

    def __call__(
        self,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        folder: str,
        filename: str,
        expansion_factor: float = 1.2,
        **kwargs,
    ):
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        folder = self._folder(folder, source)
        filename = self._filename(filename, source)

        start_times, end_times, file_times = self.file_structure.create_time_stamps(
            start_time, end_time
        )

        setup_temp_dir(DnoraDataType.WIND, self.name())
        # Define area to search in
        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        lon, lat = utils.grid.expand_area(
            grid.edges("lon"), grid.edges("lat"), expansion_factor, dlon=0.25, dlat=0.25
        )

        ds_creator_function = partial(
            ds_xarray_read,
            lon=lon,
            lat=lat,
            data_vars=["u10m", "v10m"],
        )
        wind_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
        )

        wind_forcing = xr.concat(wind_list, dim="time")

        data_dict = {"u": wind_forcing.u10m.data, "v": wind_forcing.v10m.data}
        wind_forcing.time.data
        coord_dict = {
            "time": wind_forcing.time.values,
            "lon": wind_forcing.lon.data,
            "lat": wind_forcing.lat.data,
        }
        meta_dict = wind_forcing.attrs

        return coord_dict, data_dict, meta_dict
