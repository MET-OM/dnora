import copernicusmarine
from copy import copy
import numpy as np
import xarray as xr
from subprocess import call
import os, glob

# Import objects
from dnora.grid import Grid

# Import abstract classes
from dnora.read.abstract_readers import DataReader
from dnora.type_manager.data_sources import DataSource
from dnora.read.file_structure import FileStructure

# Import aux_funcsiliry functions
from dnora import msg
from dnora.aux_funcs import (
    get_url,
)
from dnora import utils

from dnora.type_manager.dnora_types import DnoraDataType
from dnora.read.ds_read_functions import read_ds_list, setup_temp_dir
from functools import partial
from dnora.read.fimex_functions import ds_fimex_read


class Global(DataReader):
    """The Operational Mercator global ocean analysis and forecast system at 1/12 degree is providing 10 days of 3D global
    ocean forecasts updated daily. The time series is aggregated in time in order to reach a two full year's time series
    sliding window. This product includes daily and monthly mean files of temperature, salinity, currents, sea level,
    mixed layer depth and ice parameters from the top to the bottom over the global ocean. It also includes hourly
    mean surface fields for sea level height, temperature and currents. The global ocean output files are displayed
    with a 1/12 degree horizontal resolution with regular longitude/latitude equirectangular projection.
    50 vertical levels are ranging from 0 to 5500 meters. This product also delivers a special dataset for
    surface current which also includes wave and tidal drift called SMOC (Surface merged Ocean Current).

    DOI (product): https://doi.org/10.48670/moi-00016
    https://data.marine.copernicus.eu/product/GLOBAL_ANALYSISFORECAST_PHY_001_024/description
    """

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

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
        """Reads in all grid points between the given times and at for the given indeces"""
        setup_temp_dir(DnoraDataType.CURRENT, self.name())
        # Define area to search in
        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        lon, lat = utils.grid.expand_area(
            grid.edges("lon"), grid.edges("lat"), expansion_factor
        )

        cred_file = os.path.expanduser("~/.copernicusmarine/.copernicusmarine-credentials")
        if not os.path.isfile(cred_file):
            msg.advice(
                f"No credentials file {cred_file} was found. Login for the first time to create it."
            )
            copernicusmarine.login()
        copernicusmarine.subset(
            dataset_id="cmems_mod_glo_phy_anfc_0.083deg_PT1H-m",
            variables=["zos"],
            minimum_longitude=lon[0],
            maximum_longitude=lon[1],
            minimum_latitude=lat[0],
            maximum_latitude=lat[1],
            start_datetime=start_time.strftime("%Y-%m-%dT%H:%M:00"),
            end_datetime=end_time.strftime("%Y-%m-%dT%H:%M:00"),
            minimum_depth=0.49402499198913574,
            maximum_depth=0.49402499198913574,
            output_directory="dnora_current_temp",
            credentials_file=cred_file,
            force_download=True,
            output_filename=f"Global_CMEMS_temp.nc",
        )

        ds = xr.open_dataset(f"dnora_current_temp/Global_CMEMS_temp.nc")

        coord_dict = {
            "time": ds.time.data,
            "lon": ds.longitude.data,
            "lat": ds.latitude.data,
        }
        data_dict = {"eta": ds.zos.values[:, 0, :, :]}
        meta_dict = ds.attrs

        return coord_dict, data_dict, meta_dict


# copernicusmarine.login()
