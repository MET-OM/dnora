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
from dnora.read.ds_read_functions import basic_xarray_read
from dnora.read.product_readers import ProductReader
from dnora.read.product_configuration import ProductConfiguration, get_constant_url
import pandas as pd
import re
from dnora.process.gridded import FillNaNs
import geo_parameters as gp


def get_norkyst800_urls(folder: str, filename: str, file_times: list[str], **kwargs):
    """This is passed to the read_ds_list. We need it because the folder and filename that makes up th URL changes is time"""
    urls = []

    for file_time in file_times:
        if file_time >= pd.Timestamp("2018-01-01T00:00"):
            remote_folder = folder
        else:
            remote_folder = re.sub(
                "fou-hi/norkyst800m-1h", "sea/norkyst800mv0_1h/", folder
            )
        urls.append(get_url(remote_folder, filename, file_time))
    return urls


def pre_process_norkyst800_ds(ds: xr.Dataset) -> xr.Dataset:
    """Chooses surface level in ds"""
    return ds.isel(depth=0)


class NorKyst800(ProductReader):
    """Reads ocean_current data of the NorKyst800 archieve directly from MET Norways servers.

    NorKyst-800 (Norwegian Coast 800m) is a numerical, high-resolution, ocean modelling
    system covering the Norwegian Coast.

    Albretsen, J., Sperrevik, A.K., Staalstrøm, A., Sandvik, A.D., Vikebø, F., Asplin, L., 2011.
    NorKyst-800 Rapport nr. 1: Brukermanual og tekniske beskrivelser. NorKyst-800 Report
    No. 1: User Manual and technical descriptions.
    """

    product_configuration = ProductConfiguration(
        filename="NorKyst-800m_ZDEPTHS_his.an.%Y%m%d00.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h",
        },
        ds_creator_function=partial(
            ds_fimex_read,
            resolution_in_km=0.8,
            data_vars=["u", "v"],
        ),
        default_data_source=DataSource.REMOTE,
        url_function=get_norkyst800_urls,
        ds_aliases={"u": gp.ocean.XCurrent, "v": gp.ocean.YCurrent},
        ds_pre_processor=(lambda ds: ds.isel(depth=0)),
    )

    file_structure = FileStructure(
        stride=24,
        hours_per_file=24,
    )

    def post_processing(self):
        return FillNaNs(0)


class NorFjords160(DataReader):
    """ """

    _default_filename = "norfjords_160m_his_%Y%m%d01_surface_interp.nc"
    _default_folders = {
        DataSource.INTERNAL: "SWAN/Bjornafjorden2/ROMS/",
    }

    def __init__(
        self,
        stride: int = 24,
        hours_per_file: int = 24,
        last_file: str = "",
        lead_time: int = 0,
        offset: int = 1,
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

    def default_data_source(self) -> DataSource:
        return DataSource.INTERNAL

    def __call__(
        self,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        folder: str,
        filename: str,
        expansion_factor: float = 1.2,
        program: str = "pyfimex",
        **kwargs,
    ):
        """Reads in all grid points between the given times and at for the given indeces"""

        folder = self._folder(folder, source)
        filename = self._filename(filename, source)

        start_times, end_times, file_times = self.file_structure.create_time_stamps(
            start_time, end_time
        )
        setup_temp_dir(DnoraDataType.CURRENT, self.name())
        # Define area to search in
        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        lon, lat = utils.grid.expand_area(
            grid.edges("lon"), grid.edges("lat"), expansion_factor
        )
        msg.process(f"Applying {program}")
        ds_creator_function = partial(
            basic_xarray_read,
        )

        current_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
            url_function=get_constant_url,
            hours_per_file=self.file_structure.hours_per_file,
            lead_time=self.file_structure.lead_time,
        )
        msg.plain("Merging xarrays (this might take a while)...")

        ds = xr.concat(current_list, dim="ocean_time")

        lons, lats = ds.lon.values[:, 0], ds.lat.values[0, :]
        lon_mask = np.logical_and(lons >= lon[0], lons <= lon[1])
        lat_mask = np.logical_and(lats >= lat[0], lats <= lat[1])
        coord_dict = {
            "time": ds.ocean_time.data,
            "lon": lons[lon_mask],
            "lat": lats[lat_mask],
        }

        u = ds.u.fillna(0).data[:, :, lon_mask, :]
        u = u[:, :, :, lat_mask]
        v = ds.v.fillna(0).data[:, :, lon_mask, :]
        v = v[:, :, :, lat_mask]

        data_dict = {
            "u": (u, ["time", "lon", "lat"]),
            "v": (v, ["time", "lon", "lat"]),
        }

        meta_dict = ds.attrs

        return coord_dict, data_dict, meta_dict
