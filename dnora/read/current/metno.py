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


class NorKyst800(DataReader):
    """Reads ocean_current data of the NorKyst800 archieve directly from MET Norways servers.

    NorKyst-800 (Norwegian Coast 800m) is a numerical, high-resolution, ocean modelling
    system covering the Norwegian Coast.

    Albretsen, J., Sperrevik, A.K., Staalstrøm, A., Sandvik, A.D., Vikebø, F., Asplin, L., 2011.
    NorKyst-800 Rapport nr. 1 : Brukermanual og tekniske beskrivelser. NorKyst-800 Report
    No. 1 : User Manual and technical descriptions.
    """

    _default_filename = "NorKyst-800m_ZDEPTHS_his.an.%Y%m%d00.nc"

    def __init__(
        self,
        stride: int = 24,
        hours_per_file: int = 24,
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
        )

        return

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
        program: str = "pyfimex",
        **kwargs,
    ):
        """Reads in all grid points between the given times and at for the given indeces"""
        if start_time.year < 2018:
            self._default_folders = {
                DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/sea/norkyst800mv0_1h/",
            }
        else:
            self._default_folders = {
                DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h",
            }

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
            ds_fimex_read,
            lon=lon,
            lat=lat,
            resolution_in_km=0.8,
            data_vars=["u", "v"],
            data_type=DnoraDataType.CURRENT,
            name=self.name(),
            program=program,
        )
        current_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
        )

        ds = xr.concat(current_list, dim="time")
        # # Rename X/Y  to lon/lat
        # oceancurrent_forcing = oceancurrent_forcing.rename_dims(
        #     {"Y": "lat", "X": "lon"}
        # )
        # oceancurrent_forcing = oceancurrent_forcing.rename_vars(
        #     {"Y": "lat", "X": "lon"}
        # )
        # Select depth = 0 m
        ds = ds.sel(depth=0)

        # ds["u"] = ds["u"].fillna(0)
        # ds["v"] = ds["v"].fillna(0)

        data_dict = {"u": ds.u.fillna(0).data, "v": ds.v.fillna(0).data}
        coord_dict = {
            "time": ds.time.data,
            "lon": ds.X.data,
            "lat": ds.Y.data,
        }
        meta_dict = ds.attrs

        return coord_dict, data_dict, meta_dict


class NorFjords160(DataReader):
    """ """

    _default_filename = "norfjords_160m_his_%Y%m%d01_surface_interp.nc"
    _default_folders = {
        DataSource.INTERNAL: "fou/om/SWAN/Bjornafjorden2/ROMS/",
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
            start_times, end_times, file_times, folder, filename, ds_creator_function
        )
        msg.plain("Merging xarrays (this might take a while)...")
        ds = xr.concat(current_list, dim="time")

        lons, lats = ds.lon.values[:, 0], ds.lat.values[0, :]
        lon_mask = np.logical_and(lons >= lon[0], lons <= lon[1])
        lat_mask = np.logical_and(lats >= lat[0], lats <= lat[1])
        coord_dict = {
            "time": ds.time.data,
            "lon": lons[lon_mask],
            "lat": lats[lat_mask],
        }

        u = ds.u.fillna(0).values[:, :, lon_mask, :]
        u = u[:, :, :, lat_mask]
        v = ds.v.fillna(0).values[:, :, lon_mask, :]
        v = v[:, :, :, lat_mask]

        data_dict = {
            "u": (u, ["time", "lon", "lat"]),
            "v": (v, ["time", "lon", "lat"]),
        }

        meta_dict = ds.attrs

        return coord_dict, data_dict, meta_dict
