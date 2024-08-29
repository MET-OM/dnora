from dnora.type_manager.data_sources import DataSource
from dnora.type_manager.dnora_types import DnoraDataType
from dnora.read.abstract_readers import DataReader
from dnora.read.fimex_functions import ds_fimex_read
from dnora.read.ds_read_functions import read_ds_list, setup_temp_dir
from dnora.read.file_structure import FileStructure
from dnora.grid import Grid
from dnora import msg
from dnora import utils
from dnora.aux_funcs import get_url
from functools import partial
import xarray as xr
import re
import numpy as np


class NORA3(DataReader):
    """Reads ice data of the NORA3 wave hindcast directly from MET Norways servers.

    The NORA3 high-resolution (ca 3 km) hindcast for the Nordic Seas and the Arctic Ocean.

    Breivik, Ø., Carrasco, A., Haakenstad, H., Aarnes, O. J., Behrens, A., Bidlot, J.-R., et al. (2022).
    The impact of a reduced high-wind Charnock parameter on wave growth with application to the North Sea,
    the Norwegian Sea, and the Arctic Ocean.
    Journal of Geophysical Research: Oceans, 127, e2021JC018196. https://doi.org/10.1029/2021JC018196
    """

    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/%Y/%m"
    }
    _default_filename = "%Y%m%d_MyWam3km_hindcast.nc"

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
        """Reads in NORA3 ice data between the given times and area"""
        self.file_structure = FileStructure(
            stride=24,
            hours_per_file=24,
        )

        folder = self._folder(folder, source)
        filename = self._filename(filename, source)

        start_times, end_times, file_times = self.file_structure.create_time_stamps(
            start_time, end_time
        )
        msg.info(f"Getting ice forcing from NORA3 from {start_time} to {end_time}")
        setup_temp_dir(DnoraDataType.ICE, self.name())
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
            resolution_in_km=3,
            data_vars=["SIC", "SIT"],
            data_type=DnoraDataType.ICE,
            name=self.name(),
            program=program,
        )
        ds_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
        )

        ice_forcing = xr.concat(ds_list, dim="time")

        data_dict = {
            "sic": ice_forcing.SIC.data,
            "sit": ice_forcing.SIT.data,
        }
        coord_dict = {
            "time": ice_forcing.time.values,
            "lon": ice_forcing.rlon.values,
            "lat": ice_forcing.rlat.values,
        }
        meta_dict = ice_forcing.attrs

        return coord_dict, data_dict, meta_dict


def get_barents25_urls(folder, filename, file_times):
    urls = []

    for file_time in file_times:
        h6 = np.floor(file_time.hour / 6) * 6
        remote_folder = re.sub("#6HSTAMP#", f"{h6:02.0f}", folder)
        urls.append(get_url(remote_folder, filename, file_time))
    return urls


class Barents25(DataReader):
    """Reads sea ice data of the Barents 2.5 km operational ocean model directly from MET Norways servers.

    Barents-2.5km is an operational data-assimilative coupled ocean and sea ice ensemble prediction model for
    the Barents Sea and Svalbard

    Röhrs, J., Gusdal, Y., Rikardsen, E. S. U., Durán Moro, M., Brændshøi, J., Kristensen, N. M.,
    Fritzner, S., Wang, K., Sperrevik, A. K., Idžanović, M., Lavergne, T., Debernard, J. B., and Christensen, K. H.:
    Barents-2.5km v2.0: an operational data-assimilative coupled ocean and sea ice ensemble prediction model
    for the Barents Sea and Svalbard, Geosci. Model Dev., 16, 5401-5426,
    https://doi.org/10.5194/gmd-16-5401-2023, 2023.
    """

    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_surface/%Y/%m/%d/T#6HSTAMP#Z"
    }
    _default_filename = "barents_sfc_%Y%m%dT%HZm%H.nc"

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def __init__(
        self,
        stride: int = 6,
        hours_per_file: int = 67,
        last_file: str = "",
        lead_time: int = 0,
    ):
        """The data is currently in 6 hourly files. Do not change the default
        setting unless you have a good reason to do so.
        """

        self.file_structure = FileStructure(
            stride=stride,
            hours_per_file=hours_per_file,
            last_file=last_file,
            lead_time=lead_time,
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
        program: str = "pyfimex",
        **kwargs,
    ):
        """Reads in Barents 2.5 km ice data between the given times and area"""

        folder = self._folder(folder, source)
        filename = self._filename(filename, source)

        start_times, end_times, file_times = self.file_structure.create_time_stamps(
            start_time, end_time
        )
        msg.info(
            f"Getting ice forcing from Barents2.5km from {start_time} to {end_time}"
        )
        setup_temp_dir(DnoraDataType.ICE, self.name())
        # Define area to search in
        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

        msg.process(f"Applying {program}")
        ds_creator_function = partial(
            ds_fimex_read,
            lon=lon,
            lat=lat,
            resolution_in_km=2.5,
            data_vars=["ice_concentration", "ice_thickness"],
            data_type=DnoraDataType.ICE,
            name=self.name(),
            program=program,
        )
        ds_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
            url_function=get_barents25_urls,
            hours_per_file=self.file_structure.hours_per_file,
        )

        ice_forcing = xr.concat(ds_list, dim="time")
        data_dict = {
            "sic": ice_forcing.ice_concentration.data,
            "sit": ice_forcing.ice_thickness.data,
        }
        coord_dict = {
            "time": ice_forcing.time.values,
            "lon": ice_forcing.X.values,
            "lat": ice_forcing.Y.values,
        }
        meta_dict = ice_forcing.attrs

        return coord_dict, data_dict, meta_dict
