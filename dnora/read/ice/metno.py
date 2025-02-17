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
import pandas as pd
from dnora.read.product_readers import ProductReader
from dnora.read.product_configuration import ProductConfiguration, get_constant_url
from dnora.read.fimex_functions import ds_fimex_read


class NORA3(ProductReader):
    """Reads wind data (from monthly files 'arome3km_1hr_YYYYMM.nc') of the NORA3 hindcast directly from MET Norways servers.

    The NORA3 HARMONIE-AROME high-resolution (ca 3 km) hindcast for the
    North Sea, the Norwegian Sea, and the Barents Sea.

    Haakenstad, H., Breivik, Ø., Furevik, B. R., Reistad, M., Bohlinger, P., &
    Aarnes, O. J. (2021). NORA3: A Nonhydrostatic High-Resolution Hindcast of
    the North Sea, the Norwegian Sea, and the Barents Sea,
    Journal of Applied Meteorology and Climatology, 60(10), 1443-1464,
    DOI: 10.1175/JAMC-D-21-0029.1
    """

    product_configuration = ProductConfiguration(
        filename="%Y%m%d_MyWam3km_hindcast.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/%Y/%m",
        },
        ds_creator_function=partial(
            ds_fimex_read,
            resolution_in_km=3,
            data_vars=["SIC", "SIT"],
        ),
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure(
        stride=24,
        hours_per_file=24,
    )


class NORA3old(DataReader):
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
            url_function=get_constant_url,
            hours_per_file=self.file_structure.hours_per_file,
            lead_time=self.file_structure.lead_time,
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


# def get_barents25_urls(folder, filename, file_times, **kwargs):
#     """Remote folder and file name changes December 2022"""
#     urls = []

#     for file_time in file_times:
#         # if file_time > pd.Timestamp("2023-03-01 00:00"):
#         h6 = np.floor(file_time.hour / 6) * 6
#         remote_folder = re.sub("#6HSTAMP#", f"{h6:02.0f}", folder)
#         remote_filename = re.sub("#6HSTAMP#", f"{h6:02.0f}", filename)
#         # else:  # Old files
#         #     # Sub out instead of setting directly, since this should only affect the remote folder and filename!
#         #     filename = re.sub(
#         #         "barents_sfc_%Y%m%dT%HZm%H.nc",
#         #         "Barents-2.5km_ZDEPTHS_his.an.%Y%m%d06.nc",
#         #         filename,
#         #     )
#         #     remote_folder = re.sub(
#         #         "fou-hi/barents_eps_surface/%Y/%m/%d/T#6HSTAMP#Z",
#         #         "barents25km_files/",
#         #         folder,
#         #     )

#         urls.append(get_url(remote_folder, remote_filename, file_time))
#     return urls


class Barents25(ProductReader):
    """Reads sea ice data of the Barents 2.5 km operational ocean model directly from MET Norways servers.

    Barents-2.5km is an operational data-assimilative coupled ocean and sea ice ensemble prediction model for
    the Barents Sea and Svalbard

    Röhrs, J., Gusdal, Y., Rikardsen, E. S. U., Durán Moro, M., Brændshøi, J., Kristensen, N. M.,
    Fritzner, S., Wang, K., Sperrevik, A. K., Idžanović, M., Lavergne, T., Debernard, J. B., and Christensen, K. H.:
    Barents-2.5km v2.0: an operational data-assimilative coupled ocean and sea ice ensemble prediction model
    for the Barents Sea and Svalbard, Geosci. Model Dev., 16, 5401-5426,
    https://doi.org/10.5194/gmd-16-5401-2023, 2023.
    """

    product_configuration = ProductConfiguration(
        filename="barents_sfc_%Y%m%dT[6]Zm00.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_surface/%Y/%m/%d/T[6]Z"
        },
        ds_creator_function=partial(
            ds_fimex_read,
            resolution_in_km=2.5,
            data_vars=["ice_concentration", "ice_thickness"],
        ),
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure(
        stride=6,
        hours_per_file=97,
    )


# class Barents25old(DataReader):
#     """Reads sea ice data of the Barents 2.5 km operational ocean model directly from MET Norways servers.

#     Barents-2.5km is an operational data-assimilative coupled ocean and sea ice ensemble prediction model for
#     the Barents Sea and Svalbard

#     Röhrs, J., Gusdal, Y., Rikardsen, E. S. U., Durán Moro, M., Brændshøi, J., Kristensen, N. M.,
#     Fritzner, S., Wang, K., Sperrevik, A. K., Idžanović, M., Lavergne, T., Debernard, J. B., and Christensen, K. H.:
#     Barents-2.5km v2.0: an operational data-assimilative coupled ocean and sea ice ensemble prediction model
#     for the Barents Sea and Svalbard, Geosci. Model Dev., 16, 5401-5426,
#     https://doi.org/10.5194/gmd-16-5401-2023, 2023.
#     """

#     _default_folders = {
#         DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_surface/%Y/%m/%d/T#6HSTAMP#Z"
#     }
#     _default_filename = "barents_sfc_%Y%m%dT%HZm%H.nc"

#     def default_data_source(self) -> DataSource:
#         return DataSource.REMOTE

#     def __init__(
#         self,
#         last_file: str = "",
#         lead_time: int = 0,
#     ):
#         """The data is currently in 6 hourly files. Do not change the default
#         setting unless you have a good reason to do so.
#         """

#         self.last_file = last_file
#         self.lead_time = lead_time
#         return

#     def __call__(
#         self,
#         grid: Grid,
#         start_time: str,
#         end_time: str,
#         source: DataSource,
#         folder: str,
#         filename: str,
#         expansion_factor: float = 1.2,
#         program: str = "pyfimex",
#         **kwargs,
#     ):
#         """Reads in Barents 2.5 km ice data between the given times and area"""
#         # if start_time > pd.Timestamp("2023-03-01 00:00"):
#         #     self.file_structure = FileStructure(
#         #         stride=6,
#         #         hours_per_file=67,
#         #         last_file=self.last_file,
#         #         lead_time=self.lead_time,
#         #     )
#         # else:
#         #     # Old version uses daily files with no overlap
#         self.file_structure = FileStructure(stride=24, hours_per_file=24, offset=6)

#         folder = self._folder(folder, source)
#         filename = self._filename(filename, source)

#         start_times, end_times, file_times = self.file_structure.create_time_stamps(
#             start_time, end_time
#         )
#         msg.info(
#             f"Getting ice forcing from Barents2.5km from {start_time} to {end_time}"
#         )
#         setup_temp_dir(DnoraDataType.ICE, self.name())
#         # Define area to search in
#         msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
#         lon, lat = utils.grid.expand_area(
#             grid.edges("lon"), grid.edges("lat"), expansion_factor
#         )

#         msg.process(f"Applying {program}")
#         ds_creator_function = partial(
#             ds_fimex_read,
#             lon=lon,
#             lat=lat,
#             resolution_in_km=2.5,
#             data_vars=["ice_concentration", "ice_thickness"],
#             data_type=DnoraDataType.ICE,
#             name=self.name(),
#             program=program,
#         )

#         ds_list = read_ds_list(
#             start_times,
#             end_times,
#             file_times,
#             folder,
#             filename,
#             ds_creator_function,
#             url_function=get_barents25_urls,
#             hours_per_file=self.file_structure.hours_per_file,
#         )

#         ice_forcing = xr.concat(ds_list, dim="time")
#         data_dict = {
#             "sic": ice_forcing.ice_concentration.data,
#             "sit": ice_forcing.ice_thickness.data,
#         }
#         coord_dict = {
#             "time": ice_forcing.time.values,
#             "lon": ice_forcing.X.values,
#             "lat": ice_forcing.Y.values,
#         }
#         meta_dict = ice_forcing.attrs

#         return coord_dict, data_dict, meta_dict
