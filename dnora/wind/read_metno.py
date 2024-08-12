from copy import copy
import numpy as np
import xarray as xr
from subprocess import call
import os, glob
import pandas as pd
from functools import partial
from dnora.readers.file_structure import FileStructure
import re

# Import objects
from dnora.grid import Grid

# Import aux_funcsiliry functions
from dnora import msg
from dnora.aux_funcs import (
    create_time_stamps,
    u_v_from_speed_dir,
    expand_area,
    lon_in_km,
    get_url,
    create_monthly_stamps,
)

from dnora.dnora_type_manager.data_sources import DataSource
from dnora.dnora_type_manager.dnora_types import DnoraDataType
from dnora.readers.abstract_readers import DataReader
from dnora.readers.fimex_functions import ds_fimex_read
from dnora.readers.ds_read_functions import read_ds_list, setup_temp_dir
import calendar


class NORA3(DataReader):
    """Reads wind data (from monthly files 'arome3km_1hr_YYYYMM.nc') of the NORA3 hindcast directly from MET Norways servers.

    The NORA3 HARMONIE-AROME high-resolution (ca 3 km) hindcast for the
    North Sea, the Norwegian Sea, and the Barents Sea.

    Haakenstad, H., Breivik, Ø., Furevik, B. R., Reistad, M., Bohlinger, P., &
    Aarnes, O. J. (2021). NORA3: A Nonhydrostatic High-Resolution Hindcast of
    the North Sea, the Norwegian Sea, and the Barents Sea,
    Journal of Applied Meteorology and Climatology, 60(10), 1443-1464,
    DOI: 10.1175/JAMC-D-21-0029.1
    """

    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/nora3_subset_atmos/atm_hourly",
        DataSource.INTERNAL: "NORA3/atmosphere/atm_hourly",
    }
    _default_filename = "arome3km_1hr_%Y%m.nc"

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
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        folder = self._folder(folder, source)
        filename = self._filename(filename, source)

        start_times, end_times = create_monthly_stamps(start_time, end_time)
        file_times = start_times

        setup_temp_dir(DnoraDataType.WIND, self.name())
        # Define area to search in
        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

        msg.process(f"Applying {program}")
        ds_creator_function = partial(
            ds_fimex_read,
            lon=lon,
            lat=lat,
            resolution_in_km=3,
            data_vars=["wind_speed", "wind_direction"],
            data_type=DnoraDataType.WIND,
            name=self.name(),
            program=program,
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
        # Go to u and v components
        u, v = u_v_from_speed_dir(wind_forcing.wind_speed, wind_forcing.wind_direction)

        data_dict = {"u": u.fillna(0).values, "v": v.fillna(0).values}
        coord_dict = {
            "time": wind_forcing.time.values,
            "lon": wind_forcing.x.values,
            "lat": wind_forcing.y.values,
        }
        meta_dict = wind_forcing.attrs

        return coord_dict, data_dict, meta_dict


class MyWave3km(DataReader):
    """Reads wind data from the MyWave 3km hindcast directly from MET Norways
    servers. You should probably use MetNo_NORA3 because:

    The wind data is from NORA3 (see the MetNo_NORA3 reader), is taken
    from the wave model output. This means that model land points have no data.
    """

    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/%Y/%m",
    }
    _default_filename = "%Y%m%d_MyWam3km_hindcast.nc"

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

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
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        folder = self._folder(folder, source)
        filename = self._filename(filename, source)

        start_times, end_times, file_times = self.file_structure.create_time_stamps(
            start_time, end_time
        )

        setup_temp_dir(DnoraDataType.WIND, self.name())
        # Define area to search in
        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

        msg.process(f"Applying {program}")
        ds_creator_function = partial(
            ds_fimex_read,
            lon=lon,
            lat=lat,
            resolution_in_km=3,
            data_vars=["ff", "dd"],
            data_type=DnoraDataType.WIND,
            name=self.name(),
            program=program,
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

        # Go to u and v components
        u, v = u_v_from_speed_dir(wind_forcing.ff, wind_forcing.dd)  # factor 1000

        u = -1 * u.fillna(0)  # *-1 due to ocean convection in WAM!!!
        v = -1 * v.fillna(0)  # *-1 due to ocean convection in WAM!!!

        data_dict = {"u": u.data, "v": v.data}
        coord_dict = {
            "time": wind_forcing.time.data,
            "lon": wind_forcing.rlon.data,
            "lat": wind_forcing.rlat.data,
        }
        meta_dict = wind_forcing.attrs

        return coord_dict, data_dict, meta_dict


def get_meps_urls(folder, filename, file_times):
    """This is passed to the read_ds_list. We need it because the folder and filename that makes up th URL changes is time"""
    urls = []

    for file_time in file_times:
        if file_time >= pd.Timestamp("2020-02-04T12:00"):
            remote_filename = filename
        else:
            remote_filename = re.sub("det", "subset", filename)
        if file_time >= pd.Timestamp("2020-01-01T00:00"):
            remote_folder = folder
        else:
            remote_folder = re.sub("meps25epsarchive", "mepsoldarchive", folder)

        urls.append(get_url(remote_folder, remote_filename, file_time))
    return urls


def meps_extra_fimex_commands(start_time, end_time, url) -> list[str]:
    """Determines the possible extra fimex commands needed to process the MEPS netcdf.
    Repend on the url ('det'/'subset')

    start_time and end_time accepted because of standard convention"""
    if "subset" in url:
        return [
            "--extract.reduceDimension.name=ensemble_member",
            "--extract.reduceDimension.start=1",
            "--extract.reduceDimension.end=1",
        ]
    else:
        return []


class MEPS(DataReader):
    """Reads wind data from MET Norways MEPS forecast.

    The data is from a 2.5 km AROME model.
    """

    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/meps25epsarchive/%Y/%m/%d",
    }
    _default_filename = f"meps_det_2_5km_%Y%m%dT%HZ.nc"

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def __init__(
        self,
        stride: int = 6,
        hours_per_file: int = 67,
        lead_time: int = 0,
        last_file: str = "",
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
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        start_times, end_times, file_times = self.file_structure.create_time_stamps(
            start_time, end_time
        )

        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

        setup_temp_dir(DnoraDataType.WIND, self.name())

        # Check weather to use 'det' or 'subset' files
        folder = self._folder(folder, source)
        filename = self._filename(filename, source)

        setup_temp_dir(DnoraDataType.WIND, self.name())
        # Define area to search in
        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

        msg.process(f"Applying {program}")
        ds_creator_function = partial(
            ds_fimex_read,
            lon=lon,
            lat=lat,
            resolution_in_km=2.5,
            data_vars=["x_wind_10m", "y_wind_10m"],
            data_type=DnoraDataType.WIND,
            name=self.name(),
            program=program,
            extra_commands=meps_extra_fimex_commands,
        )
        wind_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
            url_function=get_meps_urls,
            hours_per_file=self.file_structure.hours_per_file,
        )

        wind_forcing = xr.concat(wind_list, dim="time", coords="minimal")

        data_dict = {
            "u": wind_forcing.x_wind_10m.values,
            "v": wind_forcing.y_wind_10m.values,
        }
        coord_dict = {
            "time": wind_forcing.time.values,
            "lon": wind_forcing.x.values,
            "lat": wind_forcing.y.values,
        }
        meta_dict = wind_forcing.attrs

        return coord_dict, data_dict, meta_dict


# class NORA3_fp(DataReader):
#     """Reads wind data of the NORA3 hindcast directly from MET Norways servers.

#     The NORA3 HARMONIE-AROME high-resolution (ca 3 km) hindcast for the
#     North Sea, the Norwegian Sea, and the Barents Sea.

#     Haakenstad, H., Breivik, Ø., Furevik, B. R., Reistad, M., Bohlinger, P., &
#     Aarnes, O. J. (2021). NORA3: A Nonhydrostatic High-Resolution Hindcast of
#     the North Sea, the Norwegian Sea, and the Barents Sea,
#     Journal of Applied Meteorology and Climatology, 60(10), 1443-1464,
#     DOI: 10.1175/JAMC-D-21-0029.1
#     """

#     def default_data_source(self) -> DataSource:
#         return DataSource.REMOTE

#     def __init__(
#         self,
#         stride: int = 1,
#         hours_per_file: int = 1,
#         last_file: str = "",
#         lead_time: int = 4,
#     ):
#         """The data is currently in hourly files. Do not change the default
#         setting unless you have a good reason to do so.
#         """

#         self.file_structure = FileStructure(
#             stride=stride,
#             hours_per_file=hours_per_file,
#             last_file=last_file,
#             lead_time=lead_time,
#         )

#     def _folder(
#         self,
#         folder: str,
#         time_stamp_file,
#         source: DataSource,
#     ) -> str:
#         if source == DataSource.REMOTE:
#             h0 = int(time_stamp_file.hour) % 6
#             subfolder = time_stamp_file.strftime("%Y/%m/%d/") + (
#                 time_stamp_file - np.timedelta64(h0, "h")
#             ).strftime("%H")
#             folder = "https://thredds.met.no/thredds/dodsC/nora3/" + subfolder

#         return folder

#     def _filename(self, filename: str, time_stamp_file, time_stamp, first_ind) -> str:
#         h0 = int(time_stamp_file.hour) % 6
#         ind = int((time_stamp.hour - first_ind) % 6) + first_ind
#         filename = (
#             time_stamp_file.strftime("fc%Y%m%d")
#             + (time_stamp_file - np.timedelta64(h0, "h")).strftime("%H")
#             + f"_{ind:03d}_fp.nc"
#         )

#         return filename

#     def __call__(
#         self,
#         grid: Grid,
#         start_time: str,
#         end_time: str,
#         source: DataSource,
#         folder: str,
#         filename: str,
#         expansion_factor: float = 1.2,
#         **kwargs,
#     ):
#         """Reads boundary spectra between given times and given area around
#         the Grid object."""

#         self.start_time = start_time
#         self.end_time = end_time

#         start_times, end_times, file_times = self.file_structure.create_time_stamps(
#             start_time,
#             end_time,
#         )

#         msg.info(
#             f"Getting wind forcing from NORA3 from {self.start_time} to {self.end_time}"
#         )

#         setup_temp_dir(DnoraDataType.WIND, self.name())

#         msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
#         lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

#         # Set resolution to about 3 km
#         x_str, y_str = create_fimex_xy_strings(lon, lat, resolution_in_km=3)

#         wnd_list = []
#         for n in range(len(file_times)):
#             folder = self._folder(folder, file_times[n], source)
#             filename = self._filename(
#                 filename, file_times[n], start_times[n], first_ind=self.lead_time
#             )
#             url = get_url(folder, filename)

#             msg.from_file(url)
#             msg.plain(f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")

#             nc_fimex = f"dnora_wind_temp/NORA3_{n:04.0f}_fimextemp.nc"

#             fimex_command = [
#                 "fimex",
#                 "--input.file=" + url,
#                 "--interpolate.method=bilinear",
#                 "--interpolate.projString=+proj=latlong +ellps=sphere +a=6371000 +e=0",
#                 "--interpolate.xAxisValues=" + x_str + "",
#                 "--interpolate.yAxisValues=" + y_str + "",
#                 "--interpolate.xAxisUnit=degree",
#                 "--interpolate.yAxisUnit=degree",
#                 "--process.rotateVector.all",
#                 "--extract.selectVariables=wind_speed",
#                 "--extract.selectVariables=wind_direction",
#                 "--extract.reduceTime.start="
#                 + start_times[n].strftime("%Y-%m-%dT%H:%M:%S"),
#                 "--extract.reduceTime.end="
#                 + end_times[n].strftime("%Y-%m-%dT%H:%M:%S"),
#                 "--process.rotateVector.direction=latlon",
#                 "--output.file=" + nc_fimex,
#             ]

#             # read_success = False
#             # for ct in range(5):  # try 6 times
#             #     try:
#             #         call(fimex_command)
#             #         read_success = True
#             #     except:
#             #         print(f'......Retry {ct}.....')
#             #         time.sleep(10)  # wait for 10 seconds before re-trying
#             #
#             # # Don't want to catch the last execption
#             # if not read_success:
#             call(fimex_command)
#             wnd_list.append(xr.open_dataset(nc_fimex).squeeze())

#         wind_forcing = xr.concat(wnd_list, dim="time")

#         # Fimex has already rotated the longitudes and latitudes, so calling them rlon/rlat is now incorrect
#         #        wind_forcing = wind_forcing.rename_dims({'y': 'lat', 'x': 'lon'})
#         # wind_forcing = wind_forcing.rename_vars({'y': 'lat', 'x': 'lon'})

#         # Go to u and v components
#         u, v = u_v_from_speed_dir(wind_forcing.wind_speed, wind_forcing.wind_direction)

#         data_dict = {"u": u.fillna(0).data, "v": v.fillna(0).data}
#         coord_dict = {
#             "time": wind_forcing.time.data,
#             "lon": wind_forcing.x.data,
#             "lat": wind_forcing.y.data,
#         }
#         meta_dict = wind_forcing.attrs

#         return coord_dict, data_dict, meta_dict
