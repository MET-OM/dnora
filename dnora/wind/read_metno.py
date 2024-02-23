from abc import ABC, abstractmethod
from copy import copy
import numpy as np
import xarray as xr
from subprocess import call
import os, glob
import pandas as pd

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
    pyfimex,
    create_monthly_stamps,
)

from dnora.dnora_types import DataSource
from dnora.readers.abstract_readers import DataReader
import calendar


def create_fimex_xy_strings(
    lon: tuple[float, float], lat: tuple[float, float], resolution_in_km: float
) -> tuple[str, str]:
    # Set resolution
    dlat = resolution_in_km / 111
    mean_lon_in_km = (lon_in_km(lat[0]) + lon_in_km(lat[-1])) * 0.5
    dlon = resolution_in_km / mean_lon_in_km

    if len(np.unique(lon)) == 1:
        x_str = str(lon[0])
    else:
        if lon[1] < lon[0] + dlon:
            lon = (lon[0] - dlon, lon[0] + dlon)
        x_str = str(lon[0]) + "," + str(lon[0] + dlon) + ",...," + str(lon[-1])
    if len(np.unique(lat)) == 1:
        y_str = str(lat[0])
    else:
        if lat[1] < lat[0] + dlat:
            lat = (lat[0] - dlat, lat[0] + dlat)
        y_str = str(lat[0]) + "," + str(lat[0] + dlat) + ",...," + str(lat[-1])

    return x_str, y_str


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

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def _folder_filename(
        self, source: DataSource, folder: str, filename: str
    ) -> tuple[str]:
        if source == DataSource.REMOTE:
            folder = (
                "https://thredds.met.no/thredds/dodsC/nora3wavesubset_files/atm_hourly"
            )
        elif source == DataSource.INTERNAL:
            folder = get_url(folder, "NORA3/atmosphere/atm_hourly")
        if filename is None:
            filename = "arome3km_1hr_%Y%m.nc"
        return folder, filename

    def __call__(
        self,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        folder: str,
        expansion_factor: float = 1.2,
        program: str = "pyfimex",
        **kwargs,
    ):
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        start_times, end_times = create_monthly_stamps(start_time, end_time)

        msg.info(f"Getting wind forcing from NORA3 from {start_time} to {end_time}")

        temp_folder = "dnora_wnd_temp"
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        msg.plain("Removing old files from temporary folder...")
        for f in glob.glob("dnora_wnd_temp/*MetNo_NORA3.nc"):
            os.remove(f)

        # Define area to search in
        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

        # Setting resolution to roughly 3 km
        dlat = 3 / 111
        mean_lon_in_km = (lon_in_km(grid.lat()[0]) + lon_in_km(grid.lat()[-1])) * 0.5
        dlon = 3 / mean_lon_in_km

        wnd_list = []
        print("Apply >>> " + program)
        for n, (t0, t1) in enumerate(zip(start_times, end_times)):
            folder, filename = self._folder_filename(source, folder, filename=None)
            url = get_url(folder, filename, t0)
            msg.from_file(url)
            msg.plain(
                f"Reading wind: {t0.strftime('%Y-%m-%d %H:%M:00')}-{t1.strftime('%Y-%m-%d %H:%M:00')}"
            )
            nc_fimex = f"dnora_wnd_temp/wind_{n:04.0f}_MetNo_NORA3.nc"
            # Apply pyfimex or fimex
            if program == "pyfimex":
                pyfimex(
                    input_file=url,
                    output_file=nc_fimex,
                    projString="+proj=latlong +ellps=sphere +a=6371000 +e=0",
                    xAxisValues=np.arange(lon[0], lon[1] + dlon, dlon),
                    yAxisValues=np.arange(lat[0], lat[1] + dlat, dlat),
                    selectVariables=["wind_speed", "wind_direction"],
                    reduceTime_start=t0.strftime("%Y-%m-%dT%H:%M:%S"),
                    reduceTime_end=t1.strftime("%Y-%m-%dT%H:%M:%S"),
                )
            elif program == "fimex":
                fimex_command = [
                    "fimex",
                    "--input.file=" + url,
                    "--interpolate.method=bilinear",
                    "--interpolate.projString=+proj=latlong +ellps=sphere +a=6371000 +e=0",
                    "--interpolate.xAxisValues="
                    + str(lon[0])
                    + ","
                    + str(lon[0] + dlon)
                    + ",...,"
                    + str(lon[1])
                    + "",
                    "--interpolate.yAxisValues="
                    + str(lat[0])
                    + ","
                    + str(lat[0] + dlat)
                    + ",...,"
                    + str(lat[1])
                    + "",
                    "--interpolate.xAxisUnit=degree",
                    "--interpolate.yAxisUnit=degree",
                    "--process.rotateVector.all",
                    "--extract.selectVariables=wind_speed",
                    "--extract.selectVariables=wind_direction",
                    "--extract.reduceTime.start=" + t0.strftime("%Y-%m-%dT%H:%M:%S"),
                    "--extract.reduceTime.end=" + t1.strftime("%Y-%m-%dT%H:%M:%S"),
                    "--process.rotateVector.direction=latlon",
                    "--output.file=" + nc_fimex,
                ]
                call(fimex_command)
            wnd_list.append(xr.open_dataset(nc_fimex).squeeze())

        wind_forcing = xr.concat(wnd_list, dim="time")
        # Go to u and v components
        u, v = u_v_from_speed_dir(wind_forcing.wind_speed, wind_forcing.wind_direction)

        data_dict = {"u": u.fillna(0).values, "v": v.fillna(0).values}
        coord_dict = {
            "time": wind_forcing.time.values,
            "lon": wind_forcing.x.values,
            "lat": wind_forcing.y.values,
        }
        meta_dict = wind_forcing.attrs

        metaparameter_dict = {}

        return coord_dict, data_dict, meta_dict, metaparameter_dict


class MyWave3km(DataReader):
    """Reads wind data from the MyWave 3km hindcast directly from MET Norways
    servers. You should probably use MetNo_NORA3 because:

    The wind data is from NORA3 (see the MetNo_NORA3 reader), is taken
    from the wave model output. This means that model land points have no data.
    """

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

        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)
        return

    def _folder_filename(
        self, source: DataSource, folder: str, filename: str
    ) -> tuple[str]:
        if source == DataSource.REMOTE:
            folder = "https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/%Y/%m"
        if filename is None:
            filename = "%Y%m%d_MyWam3km_hindcast.nc"
        return folder, filename

    def __call__(
        self,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        folder: str,
        expansion_factor: float = 1.2,
        **kwargs,
    ):
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(
            start_time,
            end_time,
            self.stride,
            self.hours_per_file,
            self.last_file,
            self.lead_time,
        )

        msg.info(
            f"Getting wind forcing from MEPS from {self.start_time} to {self.end_time}"
        )
        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        temp_folder = "dnora_wnd_temp"
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        msg.plain("Removing old files from temporary folder...")
        for f in glob.glob("dnora_wnd_temp/*MetNo_MyWave3km.nc"):
            os.remove(f)

        # Define area to search in
        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

        # Setting resolution to roughly 3 km
        dlat = 3 / 111
        mean_lon_in_km = (
            lon_in_km(grid.edges("lat")[0]) + lon_in_km(grid.edges("lat")[-1])
        ) * 0.5
        dlon = 3 / mean_lon_in_km

        wnd_list = []
        for n in range(len(file_times)):
            folder, filename = self._folder_filename(source, folder, filename)
            url = get_url(folder, filename, file_times[n])

            msg.from_file(url)
            msg.plain(f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")

            nc_fimex = f"dnora_wnd_temp/wind_{n:04.0f}_MetNo_MyWave3km.nc"

            fimex_command = [
                "fimex",
                "--input.file=" + url,
                "--interpolate.method=bilinear",
                "--interpolate.projString=+proj=latlong +ellps=sphere +a=6371000 +e=0",
                "--interpolate.xAxisValues="
                + str(lon[0])
                + ","
                + str(lon[0] + dlon)
                + ",...,"
                + str(lon[1])
                + "",
                "--interpolate.yAxisValues="
                + str(lat[0])
                + ","
                + str(lat[0] + dlat)
                + ",...,"
                + str(lat[1])
                + "",
                "--interpolate.xAxisUnit=degree",
                "--interpolate.yAxisUnit=degree",
                "--process.rotateVector.all",
                "--extract.selectVariables=ff",
                "--extract.selectVariables=dd",
                "--extract.reduceTime.start="
                + start_times[n].strftime("%Y-%m-%dT%H:%M:%S"),
                "--extract.reduceTime.end="
                + end_times[n].strftime("%Y-%m-%dT%H:%M:%S"),
                "--process.rotateVector.direction=latlon",
                "--output.file=" + nc_fimex,
            ]

            call(fimex_command)
            wnd_list.append(xr.open_dataset(nc_fimex).squeeze())

        wind_forcing = xr.concat(wnd_list, dim="time")

        # Go to u and v components
        u, v = u_v_from_speed_dir(wind_forcing.ff, wind_forcing.dd)  # factor 1000

        u = -1 * u.fillna(0)  # *-1 due to ocean convection in WAM!!!
        v = -1 * v.fillna(0)  # *-1 due to ocean convection in WAM!!!

        data_dict = {"u": u, "v": v}
        coord_dict = {
            "time": wind_forcing.time.values,
            "lon": wind_forcing.rlon.values,
            "lat": wind_forcing.rlat.values,
        }
        meta_dict = wind_forcing.attrs
        metaparameter_dict = {}
        return coord_dict, data_dict, meta_dict, metaparameter_dict


class MEPS(DataReader):
    """Reads wind data from MET Norways MEPS forecast.

    The data is from a 2.5 km AROME model.
    """

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

        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)
        return

    def _folder_filename(
        self, source: DataSource, folder: str, filename: str
    ) -> tuple[str]:
        if source == DataSource.REMOTE:
            folder = "https://thredds.met.no/thredds/dodsC/meps25epsarchive/%Y/%m/%d"
        if filename is None:
            filename = f"meps_{self._prefix}_2_5km_%Y%m%dT%HZ.nc"
        return folder, filename

    def __call__(
        self,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        folder: str,
        expansion_factor: float = 1.2,
        filename: str = None,
        **kwargs,
    ):
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(
            start_time,
            end_time,
            self.stride,
            self.hours_per_file,
            self.last_file,
            self.lead_time,
        )

        msg.info(
            f"Getting wind forcing from MEPS from {self.start_time} to {self.end_time}"
        )
        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        temp_folder = "dnora_wnd_temp"
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            msg.info("Creating folder %s..." % temp_folder)

        msg.plain("Removing old files from temporary folder...")
        for f in glob.glob("dnora_wnd_temp/*MetNo_MEPS.nc"):
            os.remove(f)

        # Check weather to use 'det' or 'subset' files
        self._prefix = "det"
        folder, filename = self._folder_filename(source, folder, filename)
        url = get_url(folder, filename, file_times[0])
        try:
            xr.open_dataset(url)
        except:
            self._prefix = "subset"

        # Define area to search in
        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

        x_str, y_str = create_fimex_xy_strings(lon, lat, resolution_in_km=2.5)
        wnd_list = []
        for n in range(len(file_times)):
            msg.plain(f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")

            nc_fimex = f"dnora_wnd_temp/wind_{n:04.0f}_MetNo_MEPS.nc"
            folder, filename = self._folder_filename(source, folder, filename)
            url = get_url(folder, filename, file_times[n])
            msg.from_file(url)
            fimex_command = [
                "fimex",
                "--input.file=" + url,
                "--interpolate.method=bilinear",
                "--interpolate.projString=+proj=latlong +ellps=sphere +a=6371000 +e=0",
                "--interpolate.xAxisValues=" + x_str + "",
                "--interpolate.yAxisValues=" + y_str + "",
                "--interpolate.xAxisUnit=degree",
                "--interpolate.yAxisUnit=degree",
                "--process.rotateVector.all",
                "--extract.selectVariables=x_wind_10m",
                "--extract.selectVariables=y_wind_10m",
                "--extract.selectVariables=latitude",
                "--extract.selectVariables=longitude",
                "--extract.reduceTime.start="
                + start_times[n].strftime("%Y-%m-%dT%H:%M:%S"),
                "--extract.reduceTime.end="
                + end_times[n].strftime("%Y-%m-%dT%H:%M:%S"),
                "--process.rotateVector.direction=latlon",
                "--output.file=" + nc_fimex,
            ]

            if self._prefix == "subset":
                fimex_command.insert(
                    -2, "--extract.reduceDimension.name=ensemble_member"
                )
                fimex_command.insert(-2, "--extract.reduceDimension.start=1")
                fimex_command.insert(-2, "--extract.reduceDimension.end=1")

            call(fimex_command)

            wnd_list.append(xr.open_dataset(nc_fimex).squeeze())

        wind_forcing = xr.concat(wnd_list, dim="time")

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

        metaparameter_dict = {}

        return coord_dict, data_dict, meta_dict, metaparameter_dict


class NORA3_fp(DataReader):
    """Reads wind data of the NORA3 hindcast directly from MET Norways servers.

    The NORA3 HARMONIE-AROME high-resolution (ca 3 km) hindcast for the
    North Sea, the Norwegian Sea, and the Barents Sea.

    Haakenstad, H., Breivik, Ø., Furevik, B. R., Reistad, M., Bohlinger, P., &
    Aarnes, O. J. (2021). NORA3: A Nonhydrostatic High-Resolution Hindcast of
    the North Sea, the Norwegian Sea, and the Barents Sea,
    Journal of Applied Meteorology and Climatology, 60(10), 1443-1464,
    DOI: 10.1175/JAMC-D-21-0029.1
    """

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def __init__(
        self,
        stride: int = 1,
        hours_per_file: int = 1,
        last_file: str = "",
        lead_time: int = 4,
    ):
        """The data is currently in hourly files. Do not change the default
        setting unless you have a good reason to do so.
        """

        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)

    def __call__(
        self,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        folder: str,
        expansion_factor: float = 1.2,
        **kwargs,
    ):
        """Reads boundary spectra between given times and given area around
        the Grid object."""

        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(
            start_time,
            end_time,
            self.stride,
            self.hours_per_file,
            self.last_file,
            self.lead_time,
        )

        msg.info(
            f"Getting wind forcing from NORA3 from {self.start_time} to {self.end_time}"
        )

        temp_folder = "dnora_wnd_temp"
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        msg.plain("Removing old files from temporary folder...")
        for f in glob.glob("dnora_wnd_temp/*MetNo_NORA3.nc"):
            os.remove(f)

        # Define area to search in
        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

        # Set resolution to about 3 km
        dlat = 3 / 111
        mean_lon_in_km = (
            lon_in_km(grid.edges("lat")[0]) + lon_in_km(grid.edges("lat")[1])
        ) * 0.5
        dlon = 3 / mean_lon_in_km

        wnd_list = []
        for n in range(len(file_times)):
            url = self.get_url(
                file_times[n], start_times[n], first_ind=self.lead_time, source=source
            )
            msg.from_file(url)
            msg.plain(f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")

            nc_fimex = f"dnora_wnd_temp/wind_{n:04.0f}_MetNo_NORA3.nc"

            fimex_command = [
                "fimex",
                "--input.file=" + url,
                "--interpolate.method=bilinear",
                "--interpolate.projString=+proj=latlong +ellps=sphere +a=6371000 +e=0",
                "--interpolate.xAxisValues="
                + str(lon[0])
                + ","
                + str(lon[0] + dlon)
                + ",...,"
                + str(lon[1])
                + "",
                "--interpolate.yAxisValues="
                + str(lat[0])
                + ","
                + str(lat[0] + dlat)
                + ",...,"
                + str(lat[1])
                + "",
                "--interpolate.xAxisUnit=degree",
                "--interpolate.yAxisUnit=degree",
                "--process.rotateVector.all",
                "--extract.selectVariables=wind_speed",
                "--extract.selectVariables=wind_direction",
                "--extract.reduceTime.start="
                + start_times[n].strftime("%Y-%m-%dT%H:%M:%S"),
                "--extract.reduceTime.end="
                + end_times[n].strftime("%Y-%m-%dT%H:%M:%S"),
                "--process.rotateVector.direction=latlon",
                "--output.file=" + nc_fimex,
            ]
            # read_success = False
            # for ct in range(5):  # try 6 times
            #     try:
            #         call(fimex_command)
            #         read_success = True
            #     except:
            #         print(f'......Retry {ct}.....')
            #         time.sleep(10)  # wait for 10 seconds before re-trying
            #
            # # Don't want to catch the last execption
            # if not read_success:
            call(fimex_command)
            wnd_list.append(xr.open_dataset(nc_fimex).squeeze())

        wind_forcing = xr.concat(wnd_list, dim="time")

        # Fimex has already rotated the longitudes and latitudes, so calling them rlon/rlat is now incorrect
        #        wind_forcing = wind_forcing.rename_dims({'y': 'lat', 'x': 'lon'})
        # wind_forcing = wind_forcing.rename_vars({'y': 'lat', 'x': 'lon'})

        # Go to u and v components
        u, v = u_v_from_speed_dir(wind_forcing.wind_speed, wind_forcing.wind_direction)

        data_dict = {"u": u.fillna(0).values, "v": v.fillna(0).values}
        coord_dict = {
            "time": wind_forcing.time.values,
            "lon": wind_forcing.x.values,
            "lat": wind_forcing.y.values,
        }
        meta_dict = wind_forcing.attrs

        metaparameter_dict = {}

        return coord_dict, data_dict, meta_dict, metaparameter_dict

    def get_url(
        self, time_stamp_file, time_stamp, first_ind, source: DataSource, folder: str
    ) -> str:
        h0 = int(time_stamp_file.hour) % 6
        subfolder = time_stamp_file.strftime("%Y/%m/%d/") + (
            time_stamp_file - np.timedelta64(h0, "h")
        ).strftime("%H")
        ind = int((time_stamp.hour - first_ind) % 6) + first_ind
        filename = (
            +time_stamp_file.strftime("fc%Y%m%d")
            + (time_stamp_file - np.timedelta64(h0, "h")).strftime("%H")
            + f"_{ind:03d}_fp.nc"
        )
        if source == DataSource.REMOTE:
            folder = get_url("https://thredds.met.no/thredds/dodsC/nora3/", subfolder)
            return get_url(folder, filename)

        if source == DataSource.INTERNAL:
            folder = get_url(folder, "aaaaa")
            folder = get_url(folder, subfolder)
            return get_url(folder, filename)
        else:
            folder = get_url(folder, subfolder)
            return get_url(folder, filename)
