from abc import ABC, abstractmethod
from copy import copy
import numpy as np
import xarray as xr
from subprocess import call
import os, glob
import time

# Import objects
from ..grd.grd_mod import Grid

# Import abstract classes
from .read import IceReader

# Import aux_funcsiliry functions
from .. import msg
from ..aux_funcs import create_time_stamps, u_v_from_speed_dir, expand_area, lon_in_km

from ..data_sources import DataSource


class NORA3(IceReader):
    """Reads ice data of the NORA3 wave hindcast directly from MET Norways servers.

    The NORA3 high-resolution (ca 3 km) hindcast for the Nordic Seas and the Arctic Ocean.

    Breivik, Ø., Carrasco, A., Haakenstad, H., Aarnes, O. J., Behrens, A., Bidlot, J.-R., et al. (2022). 
    The impact of a reduced high-wind Charnock parameter on wave growth with application to the North Sea, 
    the Norwegian Sea, and the Arctic Ocean. 
    Journal of Geophysical Research: Oceans, 127, e2021JC018196. https://doi.org/10.1029/2021JC018196
    """

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
    
    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE
    
    def __call__(
        self,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
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
            f"Getting ice forcing from NORA3 from {self.start_time} to {self.end_time}"
        )
        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        temp_folder = "dnora_ice_temp"
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        msg.plain("Removing old files from temporary folder...")
        for f in glob.glob("dnora_ice_temp/*MetNo_MyWave3km.nc"):
            os.remove(f)

        # Define area to search in
        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

        # Setting resolution to roughly 3 km
        dlat = 3 / 111
        mean_lon_in_km = (
            lon_in_km(grid.edges("lat")[0]) + lon_in_km(grid.edges("lat")[-1])
        ) * 0.5
        dlon = 3 / mean_lon_in_km

        ice_list = []
        for n in range(len(file_times)):
            url = self.get_url(file_times[n], source=source)

            msg.from_file(url)
            msg.plain(f"Reading ice forcing data: {start_times[n]}-{end_times[n]}")

            nc_fimex = f"dnora_ice_temp/ice_{n:04.0f}_MetNo_MyWave3km.nc"

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
                "--extract.selectVariables=SIC",
                "--extract.reduceTime.start="
                + start_times[n].strftime("%Y-%m-%dT%H:%M:%S"),
                "--extract.reduceTime.end="
                + end_times[n].strftime("%Y-%m-%dT%H:%M:%S"),
                "--process.rotateVector.direction=latlon",
                "--output.file=" + nc_fimex,
            ]

            call(fimex_command)
            ice_list.append(xr.open_dataset(nc_fimex).squeeze())

        ice = xr.concat(ice_list, dim="time")

        # Go to u and v components
        sic = ice.SIC.values

        time = ice.time.values
        lon = ice.rlon.values
        lat = ice.rlat.values
        x = None
        y = None

        return time, sic, lon, lat, x, y, ice.attrs

    def get_url(self, time_stamp, source):
        filename = (
            time_stamp.strftime("%Y")
            + time_stamp.strftime("%m")
            + time_stamp.strftime("%d")
            + "_MyWam3km_hindcast.nc"
        )

        if source == DataSource.REMOTE:
            return (
                "https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/"
                + time_stamp.strftime("%Y")
                + "/"
                + time_stamp.strftime("%m")
                + "/"
                + filename
            )
        else:
            pass
            return source + "/" + filename


class Barents25(IceReader):
    """Reads sea ice data of the Barents 2.5 km operational ocean model directly from MET Norways servers.

    Barents-2.5km is an operational data-assimilative coupled ocean and sea ice ensemble prediction model for 
    the Barents Sea and Svalbard

    Röhrs, J., Gusdal, Y., Rikardsen, E. S. U., Durán Moro, M., Brændshøi, J., Kristensen, N. M., 
    Fritzner, S., Wang, K., Sperrevik, A. K., Idžanović, M., Lavergne, T., Debernard, J. B., and Christensen, K. H.:
    Barents-2.5km v2.0: an operational data-assimilative coupled ocean and sea ice ensemble prediction model 
    for the Barents Sea and Svalbard, Geosci. Model Dev., 16, 5401-5426, 
    https://doi.org/10.5194/gmd-16-5401-2023, 2023.
    """

    def __init__(
        self,
        stride: int = 6,
        hours_per_file: int = 67,
        last_file: str = "",
        lead_time: int = 0,
        program: str = "pyfimex",
    ):
        """The data is currently in 6 hourly files. Do not change the default
        setting unless you have a good reason to do so.
        """

        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)
        self.program = program
        return

    def __call__(
        self,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        expansion_factor: float,
    ):
        """Reads in all grid points between the given times and at for the given indeces"""
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
            f"Getting ice forcing from Barents2.5km from {self.start_time} to {self.end_time}"
        )

        temp_folder = "dnora_ice_temp"
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        msg.plain("Removing old files from temporary folder...")
        for f in glob.glob("dnora_ice_temp/*MetNo_Barents25.nc"):
            os.remove(f)

        # Define area to search in
        lon_min, lon_max, lat_min, lat_max = expand_area(
            min(grid.lon()),
            max(grid.lon()),
            min(grid.lat()),
            max(grid.lat()),
            expansion_factor,
        )

        # Setting resolution to roughly 0.8 km
        dlat = 2.5 / 111
        mean_lon_in_km = (lon_in_km(min(grid.lat())) + lon_in_km(max(grid.lat()))) * 0.5
        dlon = 2.5 / mean_lon_in_km

        ocr_list = []
        print("Apply >>> " + self.program)
        for n in range(len(file_times)):
            url = self.get_url(file_times[n])

            msg.from_file(url)
            msg.plain(
                f"Reading ice forcing data: {start_times[n]}-{end_times[n]}"
            )

            nc_fimex = f"dnora_ice_temp/ice_{n:04.0f}_MetNo_Barents25.nc"
            # Apply pyfimex or fimex
            if self.program == "pyfimex":
                pyfimex(
                    input_file=url,
                    output_file=nc_fimex,
                    projString="+proj=latlong +ellps=sphere +a=6371000 +e=0",
                    xAxisValues=np.arange(lon_min, lon_max + dlon, dlon),
                    yAxisValues=np.arange(lat_min, lat_max + dlat, dlat),
                    selectVariables=["ice_concentration", "ice_thickness"],
                    reduceTime_start=start_times[n].strftime("%Y-%m-%dT%H:%M:%S"),
                    reduceTime_end=end_times[n].strftime("%Y-%m-%dT%H:%M:%S"),
                )
            elif self.program == "fimex":
                fimex_command = [
                    "fimex",
                    "--input.file=" + url,
                    "--interpolate.method=bilinear",
                    "--interpolate.projString=+proj=latlong +ellps=sphere +a=6371000 +e=0",
                    "--interpolate.xAxisValues="
                    + str(lon_min)
                    + ","
                    + str(lon_min + dlon)
                    + ",...,"
                    + str(lon_max)
                    + "",
                    "--interpolate.yAxisValues="
                    + str(lat_min)
                    + ","
                    + str(lat_min + dlat)
                    + ",...,"
                    + str(lat_max)
                    + "",
                    "--interpolate.xAxisUnit=degree",
                    "--interpolate.yAxisUnit=degree",
                    "--process.rotateVector.all",
                    #"--extract.selectVariables=ice_thickness",
                    "--extract.selectVariables=ice_concentration",
                    "--extract.reduceTime.start="
                    + start_times[n].strftime("%Y-%m-%dT%H:%M:%S"),
                    "--extract.reduceTime.end="
                    + end_times[n].strftime("%Y-%m-%dT%H:%M:%S"),
                    "--process.rotateVector.direction=latlon",
                    #'--extract.reduceDimension.name=depth',
                    #'--extract.reduceDimension.start=0',
                    #'--extract.reduceDimension.end=0',
                    "--output.file=" + nc_fimex,
                ]
                call(fimex_command)
            ocr_list.append(xr.open_dataset(nc_fimex).squeeze())

        ds = xr.concat(ocr_list, dim="time")

        # ice_forcing = ice_forcing.transpose("time", "lat", "lon")
        x = None
        y = None
        return (
            ds.time.values,
            ds.ice_concentration.values,
            #ds.ice_thickness.values,
            ds.X.values,
            ds.Y.values,
            x,
            y,
            ds.attrs,
        )

    @staticmethod
    def get_url(time_stamp):
          #https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_surface/2023/04/23/T00Z/barents_sfc_20230423T00Zm00.nc
        filename = (
            "barents_sfc_" + time_stamp.strftime("%Y%m%d") + "T" 
            + time_stamp.strftime("%H") + "Zm" + time_stamp.strftime("%H")+".nc"
        )
        url = "https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_surface/"
        + time_stamp.strftime("%Y")
        + "/"
        + time_stamp.strftime("%m")
        + "/"
        + time_stamp.strftime("%d")
        + "/T"
        + time_stamp.strftime("%H")
        + "Z/"
        + filename
        
        return url
