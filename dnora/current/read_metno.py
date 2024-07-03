from copy import copy
import numpy as np
import xarray as xr
from subprocess import call
import os, glob

# Import objects
from dnora.grid import Grid

# Import abstract classes
from dnora.readers.abstract_readers import DataReader
from dnora.dnora_type_manager.data_sources import DataSource

# Import aux_funcsiliry functions
from dnora import msg
from dnora.aux_funcs import (
    create_time_stamps,
    expand_area,
    lon_in_km,
    pyfimex,
    get_url,
)


class NorKyst800(DataReader):
    """Reads ocean_current data of the NorKyst800 archieve directly from MET Norways servers.

    NorKyst-800 (Norwegian Coast 800m) is a numerical, high-resolution, ocean modelling
    system covering the Norwegian Coast.

    Albretsen, J., Sperrevik, A.K., Staalstrøm, A., Sandvik, A.D., Vikebø, F., Asplin, L., 2011.
    NorKyst-800 Rapport nr. 1 : Brukermanual og tekniske beskrivelser. NorKyst-800 Report
    No. 1 : User Manual and technical descriptions.
    """

    def __init__(
        self,
        stride: int = 24,
        hours_per_file: int = 24,
        last_file: str = "",
        lead_time: int = 0,
        program: str = "pyfimex",
    ):
        """The data is currently in daily files. Do not change the default
        setting unless you have a good reason to do so.
        """

        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)
        self.program = program
        return

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def _folder_filename(
        self, source: DataSource, folder: str, filename: str
    ) -> tuple[str]:
        if source == DataSource.REMOTE:
            if self.start_time.year < 2018:
                folder = "https://thredds.met.no/thredds/dodsC/sea/norkyst800mv0_1h/"
            else:
                folder = "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h"
        if filename is None:
            filename = "NorKyst-800m_ZDEPTHS_his.an.%Y%m%d00.nc"
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
            f"Getting ocean_current forcing from Norkyst800 from {self.start_time} to {self.end_time}"
        )
        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        temp_folder = "dnora_ocr_temp"
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        msg.plain("Removing old files from temporary folder...")
        for f in glob.glob("dnora_ocr_temp/*MetNo_Norkyst800.nc"):
            os.remove(f)

        # Define area to search in
        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

        # Setting resolution to roughly 0.8 km
        dlat = 0.8 / 111
        mean_lon_in_km = (lon_in_km(min(grid.lat())) + lon_in_km(max(grid.lat()))) * 0.5
        dlon = 0.8 / mean_lon_in_km

        ocr_list = []
        print("Apply >>> " + self.program)
        for n in range(len(file_times)):
            folder, filename = self._folder_filename(source, folder, filename)
            url = get_url(folder, filename, file_times[n])

            msg.from_file(url)
            msg.plain(
                f"Reading ocean_current forcing data: {start_times[n]}-{end_times[n]}"
            )

            nc_fimex = f"dnora_ocr_temp/ocean_current_{n:04.0f}_MetNo_Norkyst800.nc"
            # Apply pyfimex or fimex
            if self.program == "pyfimex":
                pyfimex(
                    input_file=url,
                    output_file=nc_fimex,
                    projString="+proj=latlong +ellps=sphere +a=6371000 +e=0",
                    xAxisValues=np.arange(lon[0], lon[1] + dlon, dlon),
                    yAxisValues=np.arange(lat[0], lat[1] + dlat, dlat),
                    selectVariables=["u", "v"],
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
                    "--extract.selectVariables=u",
                    "--extract.selectVariables=v",
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
