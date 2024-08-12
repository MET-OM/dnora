import xarray as xr
import numpy as np
from copy import deepcopy
import pandas as pd
from abc import ABC, abstractmethod
from typing import Union
from dnora.grid import Grid
import geo_parameters as gp
from pathlib import Path

# Import aux_funcsiliry functions
from dnora import msg, aux_funcs

# from .conventions import Spectra1DlConvention
import geo_parameters as gp
from geo_parameters.metaparameter import MetaParameter

from dnora.spectra1d import Spectra1D, process
from dnora.aux_funcs import get_url, create_monthly_stamps
from dnora.readers.abstract_readers import PointDataReader
from dnora.waveseries import wave_parameters
import inspect
from dnora.dnora_type_manager.dnora_types import DnoraDataType
from dnora.dnora_type_manager.data_sources import DataSource
from dnora.wave_parameters.parameters import get_function
from dnora.readers.ds_read_functions import read_ds_list, read_first_ds


def ds_xarray_read(start_time, end_time, url):
    ds = xr.open_dataset(url).sel(time=slice(start_time, end_time))
    return ds


class Spectra1DToWaveSeries(PointDataReader):
    """Integrates Spectra1D to wave series"""

    def __init__(self, Spectra1D: Spectra1D, freq: tuple = (0, 10_000)) -> None:
        self._Spectra1D = deepcopy(Spectra1D)
        try:
            self._Spectra1D.process_Spectra1D(process.CutFrequency(freq))
        except AttributeError:
            msg.warning(
                f"Object {self.name()} does not have a process_Spectra1D method!\nNot cutting any frequencies!"
            )
        self._freq = freq

    def default_data_source(self) -> DataSource:
        return DataSource.CREATION

    def get_coordinates(
        self,
        grid,
        start_time,
        source: DataSource,
        folder: str,
        filename: str,
        **kwargs,
    ) -> dict:
        all_points = {
            "lon": self._Spectra1D.lon(strict=True),
            "lat": self._Spectra1D.lat(strict=True),
            "x": self._Spectra1D.x(strict=True),
            "y": self._Spectra1D.y(strict=True),
        }
        return all_points

    def __call__(
        self,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        filename: str,
        inds,
        parameters: list[MetaParameter] = [
            gp.wave.Hs,
            gp.wave.Tm01,
            gp.wave.Tm02,
            gp.wave.Tm_10,
            gp.wave.Dirm,
        ],
        **kwargs,
    ) -> tuple:
        self.name = self._Spectra1D.name
        time = (
            self._Spectra1D.time(data_array=True)
            .sel(time=slice(start_time, end_time))
            .values
        )
        lon = self._Spectra1D.lon(strict=True)
        lat = self._Spectra1D.lat(strict=True)
        x = self._Spectra1D.x(strict=True)
        y = self._Spectra1D.y(strict=True)

        data_dict = {}
        for wp in parameters:
            func = get_function(wp)
            data_dict[wp] = func(self._Spectra1D)

        meta_dict = self._Spectra1D.ds().attrs
        meta_dict["integration_range"] = f"{self._freq[0]}-{self._freq[-1]} Hz"

        coord_dict = {"lon": lon, "lat": lat, "x": x, "y": y, "time": time}

        return coord_dict, data_dict, meta_dict

    def name(self) -> str:
        if self._Spectra1D is None:
            return "EmptyData"
        return self._Spectra1D.name


class E39(PointDataReader):
    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/obs/buoy-svv-e39/%Y/%m"
    }

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def __init__(self, loc: str = "D", mode="wave"):
        self._loc = loc  # Given as "D", or "D_Breisundet"
        self._mode = mode  # 'wind' or 'wave'
        self._default_filename = f"%Y%m_E39_{self.loc()}_{self._mode}.nc"

    def _buoy_dict(self) -> dict:
        return {
            "A": "A_Sulafjorden",
            "B": "B_Sulafjorden",
            "B1": "B1_Sulafjorden",
            "C": "C_Sulafjorden",
            "C1": "C1_Sulafjorden",
            "D": "D_Breisundet",
            "F": "F_Vartdalsfjorden",
            "G": "G_Halsafjorden",
        }

    def loc(self) -> list[str]:
        if len(self._loc) > 2:
            return self._loc
        else:
            return self._buoy_dict()[self._loc]

    # def _folder_filename(
    #     self, source: DataSource, folder: str, filename: str
    # ) -> tuple[str]:
    #     if source == DataSource.REMOTE:
    #         folder = "https://thredds.met.no/thredds/dodsC/obs/buoy-svv-e39/%Y/%m"
    #     if filename is None:
    #         filename = f"%Y%m_E39_{self.loc()}_{self._mode}.nc"
    #     return folder, filename

    def get_coordinates(
        self,
        grid,
        start_time,
        source: DataSource,
        folder: str,
        filename: str,
        **kwargs,
    ) -> dict:
        # start_time = pd.to_datetime(start_time)
        folder = self._folder(folder, source)
        filename = self._filename(filename, source)
        ds = read_first_ds(folder, filename, start_time)
        return {"lon": ds.longitude.values, "lat": ds.latitude.values}

    def __call__(
        self,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        filename: str,
        inds,
        **kwargs,
    ) -> tuple:
        # loc = np.array(self._buoys())[inds][0]

        folder = self._folder(folder, source)
        filename = self._filename(filename, source)

        start_times, end_times = create_monthly_stamps(start_time, end_time)
        file_times = start_times

        ds_creator_function = ds_xarray_read
        ds_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
        )
        # breakpoint()
        # ds_list = []
        # for month in months:
        #     url = get_url(folder, filename, month)
        #     ds_list.append(xr.open_dataset(url))
        ds = xr.concat(ds_list, dim="time")
        ds["lon"] = np.nanmedian(ds.longitude.values)
        ds["lat"] = np.nanmedian(ds.latitude.values)

        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds)
        data_dict = {}
        for var in ds.data_vars:
            if var not in ["lon", "lat", "longitude", "latitude", "x", "y"]:
                if hasattr(ds[var], "standard_name"):
                    meta_param = gp.get(
                        ds[var].standard_name
                    )  # Find geo-parameter based on standard name
                else:
                    meta_param = None
                if meta_param is not None:
                    data_dict[meta_param] = np.swapaxes(
                        np.expand_dims(ds.get(var).values, axis=1), 0, 1
                    )

        coord_dict = {"time": ds.time.values, "lon": lon, "lat": lat}
        meta_dict = ds.attrs

        return coord_dict, data_dict, meta_dict


# class WW3Nc(WaveSeriesReader):
#     def __init__(
#         self, filename: str = "ww3.%Y%m.nc", folder: str = "", mode: str = "single"
#     ):
#         """Mode can be 'single' (one file), 'monthly'"""
#         self._filename = filename
#         self._mode = mode
#         self._folder = folder

#     def _filenames(self, start_time, end_time, folder):
#         filenames = []
#         if self._mode == "single":
#             filenames.append(f"{folder}/{self._filename}")
#         else:
#             for file in aux_funcs.month_list(start_time, end_time, fmt=self._filename):
#                 filenames.append(f"{folder}/{file}")
#         return filenames

#     def _convert_var(self, var: str) -> str:
#         var_dict = {
#             "hs": "hs",
#             "t02": "tm02",
#             "fp": "fp",
#             "t0m1": "tm_10",
#             "t01": "tm01",
#             "spr": "sprm",
#             "dir": "dirm",
#             "dp": "dirp",
#             "wnd": "ff",
#             "wnddir": "dd",
#         }
#         return var_dict.get(var)

#     def get_coordinates(
#         self, grid, start_time, source: DataSource, folder: str
#     ) -> tuple:
#         """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""
#         # day = pd.date_range(start_time, start_time,freq='D')
#         filename = self._filenames(start_time, start_time, folder=self._folder)[0]
#         aux_funcs.check_if_file(filename, halt=True)
#         ds = xr.open_dataset(filename).isel(time=[0])
#         lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds)

#         return lon, lat, x, y

#     def __call__(
#         self, grid, start_time, end_time, inds, source: DataSource, **kwargs
#     ) -> tuple:
#         """Reads in all wave data between the given times and at for the given indeces"""

#         msg.info(
#             f"Getting wave data from WW3 netcdf files from {start_time} to {end_time}"
#         )

#         # for file in self._filenames(start_time, end_time, self._folder):
#         msg.from_file(self._filenames(start_time, end_time, self._folder))

#         def _crop(ds):
#             """Crop to given time and index."""

#             try:  # Unstructured mesh output
#                 cropped_ds = ds.sel(
#                     time=slice(start_time, end_time), node=inds
#                 ).drop_dims(["noel", "element"])
#                 if not self._file_identification_message_printed:
#                     msg.info("Identified file as WW3 unstructured output")
#                     self._file_identification_message_printed = True
#                 return cropped_ds
#             except KeyError:
#                 try:  # WW3 Spectra1Dl output (get e.g. wind, hs from those files)
#                     cropped_ds = ds.sel(
#                         time=slice(start_time, end_time), station=(inds + 1)
#                     ).drop_dims(["frequency", "direction", "string40"])
#                     if not self._file_identification_message_printed:
#                         msg.info("Identified file as WW3 Spectra1Dl output")
#                         self._file_identification_message_printed = True
#                     return cropped_ds
#                 except:
#                     raise KeyError(
#                         "Could not identify the file as a known format (WW3 Unstructured or Spectra1Dl output)!"
#                     )

#         import dask

#         with dask.config.set(**{"array.slicing.split_large_chunks": True}):
#             self._file_identification_message_printed = False
#             with xr.open_mfdataset(
#                 self._filenames(start_time, end_time, self._folder), preprocess=_crop
#             ) as ds:
#                 lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds)

#                 data = {}
#                 for var in ds.data_vars:
#                     # if var not in ['lon', 'lat', 'longitude', 'latitude', 'x', 'y']:
#                     dnora_var = self._convert_var(var)
#                     if dnora_var is not None:
#                         # if dnora_var == 'fp':
#                         #    data[get_wave_parameter('tp')] = np.swapaxes(ds.get(var).values**-1,0,1)
#                         # else:
#                         msg.info(f"Importing {dnora_var} << {var}")
#                         data[get_wave_parameter(dnora_var)] = np.swapaxes(
#                             ds.get(var).values, 0, 1
#                         )
#                     else:
#                         msg.info(
#                             f"Could not identify {var} as a DNORA WaveSeries variable"
#                         )

#                 return ds.time.values, data, lon, lat, x, y, ds.attrs


# class WW3Nc_old(WaveSeriesReader):
#     def __init__(self, filename: str, parameters: list[str] = ["hs", "tp"]):
#         self._filename = filename

#     def _convert_var(self, var: str) -> str:
#         var_dict = {
#             "hs": "hs",
#             "t02": "tm02",
#             "fp": "fp",
#             "t0m1": "tm_10",
#             "t01": "tm01",
#             "spr": "sprm",
#             "dir": "dirm",
#             "dp": "dirp",
#         }
#         return var_dict.get(var)

#     def get_coordinates(
#         self, grid, start_time, source: DataSource, folder: str
#     ) -> tuple:
#         start_time = pd.to_datetime(start_time)
#         ds = xr.open_dataset(self._filename).isel(time=[0])
#         return ds.longitude.values, ds.latitude.values, None, None

#     def __call__(
#         self, grid, start_time, end_time, inds, source: DataSource, **kwargs
#     ) -> tuple:
#         ds = xr.open_dataset(self._filename).sel(node=inds)

#         ds["lon"] = ds.longitude.values
#         ds["lat"] = ds.latitude.values
#         lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds)

#         data = {}
#         for var in ds.data_vars:
#             if var not in ["lon", "lat", "longitude", "latitude", "x", "y"]:
#                 dnora_var = self._convert_var(var)
#                 if dnora_var is not None:
#                     if dnora_var == "fp":
#                         data[get_wave_parameter("tp")] = np.swapaxes(
#                             ds.get(var).values ** -1, 0, 1
#                         )
#                     else:
#                         data[get_wave_parameter(dnora_var)] = np.swapaxes(
#                             ds.get(var).values, 0, 1
#                         )
#         return ds.time.values, data, lon, lat, x, y, ds.attrs
