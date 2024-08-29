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
from dnora.read.abstract_readers import PointDataReader
from dnora.waveseries import wave_parameters
import inspect
from dnora.dnora_type_manager.dnora_types import DnoraDataType
from dnora.dnora_type_manager.data_sources import DataSource
from dnora.wave_parameters.parameters import get_function
from dnora.read.ds_read_functions import read_ds_list, read_first_ds
from functools import partial


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

        meta_dict = self._Spectra1D.meta.get()
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


def ds_unstruct_ww3_read(
    start_time: pd.Timestamp,
    end_time: pd.Timestamp,
    url: str,
    inds: np.ndarray,
):
    with xr.open_dataset(url) as f:
        ds = f.sel(
            time=slice(start_time, end_time),
            node=inds,
        )
    return ds


WW3_ALIASES = {
    gp.wave.Hs.standard_name(): "hs",
    gp.wave.Tm01.standard_name(): "t01",
    gp.wave.Tm02.standard_name(): "t02",
    gp.wave.Tm_10.standard_name(): "t0m1",
    gp.wave.Dirm.standard_name(): "dir",
    gp.wave.Spr.standard_name(): "spr",
    gp.wave.Dirp.standard_name(): "dp",
    gp.wave.Tp.standard_name(): "tp",
    gp.ocean.WaterDepth.standard_name(): "dpt",
    gp.wind.XWind.standard_name(): "uwnd",
    gp.wind.YWind.standard_name(): "vwnd",
    gp.wave.HsSea.standard_name(): "phs0",
    gp.wave.HsSwell.standard_name(): "phs1",
    # gp.wave.TpSea.standard_name(): "ptp0",  ## Doesn'r exist!
    gp.wave.TpSwell.standard_name(): "ptp1",
    gp.wave.DirmSea.standard_name(): "pdir0",
    gp.wave.DirmSwell.standard_name(): "pdir1",
    gp.wave.DirpSea.standard_name(): "pdp0",
    gp.wave.DirpSwell.standard_name(): "pdp1",
    gp.wave.Tm01Sea.standard_name(): "pt01c0",
    gp.wave.Tm01Swell.standard_name(): "pt01c1",
    gp.wave.Tm02Sea.standard_name(): "pt02c0",
    gp.wave.Tm02Swell.standard_name(): "pt02c1",
    gp.wave.Tm_10Sea.standard_name(): "ptm10c0",
    gp.wave.Tm_10Swell.standard_name(): "ptm10c1",
}


class WW3Unstruct(PointDataReader):
    stride = "month"  # int (for hourly), or 'month'
    hours_per_file = None  # int (if not monthly files)
    offset = 0  # int
    _keep_source_names = False
    _keep_gp_names = False
    _decode_cf = False
    _data_vars = None

    def __init__(
        self,
        stride: (
            int | str | None
        ) = None,  # Integer is number of hours, 'month' for monthly files
        hours_per_file: int | None = None,  # None for stride = 'month'
        last_file: str = "",
        lead_time: int = 0,
        offset: int | None = None,
    ) -> None:
        if stride is not None:
            self.stride = stride

        if hours_per_file is not None:
            self.hours_per_file = hours_per_file

        if offset is not None:
            self.offset = offset

        if self.hours_per_file is not None:
            self.file_structure = FileStructure(
                stride=self.stride,
                hours_per_file=self.hours_per_file,
                last_file=last_file,
                lead_time=lead_time,
                offset=self.offset,
            )
        else:
            # This assumes monthly files!
            self.file_structure = None

    def default_data_source(self) -> DataSource:
        return DataSource.LOCAL

    def get_coordinates(
        self,
        grid,
        start_time,
        source: DataSource,
        folder: str,
        filename: str,
        **kwargs,
    ) -> dict:
        """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""
        folder = self._folder(folder, source)
        filename = self._filename(filename, source)
        ds = read_first_ds(folder, filename, start_time, self.file_structure)

        all_points = {"lon": ds.longitude.values, "lat": ds.latitude.values}
        return all_points

    def __call__(
        self,
        obj_type,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        filename: str,
        inds,
        data_vars: list[str] = None,
        keep_source_names: bool = None,
        keep_gp_names: bool = None,
        decode_cf: bool = None,
        **kwargs,
    ) -> tuple[dict]:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        if keep_source_names is None:
            keep_source_names = self._keep_source_names
        if keep_gp_names is None:
            keep_gp_names = self._keep_gp_names
        if decode_cf is None:
            decode_cf = self._decode_cf
        if data_vars is None:
            data_vars = self._data_vars

        folder = self._folder(folder, source)
        filename = self._filename(filename, source)
        if self.file_structure is None:
            start_times, end_times = create_monthly_stamps(start_time, end_time)
            file_times = start_times
            hours_per_file = None
        else:
            start_times, end_times, file_times = self.file_structure.create_time_stamps(
                start_time, end_time
            )
            hours_per_file = self.file_structure.hours_per_file

        msg.info(
            f"Getting waveseries data from {self.name()} from {start_time} to {end_time}"
        )

        msg.blank()
        msg.info("Compiling list of parameters:")
        data_vars = compile_data_vars(data_vars, aliases=WW3_ALIASES)
        msg.blank()

        ds_creator_function = partial(ds_unstruct_ww3_read, inds=inds)
        ds_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
            hours_per_file=hours_per_file,
        )

        msg.info("Merging dataset together (this might take a while)...")
        ds = xr.concat(ds_list, dim="time")

        coord_dict = {
            "time": ds.time.values,
            "lon": ds.longitude.values[0],
            "lat": ds.latitude.values[0],
        }
        data_dict = {}

        msg.blank()
        msg.info("Reading parameters:")

        if data_vars:
            msg.plain("Data variables specified. Reading only those!")
            data_dict = read_data_vars(
                data_vars, ds, keep_gp_names, keep_source_names, decode_cf
            )
        else:
            msg.plain("No data variables specified. Decoding from dataset!")
            data_vars = list(set(list(ds.data_vars)) - {"longitude", "latitude", "tri"})
            data_dict = read_data_vars(
                data_vars, ds, keep_gp_names, keep_source_names, decode_cf
            )

        meta_dict = ds.attrs

        return coord_dict, data_dict, meta_dict


def compile_data_vars(data_vars, aliases: dict = WW3_ALIASES):
    """Compiles data_vars accounting for possible aliases etc."""
    new_vars = []
    for var in data_vars:
        name, param = gp.decode(var)
        if param is not None:
            name = param.standard_name()
        if aliases.get(name) is not None:
            if gp.is_gp_class(param):
                param = param(aliases.get(name))  # Create e.g. Dirp('dp')
            else:
                param.name = aliases.get(name)
            msg.plain(
                f"Identified that {type(param).__name__} is called {aliases.get(name)} in the source!"
            )
            new_vars.append(param)
        else:
            new_vars.append(name)
    return new_vars


def read_data_vars(
    data_vars: list[str, gp.metaparameter.MetaParameter],
    ds: xr.Dataset,
    keep_gp_names: bool,
    keep_source_names: bool,
    decode_cf: bool,
) -> dict:
    """Read the given variables from the xarray dataset

    Uninitilized geo-parameter, e.g. gp.wave.Hs
        - Decode using standard-name. If not found, try default name of parameter (e.g. 'hs').
    Uninitilized geo-parameter, e.g. gp.wave.Hs('hsig')
        - Decode using standard-name. If not found, try given name of parameter ('hsig').
    String, e.g. 'hs'
        - Try to read that variable from dataset

    keep_source_names = True
        - Keep the variable name that is used by the data source even if a geo-parameter is given
        - E.g. gp.wave.Hs can match to 'hsig' by standard_name and 'hsig' will be used ('hs' if false)

    keep_gp_names = True
        - Create the variable names using the default class name even if it is initialized
        - E.g. gp.wave.Hs('swh') matches 'swh' variable with no metadata, but the variable created is still called 'hs'

    decode_cf = True
        - Try to decode standard_name from source and create geo-parameters if only string given
        - E.g. given 'hs' will try to find standard name in metadata of 'hs' and create gp.wave.Hs
    """

    if keep_gp_names and keep_source_names:
        raise ValueError("keep_source_names and keep_gp_names cannot both be True!")
    data_dict = {}
    for var in data_vars:
        var, param = gp.decode(var, init=True)
        if param is not None:
            var_in_ds = param.find_me_in_ds(ds)
            if var_in_ds is None:
                var_in_ds = var
            if keep_gp_names:
                param.name = gp.get(param.standard_name()).name
            elif keep_source_names:
                param.name = var_in_ds

            if hasattr(ds[var_in_ds], "standard_name"):
                std_name = ds[var_in_ds].standard_name
            else:
                std_name = ""

            msg.plain(
                f"[{type(param).__name__}('{param.name}')] << '{var_in_ds}' {std_name}"
            )
            data_dict[param] = ds.get(var_in_ds).data
        elif ds.get(var) is not None:
            if hasattr(ds[var], "standard_name"):
                std_name = ds[var].standard_name
            else:
                std_name = ""
            if std_name and decode_cf:
                param = gp.get(std_name)
            else:
                param = None
            if param is not None:
                if keep_gp_names:
                    param = param()
                else:
                    param = param(var)
                msg.plain(
                    f"[{type(param).__name__}('{param.name}')] << '{var}' {std_name}"
                )
                data_dict[param] = ds.get(var).data
            else:
                msg.plain(f"'{var}' << {var}")
                data_dict[var] = ds.get(var).data
        else:
            msg.plain(f"'{var}' unknown. Skipping.")
    return data_dict
