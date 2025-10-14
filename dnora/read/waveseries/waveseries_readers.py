import xarray as xr
import numpy as np
from copy import deepcopy
import pandas as pd
from abc import ABC
from dnora.grid import Grid
import geo_parameters as gp
from pathlib import Path
from geo_skeletons import GriddedSkeleton
from dnora import msg, utils

from geo_parameters.metaparameter import MetaParameter

from dnora.spectra1d import Spectra1D
from dnora.process.spectra import CutFrequency
from dnora.utils.time import create_monthly_stamps
from dnora.read.abstract_readers import PointDataReader

from dnora.type_manager.dnora_types import DnoraDataType
from dnora.type_manager.data_sources import DataSource
from dnora.wave_parameters.parameters import get_function
from dnora.read.ds_read_functions import read_ds_list, read_first_ds
from functools import partial
from dnora.read.data_var_decoding import read_data_vars, compile_data_vars
from dnora.read.file_structure import FileStructure
from typing import Union, Optional


def ds_xarray_read(start_time, end_time, url):
    ds = xr.open_dataset(url).sel(time=slice(start_time, end_time))
    return ds


class Spectra1DToWaveSeries(PointDataReader):
    """Integrates Spectra1D to wave series"""

    def __init__(self, Spectra1D: Spectra1D, freq: tuple = (0, 10_000)) -> None:
        self._Spectra1D = deepcopy(Spectra1D)
        if Spectra1D is not None:
            self._freq = None
            if Spectra1D.freq()[0] < freq[0] or Spectra1D.freq()[-1] > freq[-1]:
                try:
                    self._Spectra1D.process(CutFrequency(freq))
                    self._freq = freq
                except AttributeError:
                    pass

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
        obj_type,
        grid,
        start_time,
        end_time,
        source,
        folder,
        filename,
        inds,
        dnora_class,
        parameters: list[MetaParameter] = None,
        **kwargs,
    ) -> tuple:
        time = (
            self._Spectra1D.time(data_array=True)
            .sel(time=slice(start_time, end_time))
            .values
        )
        lon = self._Spectra1D.lon(strict=True)
        lat = self._Spectra1D.lat(strict=True)
        x = self._Spectra1D.x(strict=True)
        y = self._Spectra1D.y(strict=True)

        obj = dnora_class(time=time, lon=lon, lat=lat, x=x, y=y)

        parameters = parameters or obj.core.data_vars()
        for var in parameters:
            if not gp.is_gp(var):
                meta = obj.core.meta_parameter(var)
            else:
                meta = var

            if meta is not None:
                func = get_function(meta)
                if func is not None:
                    msg.plain(f"{meta.__name__}...")
                    obj.set(var, func(self._Spectra1D))
                else:
                    msg.plain(f"Skipping {meta.__name__}...")

        if self._freq is not None:
            obj.meta.set({"integration_range": f"{self._freq[0]}-{self._freq[-1]} Hz"})

        return obj.ds()

    def name(self) -> str:
        if self._Spectra1D is None:
            return "EmptyData"
        return self._Spectra1D.name


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
    gp.wave.TpSea.standard_name(): "ptp0",
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
    _force_names: str = "gp"  #'gp' or 'source'
    _data_vars = []

    def __init__(
        self,
        stride: Union[
            int, str, None
        ] = None,  # Integer is number of hours, 'month' for monthly files
        hours_per_file: Optional[int] = None,  # None for stride = 'month'
        last_file: str = "",
        lead_time: int = 0,
        offset: Optional[int] = None,
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
        inds: list[int],
        obj_data_vars: list[str],
        dnora_class=None,
        force_names: str = None,
        **kwargs,
    ) -> tuple[dict]:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        if force_names is None:
            force_names = self._force_names

        keep_gp_names = False
        keep_source_names = False
        if force_names == "gp":
            keep_gp_names = True
        elif force_names == "source":
            keep_source_names = True

        if data_vars is None:
            data_vars = self._data_vars

        # If no data variables have been provided, read the ones that might be prsent in the class
        if data_vars:
            dynamic = True
        else:
            dynamic = False

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

        breakpoint()
        cls = dnora_objects.get(obj_type)
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


SWAN_ALIASES = {
    gp.wave.Hs.standard_name(): "hs",
    gp.wave.Dirp.standard_name(): "thetap",
    gp.wave.Tp.standard_name(): "tp",
}


def ds_xarray_swan_read(start_time, end_time, url, inds):
    ds = xr.open_dataset(url).sel(time=slice(start_time, end_time))
    ds = ds.stack(inds=("latitude", "longitude")).isel(inds=inds)
    return ds


class SWANnc(PointDataReader):
    stride = "month"  # int (for hourly), or 'month'
    hours_per_file = None  # int (if not monthly files)
    offset = 0  # int
    _force_names: str = "gp"  #'gp' or 'source'
    _decode_cf = False
    _data_vars = [
        gp.wave.Hs("hs"),
        gp.wave.Tp("tp"),
        gp.wave.Dirp("thetap"),
    ]

    def __init__(
        self,
        stride: Union[
            int, str, None
        ] = None,  # Integer is number of hours, 'month' for monthly files
        hours_per_file: Optional[int] = None,  # None for stride = 'month'
        last_file: str = "",
        lead_time: int = 0,
        offset: Optional[int] = None,
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
        points = GriddedSkeleton(lon=ds.longitude.values, lat=ds.latitude.values)
        lon, lat = points.lonlat()
        all_points = {"lon": lon, "lat": lat}
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
        inds: list[int],
        obj_data_vars: list[str],
        data_vars: list[str] = None,
        force_names: str = None,
        decode_cf: bool = None,
        **kwargs,
    ) -> tuple[dict]:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        if force_names is None:
            force_names = self._force_names

        keep_gp_names = False
        keep_source_names = False
        if force_names == "gp":
            keep_gp_names = True
        elif force_names == "source":
            keep_source_names = True

        if decode_cf is None:
            decode_cf = self._decode_cf
        if data_vars is None:
            data_vars = self._data_vars

        # If no data variables have been provided, read the ones that might be prsent in the class
        if not data_vars:
            data_vars = obj_data_vars

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
        msg.process("Compiling list of parameters accounting for known WW3 names:")
        data_vars = compile_data_vars(data_vars, aliases=SWAN_ALIASES)
        msg.blank()

        ds_creator_function = partial(ds_xarray_swan_read, inds=inds)
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
            "lon": ds.longitude.values,
            "lat": ds.latitude.values,
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
