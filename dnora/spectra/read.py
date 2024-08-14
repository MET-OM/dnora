from dnora.spectral_conventions import SpectralConvention
from dnora import msg
from dnora.aux_funcs import create_time_stamps
from dnora.dnora_type_manager.data_sources import DataSource
from dnora.readers.abstract_readers import SpectralDataReader
from dnora.aux_funcs import get_url
from .swan_ascii import decode_lonlat, read_swan_ascii_spec
import pandas as pd
import numpy as np
from dnora.readers.ds_read_functions import read_ds_list, read_first_ds, create_dicts
from dnora.readers.file_structure import FileStructure
from dnora.aux_funcs import create_monthly_stamps
from functools import partial
import xarray as xr


class SWAN_Ascii(SpectralDataReader):
    def convention(self) -> SpectralConvention:
        return SpectralConvention.MET

    def default_data_source(self) -> DataSource:
        return DataSource.LOCAL

    def _folder_filename(
        self, source: DataSource, folder: str, filename: str
    ) -> tuple[str]:
        if source == DataSource.INTERNAL:
            folder = get_url(folder, "SWAN")
        return folder, filename

    def get_coordinates(
        self,
        grid,
        start_time,
        source: DataSource,
        folder: str,
        filename: str = None,
        **kwargs,
    ) -> dict:
        """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""

        folder, filename = self._folder_filename(source, folder, filename)
        url = get_url(folder, filename)
        with open(url, "r") as file:
            line = file.readline()
            while line:
                line = file.readline()
                if "LONLAT" in line:
                    lon, lat, file = decode_lonlat(file)
                    break

        all_points = {"lon": lon, "lat": lat}
        return all_points

    def __call__(
        self,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        inds,
        filename: str = None,
        **kwargs,
    ) -> tuple[dict]:
        """Reads in all boundary spectra between the given times and at for the given indeces"""

        msg.info(
            f"Getting boundary spectra from SWAN ascii file from {start_time} to {end_time}"
        )
        folder, filename = self._folder_filename(source, folder, filename)
        url = get_url(folder, filename)
        time, lon, lat, spec, freq, dirs = read_swan_ascii_spec(
            url, start_time=start_time, end_time=end_time
        )

        coord_dict = {
            "lon": lon,
            "lat": lat,
            "time": time,
            "freq": freq,
            "dirs": dirs,
        }
        data_dict = {"spec": spec}
        meta_dict = {"source": "Spectral wave data from a SWAN run"}

        return coord_dict, data_dict, meta_dict


WAM_VAR_MAPPING = {
    "spec": "SPEC",
    "lon": "longitude",
    "lat": "latitude",
    "time": "time",
    "dirs": "direction",
    "freq": "freq",
}


def ds_wam_xarray_read(
    start_time: pd.Timestamp,
    end_time: pd.Timestamp,
    url: str,
    inds: np.ndarray,
    data_vars: list[str],
):
    with xr.open_dataset(url) as f:
        ds = f.sel(
            time=slice(start_time, end_time),
            x=inds + 1,
        )[data_vars]
    return ds


class WAM(SpectralDataReader):

    stride = "month"  # int (for hourly), or 'month'
    hours_per_file = None  # int (if not monthly files)
    offset = 0  # int

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

    def convention(self) -> SpectralConvention:
        return SpectralConvention.WW3

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

        all_points = {"lon": ds.longitude.values[0], "lat": ds.latitude.values[0]}
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
        **kwargs,
    ) -> tuple[dict]:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
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
            f"Getting boundary spectra from {self.name()} from {start_time} to {end_time}"
        )

        data_vars = list(WAM_VAR_MAPPING.values())
        ds_creator_function = partial(
            ds_wam_xarray_read, inds=inds, data_vars=data_vars
        )
        bnd_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
            hours_per_file=hours_per_file,
        )

        msg.info("Merging dataset together (this might take a while)...")
        bnd = xr.concat(bnd_list, dim="time").squeeze("y")

        if "time" in list(bnd.longitude.coords):
            bnd["longitude"] = bnd.longitude[0, :]
            bnd["latitude"] = bnd.latitude[0, :]
        coord_dict, data_dict, meta_dict = create_dicts(bnd, WAM_VAR_MAPPING)

        return coord_dict, data_dict, meta_dict


WW3_VAR_MAPPING = {
    "spec": "efth",
    "lon": "longitude",
    "lat": "latitude",
    "time": "time",
    "dirs": "direction",
    "freq": "frequency",
}


def ds_ww3_xarray_read(
    start_time: pd.Timestamp,
    end_time: pd.Timestamp,
    url: str,
    inds: np.ndarray,
    data_vars: list[str],
):
    with xr.open_dataset(url) as f:
        ds = f.sel(
            time=slice(start_time, end_time),
            station=inds + 1,
        )[data_vars]
    return ds


class WW3(SpectralDataReader):
    _default_filename = "ww3_spec.%Y%m.nc"

    stride = "month"  # int (for hourly), or 'month'
    hours_per_file = None  # int (if not monthly files)
    offset = 0  # int

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

    def convention(self) -> SpectralConvention:
        return SpectralConvention.WW3

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

        all_points = {"lon": ds.longitude.values[0], "lat": ds.latitude.values[0]}
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
        **kwargs,
    ) -> tuple[dict]:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        folder = self._folder(folder, source)
        filename = self._filename(filename, source)

        if self.file_structure is None:
            start_times, end_times = create_monthly_stamps(start_time, end_time)
            file_times = start_times
        else:
            start_times, end_times, file_times = self.file_structure.create_time_stamps(
                start_time, end_time
            )
        msg.info(
            f"Getting boundary spectra from {self.name()} from {start_time} to {end_time}"
        )

        data_vars = list(WW3_VAR_MAPPING.values())
        ds_creator_function = partial(
            ds_ww3_xarray_read, inds=inds, data_vars=data_vars
        )
        bnd_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
        )

        msg.info("Merging dataset together (this might take a while)...")
        bnd = xr.concat(bnd_list, dim="time")

        if "time" in list(bnd.longitude.coords):
            bnd["longitude"] = bnd.longitude[0, :]
            bnd["latitude"] = bnd.latitude[0, :]
        coord_dict, data_dict, meta_dict = create_dicts(bnd, WW3_VAR_MAPPING)

        return coord_dict, data_dict, meta_dict
