import xarray as xr
from copy import copy
from typing import Tuple
import pandas as pd
import numpy as np
from functools import partial

# Import abstract classes and needed instances of them
from .process import RemoveEmpty
from dnora.spectral_conventions import SpectralConvention

# Import aux_funcsiliry functions
from dnora import msg
from dnora.aux_funcs import create_time_stamps
from dnora.dnora_type_manager.data_sources import DataSource
from dnora.readers.abstract_readers import SpectralDataReader
from dnora.aux_funcs import get_url


from dnora.readers.ds_read_functions import read_ds_list, read_first_ds, create_dicts
from dnora.readers.file_structure import FileStructure

# mapping dnora variable name to the variable name found in MET Norway wave model netcdf's
VAR_MAPPING = {
    "spec": "SPEC",
    "lon": "longitude",
    "lat": "latitude",
    "time": "time",
    "dirs": "direction",
    "freq": "freq",
}


def ds_xarray_read(start_time, end_time, url, inds, data_vars):
    with xr.open_dataset(url) as f:
        ds = f.sel(
            time=slice(start_time, end_time),
            x=inds + 1,
        )[data_vars]
    return ds


class WAM4km(SpectralDataReader):
    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/%Y/%m/%d",
    }
    _default_filename = "MyWave_wam4_SPC_%Y%m%dT%HZ.nc"

    def __init__(
        self,
        stride: int = 6,
        hours_per_file: int = 73,
        last_file: str = "",
        lead_time: int = 0,
    ) -> None:
        self.file_structure = FileStructure(
            stride=stride,
            hours_per_file=hours_per_file,
            last_file=last_file,
            lead_time=lead_time,
        )

    def convention(self) -> SpectralConvention:
        return SpectralConvention.OCEAN

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def post_processing(self):
        return RemoveEmpty()

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
        ds = read_first_ds(self.file_structure, start_time, folder, filename)

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
        start_times, end_times, file_times = self.file_structure.create_time_stamps(
            start_time, end_time
        )
        msg.info(f"Getting boundary spectra from WAM4 from {start_time} to {end_time}")

        data_vars = list(VAR_MAPPING.values())
        ds_creator_function = partial(ds_xarray_read, inds=inds, data_vars=data_vars)
        bnd_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
            self.file_structure.hours_per_file,
        )

        msg.info("Merging dataset together (this might take a while)...")
        bnd = xr.concat(bnd_list, dim="time").squeeze("y")

        coord_dict, data_dict, meta_dict = create_dicts(bnd, VAR_MAPPING)

        return coord_dict, data_dict, meta_dict


class NORA3(SpectralDataReader):
    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/%Y/%m",
        DataSource.INTERNAL: "WINDSURFER/mw3hindcast/spectra/%Y/%m",
    }
    _default_filename = "SPC%Y%m%d00.nc"

    def __init__(
        self,
        stride: int = 24,
        hours_per_file: int = 24,
        last_file: str = "",
        lead_time: int = 0,
    ) -> None:
        self.file_structure = FileStructure(
            stride=stride,
            hours_per_file=hours_per_file,
            last_file=last_file,
            lead_time=lead_time,
        )

    def convention(self) -> SpectralConvention:
        return SpectralConvention.OCEAN

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

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
        ds = read_first_ds(self.file_structure, start_time, folder, filename)

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

        start_times, end_times, file_times = self.file_structure.create_time_stamps(
            start_time, end_time
        )

        msg.info(f"Getting boundary spectra from NORA3 from {start_time} to {end_time}")

        data_vars = list(VAR_MAPPING.values())

        ds_creator_function = partial(ds_xarray_read, inds=inds, data_vars=data_vars)
        bnd_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
            self.file_structure.hours_per_file,
        )

        msg.info("Merging dataset together (this might take a while)...")
        bnd = xr.concat(bnd_list, dim="time").squeeze("y")

        # Longitude and latitude defined over time also
        bnd["longitude"] = bnd.longitude[0, :]
        bnd["latitude"] = bnd.latitude[0, :]
        coord_dict, data_dict, meta_dict = create_dicts(bnd, VAR_MAPPING)
        meta_dict.pop("direction_convention")

        return coord_dict, data_dict, meta_dict


class WW3_4km(SpectralDataReader):

    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/ww3_4km_archive_files/%Y/%m/%d",
        DataSource.IMMUTABLE: "DNMI_WAVE/%Y/%m/%d",
    }

    def __init__(
        self,
        stride: int = 6,
        hours_per_file: int = 73,
        last_file: str = "",
        lead_time: int = 0,
        tile="POI",
    ) -> None:
        self.file_structure = FileStructure(
            stride=stride,
            hours_per_file=hours_per_file,
            last_file=last_file,
            lead_time=lead_time,
        )
        self.tile = tile
        self._default_filename = f"ww3_4km_{self.tile}_SPC_%Y%m%dT%HZ.nc"

    def convention(self) -> SpectralConvention:
        return SpectralConvention.OCEAN

    def post_processing(self):
        return RemoveEmpty()

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

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
        ds = read_first_ds(self.file_structure, start_time, folder, filename)

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
        start_times, end_times, file_times = self.file_structure.create_time_stamps(
            start_time, end_time
        )
        msg.info(
            f"Getting boundary spectra from WW3-4km from {start_time} to {end_time}"
        )

        data_vars = list(VAR_MAPPING.values())
        ds_creator_function = partial(ds_xarray_read, inds=inds, data_vars=data_vars)
        bnd_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
            self.file_structure.hours_per_file,
        )

        msg.info("Merging dataset together (this might take a while)...")
        bnd = xr.concat(bnd_list, dim="time").squeeze("y")

        coord_dict, data_dict, meta_dict = create_dicts(bnd, VAR_MAPPING)

        return coord_dict, data_dict, meta_dict


class WAM3(SpectralDataReader):
    """covers Nordic Seas and the Arctic"""

    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam3_latest/",
        DataSource.IMMUTABLE: "DNMI_WAVE/%Y/%m/%d",
    }
    _default_filenames = {DataSource.REMOTE: "MyWave_wam3_WAVE_%Y%m%dT%HZ.nc"}
    _default_filename = "MyWave_wam3_SPC_%Y%m%dT%HZ.nc"

    def __init__(
        self,
        stride: int = 12,
        hours_per_file: int = 121,
        last_file: str = "",
        lead_time: int = 0,
        offset: int = 6,
    ) -> None:
        self.file_structure = FileStructure(
            stride=stride,
            hours_per_file=hours_per_file,
            last_file=last_file,
            lead_time=lead_time,
            offset=offset,
        )

    def convention(self) -> SpectralConvention:
        return SpectralConvention.OCEAN

    def default_data_source(self) -> DataSource:
        return DataSource.IMMUTABLE

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

        ds = read_first_ds(self.file_structure, start_time, folder, filename)

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
        start_times, end_times, file_times = self.file_structure.create_time_stamps(
            start_time, end_time
        )
        msg.info(f"Getting boundary spectra from WAM3 from {start_time} to {end_time}")

        data_vars = list(VAR_MAPPING.values())
        ds_creator_function = partial(ds_xarray_read, inds=inds, data_vars=data_vars)
        bnd_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
            self.file_structure.hours_per_file,
        )

        msg.info("Merging dataset together (this might take a while)...")
        bnd = xr.concat(bnd_list, dim="time").squeeze("y")

        coord_dict, data_dict, meta_dict = create_dicts(bnd, VAR_MAPPING)
        return coord_dict, data_dict, meta_dict


class WAM800(SpectralDataReader):
    """c0 covers Finnmark, c1 covers NordNorge, c2 covers MidtNorge, c3 covers Vestlandet and c4 covers Skagerrak (these are the names of the different domains as used below)."""

    tile_names = {
        "c0": "Finnmark",
        "c1": "NordNorge",
        "c2": "MidtNorge",
        "c3": "Vestlandet",
        "c4": "Skagerrak",
    }

    def __init__(
        self,
        stride: int = 12,
        hours_per_file: int = 73,
        last_file: str = "",
        lead_time: int = 0,
        offset: int = 0,
        tile="c3",
    ) -> None:
        self.file_structure = FileStructure(
            stride=stride,
            hours_per_file=hours_per_file,
            last_file=last_file,
            lead_time=lead_time,
            offset=offset,
        )
        self.tile = tile
        self._default_folders = {
            DataSource.REMOTE: f"https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800{self.tile_names[self.tile][0].lower()}"
        }
        self._default_filenames = {
            DataSource.REMOTE: f"MyWave_wam800_{self.tile}SPC%H.nc"
        }
        self._default_filename = f"MyWave_wam800_{self.tile}SPC_%Y%m%dT%HZ.nc"

    def convention(self) -> SpectralConvention:
        return SpectralConvention.OCEAN

    def default_data_source(self) -> DataSource:
        return DataSource.IMMUTABLE

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
        ds = read_first_ds(self.file_structure, start_time, folder, filename)

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
        start_times, end_times, file_times = self.file_structure.create_time_stamps(
            start_time, end_time
        )
        msg.info(
            f"Getting boundary spectra from WAM800 {self.tile} from {start_time} to {end_time}"
        )
        data_vars = list(VAR_MAPPING.values())
        ds_creator_function = partial(ds_xarray_read, inds=inds, data_vars=data_vars)
        bnd_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
            self.file_structure.hours_per_file,
        )

        msg.info("Merging dataset together (this might take a while)...")
        bnd = xr.concat(bnd_list, dim="time").squeeze("y")

        coord_dict, data_dict, meta_dict = create_dicts(bnd, VAR_MAPPING)
        return coord_dict, data_dict, meta_dict
