import xarray as xr
from copy import copy
from typing import Tuple
import pandas as pd
import numpy as np

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


class WAM4km(SpectralDataReader):
    def __init__(
        self,
        stride: int = 6,
        hours_per_file: int = 73,
        last_file: str = "",
        lead_time: int = 0,
    ) -> None:
        self.file_structure = FileStructure(
            stride=stride, hours_per_file=hours_per_file, last_file=last_file
        )

    def convention(self) -> SpectralConvention:
        return SpectralConvention.OCEAN

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def post_processing(self):
        return RemoveEmpty()

    def _folder_filename(
        self, source: DataSource, folder: str, filename: str
    ) -> tuple[str]:
        if source == DataSource.REMOTE:
            folder = (
                "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/%Y/%m/%d"
            )
        if filename is None:
            filename = f"MyWave_wam4_SPC_%Y%m%dT%HZ.nc"
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
        inds,
        filename: str = None,
        **kwargs,
    ) -> tuple[dict]:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        folder, filename = self._folder_filename(source, folder, filename)

        msg.info(f"Getting boundary spectra from WAM4 from {start_time} to {end_time}")

        data_vars = list(VAR_MAPPING.values())

        bnd_list = read_ds_list(
            start_time, end_time, inds, data_vars, folder, filename, self.file_structure
        )

        msg.info("Merging dataset together (this might take a while)...")
        bnd = xr.concat(bnd_list, dim="time").squeeze("y")

        coord_dict, data_dict, meta_dict = create_dicts(bnd, VAR_MAPPING)

        return coord_dict, data_dict, meta_dict


class NORA3(SpectralDataReader):
    def __init__(
        self,
        stride: int = 24,
        hours_per_file: int = 24,
        last_file: str = "",
        lead_time: int = 0,
    ) -> None:
        self.file_structure = FileStructure(
            stride=stride, hours_per_file=hours_per_file, last_file=last_file
        )

    def convention(self) -> SpectralConvention:
        return SpectralConvention.OCEAN

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def _folder_filename(
        self, source: DataSource, folder: str, filename: str
    ) -> tuple[str]:
        if source == DataSource.REMOTE:
            folder = "https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/%Y/%m"
        elif source == DataSource.INTERNAL:
            folder = get_url(folder, "WINDSURFER/mw3hindcast/spectra/%Y/%m")
        if filename is None:
            filename = "SPC%Y%m%d00.nc"
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
        inds,
        filename: str = None,
        **kwargs,
    ) -> tuple[dict]:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        folder, filename = self._folder_filename(source, folder, filename)

        start_times, end_times, file_times = self.file_structure.create_time_stamps(
            start_time, end_time
        )

        msg.info(f"Getting boundary spectra from NORA3 from {start_time} to {end_time}")

        data_vars = list(VAR_MAPPING.values())

        bnd_list = read_ds_list(
            start_time, end_time, inds, data_vars, folder, filename, self.file_structure
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
    def __init__(
        self,
        stride: int = 6,
        hours_per_file: int = 73,
        last_file: str = "",
        lead_time: int = 0,
        tile="POI",
    ) -> None:
        self.file_structure = FileStructure(
            stride=stride, hours_per_file=hours_per_file, last_file=last_file
        )
        self.tile = tile

    def convention(self) -> SpectralConvention:
        return SpectralConvention.OCEAN

    def post_processing(self):
        return RemoveEmpty()

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def _folder_filename(
        self, source: DataSource, folder: str, filename: str
    ) -> tuple[str]:
        if source == DataSource.REMOTE:
            folder = (
                "https://thredds.met.no/thredds/dodsC/ww3_4km_archive_files/%Y/%m/%d"
            )
        elif source == DataSource.IMMUTABLE:
            folder = get_url(folder, "DNMI_WAVE/%Y/%m/%d")
        if filename is None:
            filename = f"ww3_4km_{self.tile}_SPC_%Y%m%dT%HZ.nc"
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
        inds,
        filename: str = None,
        **kwargs,
    ) -> tuple[dict]:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        folder, filename = self._folder_filename(source, folder, filename)

        msg.info(
            f"Getting boundary spectra from WW3-4km from {start_time} to {end_time}"
        )

        data_vars = list(VAR_MAPPING.values())

        bnd_list = read_ds_list(
            start_time, end_time, inds, data_vars, folder, filename, self.file_structure
        )

        msg.info("Merging dataset together (this might take a while)...")
        bnd = xr.concat(bnd_list, dim="time").squeeze("y")

        coord_dict, data_dict, meta_dict = create_dicts(bnd, VAR_MAPPING)

        return coord_dict, data_dict, meta_dict


class WAM3(SpectralDataReader):
    """covers Nordic Seas and the Arctic"""

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
            offset=offset,
        )

    def convention(self) -> SpectralConvention:
        return SpectralConvention.OCEAN

    def default_data_source(self) -> DataSource:
        return DataSource.IMMUTABLE

    def _folder_filename(
        self, source: DataSource, folder: str, filename: str
    ) -> tuple[str]:
        if source == DataSource.REMOTE:
            folder = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam3_latest/"
            if filename is None:
                filename = "MyWave_wam3_WAVE_%Y%m%dT%HZ.nc"
        elif source == DataSource.IMMUTABLE:
            folder = get_url(folder, "DNMI_WAVE/%Y/%m/%d")
            if filename is None:
                filename = "MyWave_wam3_SPC_%Y%m%dT%HZ.nc"
        if filename is None:
            filename = "MyWave_wam3_SPC_%Y%m%dT%HZ.nc"
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
        inds,
        filename: str = None,
        **kwargs,
    ) -> tuple[dict]:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        folder, filename = self._folder_filename(source, folder, filename)

        msg.info(f"Getting boundary spectra from WAM3 from {start_time} to {end_time}")

        data_vars = list(VAR_MAPPING.values())

        bnd_list = read_ds_list(
            start_time, end_time, inds, data_vars, folder, filename, self.file_structure
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
            offset=offset,
        )
        self.tile = tile

    def convention(self) -> SpectralConvention:
        return SpectralConvention.OCEAN

    def default_data_source(self) -> DataSource:
        return DataSource.IMMUTABLE

    def _folder_filename(
        self, source: DataSource, folder: str, filename: str
    ) -> tuple[str]:
        if source == DataSource.REMOTE:
            folder = f"https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800{self.tile_names[self.tile][0].lower()}"
            if filename is None:
                filename = f"MyWave_wam800_{self.tile}SPC%H.nc"
        elif source == DataSource.IMMUTABLE:
            folder = get_url(folder, "DNMI_WAVE/%Y/%m/%d")
            if filename is None:
                filename = f"MyWave_wam800_{self.tile}SPC_%Y%m%dT%HZ.nc"
        if filename is None:
            filename = f"MyWave_wam800_{self.tile}SPC_%Y%m%dT%HZ.nc"
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
        inds,
        filename: str = None,
        **kwargs,
    ) -> tuple[dict]:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        folder, filename = self._folder_filename(source, folder, filename)

        msg.info(
            f"Getting boundary spectra from WAM800 {self.tile} from {start_time} to {end_time}"
        )
        data_vars = list(VAR_MAPPING.values())

        bnd_list = read_ds_list(
            start_time, end_time, inds, data_vars, folder, filename, self.file_structure
        )

        msg.info("Merging dataset together (this might take a while)...")
        bnd = xr.concat(bnd_list, dim="time").squeeze("y")

        coord_dict, data_dict, meta_dict = create_dicts(bnd, VAR_MAPPING)
        return coord_dict, data_dict, meta_dict
