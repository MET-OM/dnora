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

from dataclasses import dataclass

# mapping dnora variable name to the variable name found in MET Norway wave model netcdf's
VAR_MAPPING = {
    "spec": "SPEC",
    "lon": "longitude",
    "lat": "latitude",
    "time": "time",
    "dirs": "direction",
    "freq": "freq",
}


@dataclass
class FileStructure:
    stride: int
    hours_per_file: int
    lead_time: int = 0
    last_file: str = ""
    offset: int = 0

    def create_time_stamps(self, start_time: str, end_time: str):
        start_times, end_times, file_times = create_time_stamps(
            start_time,
            end_time,
            stride=self.stride,
            hours_per_file=self.hours_per_file,
            last_file=self.last_file,
            lead_time=self.lead_time,
            offset=self.offset,
        )
        return start_times, end_times, file_times


def get_first_ds(
    file_structure: FileStructure, start_time: str, folder: str, filename: str
) -> xr.Dataset:
    start_times, end_times, file_times = file_structure.create_time_stamps(
        start_time, start_time
    )
    url = get_url(folder, filename, file_times[0])
    ds = xr.open_dataset(url).isel(time=[0])

    return ds


def read_one_ds(
    start_time: str,
    end_time: str,
    file_times: list[str],
    urls: list[str],
    hours_per_file: int,
    n: int,
    inds: np.ndarray,
    data_vars: list[str],
    expected_lon: np.ndarray,
    expected_lat: np.ndarray,
) -> tuple[xr.Dataset, str]:
    """This functions reads one Dataset and crops it.
    If the expected file is not found, it goes to the previous one nad reads data with a longer lead time.

    This function can be used for e.g. forecast data where we have overlapping data in time in different files
    """

    file_time = file_times[n]
    ct = 1
    keep_trying = True
    try_next_file = False
    while keep_trying:

        try:
            url = urls[n]
            with xr.open_dataset(url) as f:
                ds = f.sel(
                    time=slice(start_time, end_time),
                    x=inds + 1,
                )[data_vars].copy()
                if file_is_consistent(ds, expected_lon, expected_lat, url):
                    try_next_file = False
                    keep_trying = False
                else:
                    try_next_file = True
        except OSError:
            try_next_file = True

        if try_next_file:
            if data_left_to_try_with(hours_per_file, n, ct, file_times, end_times):
                file_time = file_times[n - ct]
                ct += 1
            else:
                ds = None
                keep_trying = False

    return ds, url


def read_ds_list(
    start_time: str,
    end_time: str,
    inds,
    data_vars,
    folder: str,
    filename: str,
    file_structure: FileStructure,
) -> list[xr.Dataset]:
    start_times, end_times, file_times = file_structure.create_time_stamps(
        start_time, end_time
    )
    urls = [get_url(folder, filename, file_time) for file_time in file_times]
    ds_list = []
    expected_lon, expected_lat = None, None
    for n in range(len(file_times)):
        msg.plain(f"Reading data for: {start_times[n]}-{end_times[n]}")

        ds, url = read_one_ds(
            start_times[n],
            end_times[n],
            file_times,
            urls,
            file_structure.hours_per_file,
            n,
            inds,
            data_vars,
            expected_lon,
            expected_lat,
        )

        if ds is not None:
            msg.from_file(url)
            if not ds_list:
                expected_lon = ds.longitude.values
                expected_lat = ds.latitude.values
            ds_list.append(ds)

    return ds_list


def create_dicts(ds):
    coord_dict = {
        "time": ds[VAR_MAPPING["time"]].values.squeeze(),
        "lon": ds[VAR_MAPPING["lon"]].values.squeeze(),
        "lat": ds[VAR_MAPPING["lat"]].values.squeeze(),
        "freq": ds[VAR_MAPPING["freq"]].values.squeeze(),
        "dirs": ds[VAR_MAPPING["dirs"]].values.squeeze(),
    }
    data_dict = {"spec": ds[VAR_MAPPING["spec"]].data}

    meta_dict = ds.attrs
    return coord_dict, data_dict, meta_dict


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
        ds = get_first_ds(self.file_structure, start_time, folder, filename)

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

        coord_dict, data_dict, meta_dict = create_dicts(bnd)

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
        ds = get_first_ds(self.file_structure, start_time, folder, filename)

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
        coord_dict, data_dict, meta_dict = create_dicts(bnd)
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
        ds = get_first_ds(self.file_structure, start_time, folder, filename)

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

        coord_dict, data_dict, meta_dict = create_dicts(bnd)

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
        ds = get_first_ds(self.file_structure, start_time, folder, filename)

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

        coord_dict, data_dict, meta_dict = create_dicts(bnd)
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
        ds = get_first_ds(self.file_structure, start_time, folder, filename)

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

        coord_dict, data_dict, meta_dict = create_dicts(bnd)
        return coord_dict, data_dict, meta_dict


def data_left_to_try_with(hours_per_file, n, ct, file_times, end_times) -> bool:
    """Checks if we can go back one file and still have data covering the entire period.

    E.g. Each files conains 72 hours but we want to read in 6 hour chunks.
    If one ifle is missing we can use hours 7-12 in the previous file etc."""
    if n - ct < 0:
        return False

    if pd.Timestamp(end_times[n]) - pd.Timestamp(file_times[n - ct]) > pd.Timedelta(
        hours_per_file, "hours"
    ):
        return False

    return True


def file_is_consistent(
    ds: xr.Dataset, expected_lon: np.ndarray, expected_lat: np.ndarray, url: str
) -> bool:
    """Checks if dataset is consistent with the expected points"""
    if ds is None:
        msg.plain(f"SKIPPING, file not found: {url}")
        return False

    if (
        expected_lon is None  # always trust the first file that is read
        or (  # for the rest, check for consistency with first file
            (ds[VAR_MAPPING["lon"]] == expected_lon).all()
            and (ds[VAR_MAPPING["lat"]] == expected_lat).all()
        )
    ):
        return True

    msg.plain(f"SKIPPING, file inconsistent: {url}")
    return False
