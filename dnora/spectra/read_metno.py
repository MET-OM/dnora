import xarray as xr
from copy import copy
from typing import Tuple
import pandas as pd

# Import abstract classes and needed instances of them
from .process import RemoveEmpty
from dnora.spectral_conventions import SpectralConvention

# Import aux_funcsiliry functions
from dnora import msg
from dnora.aux_funcs import create_time_stamps
from dnora.dnora_type_manager.data_sources import DataSource
from dnora.readers.abstract_readers import SpectralDataReader
from dnora.aux_funcs import get_url


class WAM4km(SpectralDataReader):
    def __init__(
        self,
        ignore_nan: bool = True,
        stride: int = 6,
        hours_per_file: int = 73,
        last_file: str = "",
        lead_time: int = 0,
    ) -> None:
        self.ignore_nan = ignore_nan
        self.stride = stride
        self.hours_per_file = hours_per_file
        self.lead_time = lead_time
        self.last_file = last_file

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
                "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/%Y/%m"
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

        start_times, end_times, file_times = create_time_stamps(
            start_time,
            start_time,
            stride=self.stride,
            hours_per_file=self.hours_per_file,
            last_file=self.last_file,
            lead_time=self.lead_time,
        )
        folder, filename = self._folder_filename(source, folder, filename)
        url = self.get_url(folder, filename, file_times[0])

        data = xr.open_dataset(url).isel(time=[0])

        all_points = {"lon": data.longitude.values[0], "lat": data.latitude.values[0]}
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
        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(
            start_time,
            end_time,
            stride=self.stride,
            hours_per_file=self.hours_per_file,
            last_file=self.last_file,
            lead_time=self.lead_time,
        )

        msg.info(
            f"Getting boundary spectra from WAM4 from {self.start_time} to {self.end_time}"
        )
        bnd_list = []
        for n in range(len(file_times)):
            msg.plain(f"Reading boundary spectra: {start_times[n]}-{end_times[n]}")

            file_time = file_times[n]
            ct = 1
            keep_trying = True
            while keep_trying:
                folder, filename = self._folder_filename(source, folder, filename)
                url = self.get_url(folder, filename, file_time)

                try:
                    with xr.open_dataset(url) as f:
                        this_ds = f.sel(
                            time=slice(start_times[n], end_times[n]),
                            x=inds + 1,
                        )[
                            [
                                "SPEC",
                                "longitude",
                                "latitude",
                                "time",
                                "freq",
                                "direction",
                            ]
                        ].copy()
                except OSError:
                    this_ds = None

                file_consistent = self.file_is_consistent(this_ds, bnd_list, url)

                if file_consistent:
                    keep_trying = False

                elif self.data_left_to_try_with(n, ct, file_times, end_times):
                    file_time = file_times[n - ct]
                    ct += 1
                else:
                    this_ds = None
                    keep_trying = False

            # We are now out of the while loop
            if this_ds is not None:
                msg.from_file(url)
                bnd_list.append(this_ds)
                this_ds = None

        msg.info("Merging dataset together (this might take a while)...")
        bnd = xr.concat(bnd_list, dim="time").squeeze("y")

        coord_dict = {
            "time": bnd.time.values,
            "lon": bnd.longitude.values,
            "lat": bnd.latitude.values,
            " freq": bnd.freq.values,
            "dirs": bnd.direction.values,
        }
        data_dict = {"spec": bnd.SPEC.values}

        meta_dict = bnd.attrs

        return coord_dict, data_dict, meta_dict

    def file_is_consistent(self, this_ds, bnd_list, url) -> bool:
        if this_ds is None:
            msg.plain(f"SKIPPING, file not found: {url}")
            return False

        if (
            not bnd_list  # always trust the first file that is read
            or (  # for the rest, check for consistency with first file
                (this_ds.longitude == bnd_list[0].longitude).all()
                and (this_ds.latitude == bnd_list[0].latitude).all()
            )
        ):
            return True

        msg.plain(f"SKIPPING, file inconsistent: {url}")
        return False

    def data_left_to_try_with(self, n, ct, file_times, end_times) -> bool:
        if n - ct < 0:
            return False

        if pd.Timestamp(end_times[n]) - pd.Timestamp(file_times[n - ct]) > pd.Timedelta(
            self.hours_per_file, "hours"
        ):
            return False

        return True


class NORA3(SpectralDataReader):
    def __init__(
        self,
        stride: int = 24,
        hours_per_file: int = 24,
        last_file: str = "",
        lead_time: int = 0,
    ) -> None:
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)

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
            folder = get_url(folder, "WINDSURFER/mw3hindcast/spectra")
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
        start_times, end_times, file_times = create_time_stamps(
            start_time,
            start_time,
            stride=self.stride,
            hours_per_file=self.hours_per_file,
            last_file=self.last_file,
            lead_time=self.lead_time,
        )

        folder, filename = self._folder_filename(source, folder, filename)
        url = get_url(folder, filename, file_times[0])
        data = xr.open_dataset(url).isel(time=[0])

        all_points = {"lon": data.longitude.values[0], "lat": data.latitude.values[0]}
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
        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(
            start_time,
            end_time,
            stride=self.stride,
            hours_per_file=self.hours_per_file,
            last_file=self.last_file,
            lead_time=self.lead_time,
        )

        msg.info(
            f"Getting boundary spectra from NORA3 from {self.start_time} to {self.end_time}"
        )
        bnd_list = []
        for n in range(len(file_times)):
            folder, filename = self._folder_filename(source, folder, filename)
            url = get_url(folder, filename, file_times[n])
            msg.from_file(url)
            msg.plain(f"Reading boundary spectra: {start_times[n]}-{end_times[n]}")
            with xr.open_dataset(url) as f:
                this_ds = f.sel(time=slice(start_times[n], end_times[n]), x=(inds + 1))[
                    ["SPEC", "longitude", "latitude", "time", "freq", "direction"]
                ].copy()
            bnd_list.append(this_ds)
        bnd = xr.concat(bnd_list, dim="time").squeeze("y")

        coord_dict = {
            "lon": bnd.longitude.values[0, :],
            "lat": bnd.latitude.values[0, :],
            "time": bnd.time.values,
            "freq": bnd.freq.values,
            "dirs": bnd.direction.values,
        }
        data_dict = {"spec": bnd.SPEC.values}
        meta_dict = bnd.attrs
        meta_dict.pop("direction_convention")
        return coord_dict, data_dict, meta_dict


class WW3_4km(SpectralDataReader):
    def __init__(
        self,
        ignore_nan: bool = True,
        stride: int = 6,
        hours_per_file: int = 73,
        last_file: str = "",
        lead_time: int = 0,
        cache: bool = True,
        clean_cache: bool = False,
        tile="POI",
    ) -> None:
        self.ignore_nan = copy(ignore_nan)
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)
        self.tile = (
            tile  # 'POI for all the domain', 'C0', ... 'C5' for WAM800 coastal areas
        )
        return

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
        if filename is None:
            filename = f"ww3_4km_POI_SPC_%Y%m%dT%HZ.nc"
        return folder, filename

    def get_coordinates(
        self,
        start_time,
        source: DataSource,
        folder: str,
        filename: str = None,
        **kwargs,
    ) -> dict:
        """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""

        start_times, end_times, file_times = create_time_stamps(
            start_time,
            start_time,
            stride=self.stride,
            hours_per_file=self.hours_per_file,
            last_file=self.last_file,
            lead_time=self.lead_time,
        )
        folder, filename = self._folder_filename(source, folder, filename)
        url = get_url(folder, filename, file_times[0])

        data = xr.open_dataset(url).isel(time=[0])

        all_points = {"lon": data.longitude.values[0], "lat": data.latitude.values[0]}

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
        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(
            start_time,
            end_time,
            stride=self.stride,
            hours_per_file=self.hours_per_file,
            last_file=self.last_file,
            lead_time=self.lead_time,
        )

        msg.info(
            f"Getting boundary spectra from WW3_4km from {self.start_time} to {self.end_time}"
        )
        bnd_list = []
        for n in range(len(file_times)):
            msg.plain(f"Reading boundary spectra: {start_times[n]}-{end_times[n]}")

            file_time = file_times[n]
            ct = 1
            keep_trying = True
            while keep_trying:
                folder, filename = self._folder_filename(source, folder, filename)
                url = get_url(folder, filename, file_time)
                try:
                    with xr.open_dataset(url) as f:
                        this_ds = f.sel(
                            time=slice(start_times[n], end_times[n]),
                            x=inds + 1,
                        )[
                            [
                                "SPEC",
                                "longitude",
                                "latitude",
                                "time",
                                "freq",
                                "direction",
                            ]
                        ].copy()
                except OSError:
                    this_ds = None

                file_consistent = self.file_is_consistent(this_ds, bnd_list, url)

                if file_consistent:
                    keep_trying = False

                elif self.data_left_to_try_with(n, ct, file_times, end_times):
                    file_time = file_times[n - ct]
                    ct += 1
                else:
                    this_ds = None
                    keep_trying = False

            # We are now out of the while loop
            if this_ds is not None:
                msg.from_file(url)
                bnd_list.append(this_ds)
                this_ds = None

        msg.info("Merging dataset together (this might take a while)...")
        bnd = xr.concat(bnd_list, dim="time").squeeze("y")

        coord_dict = {
            "lon": bnd.longitude.values,
            "lat": bnd.latitude.values,
            "time": bnd.time.values,
            "freq": bnd.freq.values,
            "dirs": bnd.direction.values,
        }

        data_dict = {"spec": bnd.SPEC.values}

        meta_dict = bnd.attrs

        return coord_dict, data_dict, meta_dict

    def file_is_consistent(self, this_ds, bnd_list, url) -> bool:
        if this_ds is None:
            msg.plain(f"SKIPPING, file not found: {url}")
            return False

        if (
            not bnd_list  # always trust the first file that is read
            or (  # for the rest, check for consistency with first file
                (this_ds.longitude == bnd_list[0].longitude).all()
                and (this_ds.latitude == bnd_list[0].latitude).all()
            )
        ):
            return True

        msg.plain(f"SKIPPING, file inconsistent: {url}")
        return False

    def data_left_to_try_with(self, n, ct, file_times, end_times) -> bool:
        if n - ct < 0:
            return False

        if pd.Timestamp(end_times[n]) - pd.Timestamp(file_times[n - ct]) > pd.Timedelta(
            self.hours_per_file, "hours"
        ):
            return False

        return True
