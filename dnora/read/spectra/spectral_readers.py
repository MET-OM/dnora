from dnora.type_manager.spectral_conventions import SpectralConvention
from dnora import msg
from dnora.utils.time import create_time_stamps
from dnora.type_manager.data_sources import DataSource
from dnora.type_manager.dnora_types import DnoraDataType
from dnora.read.abstract_readers import SpectralDataReader
from dnora.aux_funcs import get_url
from .swan_ascii import decode_lonlat, read_swan_ascii_spec
import pandas as pd
import numpy as np
from dnora.process.spectra import RemoveEmpty


import xarray as xr
import geo_parameters as gp
from geo_skeletons import PointSkeleton


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

    def __init__(self, keep_empty: bool = False):
        if keep_empty:
            self._post_processing = None
        else:
            self._post_processing = RemoveEmpty()

    def post_processing(self):
        return self._post_processing

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
        data_dict = {
            "spec": spec * 180 / np.pi
        }  # SWAN normalizes using degrees, we want normal radians normalization
        meta_dict = {"source": "Spectral wave data from a SWAN run"}

        return coord_dict, data_dict, meta_dict
