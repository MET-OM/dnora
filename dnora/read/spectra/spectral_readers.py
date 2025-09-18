from dnora.type_manager.spectral_conventions import SpectralConvention
from dnora import msg
from dnora.utils.time import create_time_stamps
from dnora.type_manager.data_sources import DataSource
from dnora.type_manager.dnora_types import DnoraDataType
from dnora.read.abstract_readers import SpectralDataReader
from .swan_ascii import decode_lonlat, read_swan_ascii_spec
import pandas as pd
import numpy as np
from dnora.process.spectra import RemoveEmpty
import glob
from dnora.utils.spec import expand_to_directional_spectrum
import xarray as xr
import geo_parameters as gp
from geo_skeletons import PointSkeleton
from dnora.utils.io import get_url
from dnora.read.generic import PointNetcdf


class SWAN_Nc(PointNetcdf):
    def convention(self) -> str:
        return SpectralConvention.MET

    def __call__(self, *args, **kwargs):
        ds = super().__call__(*args, **kwargs)
        theta = ds.dirs.values
        dd = np.mod(np.rad2deg(theta) + 360, 360)
        ds["dirs"] = np.round(dd, 2)
        ds = ds.sortby("dirs")
        return ds


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

        folder = self._folder(folder, source)
        filename = self._filename(filename, source)

        url = get_url(folder, filename)
        url = glob.glob(url)
        if len(url) > 1:
            raise ValueError(f"Cannot read several ascii files at once!")
        else:
            url = url[0]

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
        folder = self._folder(folder, source)
        filename = self._filename(filename, source)
        url = get_url(folder, filename)
        url = glob.glob(url)
        if len(url) > 1:
            raise ValueError(f"Cannot read several ascii files at once!")
        else:
            url = url[0]
        time, lon, lat, spec, freq, dirs = read_swan_ascii_spec(
            url, start_time=start_time, end_time=end_time
        )

        coord_dict = {
            "lon": lon[inds],
            "lat": lat[inds],
            "time": time,
            "freq": freq,
            "dirs": dirs,
        }

        data_dict = {
            "spec": spec[:, inds, :, :] * 180 / np.pi
        }  # SWAN normalizes using degrees, we want normal radians normalization
        meta_dict = {"source": "Spectral wave data from a SWAN run"}

        return coord_dict, data_dict, meta_dict


class Spectra1DToSpectra(SpectralDataReader):
    """Integrates boundary spectra to omnidairectional spectra"""

    def __init__(self, spectra1d, dirs, dirp=None) -> None:
        self._spectra1d = spectra1d
        self._dirs = dirs
        self._dirp = dirp

    def convention(self):
        return SpectralConvention.MET

    def post_processing(self):
        return RemoveEmpty(0)

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
            "lon": self._spectra1d.lon(strict=True),
            "lat": self._spectra1d.lat(strict=True),
            "x": self._spectra1d.x(strict=True),
            "y": self._spectra1d.y(strict=True),
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
        **kwargs,
    ) -> tuple:
        time = (
            self._spectra1d.time(data_array=True)
            .sel(time=slice(start_time, end_time))
            .values
        )
        lon = self._spectra1d.lon(strict=True)
        lat = self._spectra1d.lat(strict=True)
        x = self._spectra1d.x(strict=True)
        y = self._spectra1d.y(strict=True)
        freq = self._spectra1d.freq()
        dirs = self._dirs
        obj = dnora_class(time=time, lon=lon, lat=lat, x=x, y=y, freq=freq, dirs=dirs)

        spec1d = self._spectra1d.spec(squeeze=False)
        dirp = self._dirp
        msg.plain(
            "Calculating directional spectrum from Dirp with spreading using a cos**(2s) distribution..."
        )
        spec2d = expand_to_directional_spectrum(spec1d, freq=freq, dirs=dirs, dirp=dirp)

        obj.set_spec(spec2d)
        obj = obj.sel(inds=inds)
        return obj.ds()

    def name(self):
        readername = "cos2s"
        if self._spectra1d is not None:
            readername += f"_{self._spectra1d.name}"
        return readername
