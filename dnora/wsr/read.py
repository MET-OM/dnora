import xarray as xr
import numpy as np
from copy import copy, deepcopy
import pandas as pd
from abc import ABC, abstractmethod
from typing import Union
from ..grd.grd_mod import Grid

# Import aux_funcsiliry functions
from .. import file_module
from .. import msg
from .. import aux_funcs

# from .conventions import SpectralConvention
from .wave_parameters import (
    WaveParameter,
    Dirp,
    Tp,
    Hs,
    Tm01,
    Dirm,
    Sprm,
    Tm_10,
    Tm02,
    TpI,
)
from ..spc import Spectra, process
from ..bnd.conventions import SpectralConvention


from dnora.wsr import wave_parameters
import inspect
from ..data_sources import DataSource


def dict_of_wave_parameters():
    list_of_members = inspect.getmembers(wave_parameters)
    dict_of_wps = {}
    for member in list_of_members:
        if inspect.isclass(member[1]):
            try:
                wps = [member[1]()]
            except:
                wps = None

            if wps is None:
                try:
                    wps = []
                    for n in np.linspace(-10, 10, 41):
                        wps.append(member[1](n))
                except:
                    wps = None

            if wps is None:
                try:
                    wps = []
                    for n in np.linspace(-10, 10, 41):
                        for m in np.linspace(-10, 10, 41):
                            wps.append(member[1](n, m))
                except:
                    wps = None

            if wps is not None and isinstance(wps[0], WaveParameter):
                for wp in wps:
                    dict_of_wps[wp.name()] = wp
    return dict_of_wps


def get_wave_parameter(parameter: Union[str, WaveParameter]) -> WaveParameter:
    if isinstance(parameter, str):
        return dict_of_wave_parameters().get(parameter.lower())
    return parameter


class WaveSeriesReader(ABC):
    """Reads boundary spectra from some source and provide it to the object."""

    def __init__(self):
        pass

    @abstractmethod
    def get_coordinates(
        self,
        grid: Grid,
        start_time: str,
        source: DataSource,
    ):
        """Return a list of all the available coordinated in the source.

        These are needed fo the PointPicker object to choose the relevant
        point to actually read in.

        Provide the result as four equally long nump arrays (None for missing values).

        return lon, lat, x, y
        """
        pass

    @abstractmethod
    def __call__(
        self, grid, start_time, end_time, inds, source: DataSource, **kwargs
    ) -> tuple:
        """Reads in the spectra from inds between start_time and end_time.

        The variables needed to be returned are:

        time:   Time stamps as numpy.datetime64 array
        data: dict{WaveParameter, np.ndarray [time, inds]}
        lon:    Longitude vector as numpy array (None if Cartesian)
        lat:    Latitude vector as numpy array (None if Cartesian)
        x:    Longitude vector as numpy array (None if Spherical)
        y:    Latitude vector as numpy array (None if Spherical)
        attributes: metadata: dict{key, value} will be set as attributes of the xr.Dataset

        return time, data, lon, lat, x, y, attributes
        """
        pass

    def name(self) -> str:
        return type(self).__name__

    # def __str__(self):
    # return (f"{self.start_time} - {self.end_time}")

    def default_data_source(self) -> DataSource:
        return DataSource.UNDEFINED


class SpectraToWaveSeries(WaveSeriesReader):
    """Integrates spectra to wave series"""

    def __init__(self, spectra: Spectra, freq: tuple = (0, 10_000)) -> None:
        self._spectra = deepcopy(spectra)
        try:
            self._spectra.process_spectra(process.CutFrequency(freq))
        except AttributeError:
            msg.warning(
                f"Object {self.name()} does not have a process_spectra method!\nNot cutting any frequencies!"
            )
        self._freq = freq

    def get_coordinates(
        self,
        grid,
        start_time: str,
        source: DataSource,
    ) -> tuple[np.ndarray, np.ndarray]:
        return (
            self._spectra.lon(strict=True),
            self._spectra.lat(strict=True),
            self._spectra.x(strict=True),
            self._spectra.y(strict=True),
        )
        # lon, lat, x, y = aux_funcs.get_coordinates_from_ds(self._spectra.ds())
        # return lon, lat, x, y

    # Hs(), Tp(), Dirp(), TpI(), Dirm(), Sprm(), Tm_10(), Tm01(), Tm02()
    def __call__(
        self,
        grid,
        start_time,
        end_time,
        inds,
        source: str,
        parameters: list[str] = [
            Hs(),
            Tp(),
            Dirp(),
            TpI(),
            Dirm(),
            Sprm(),
            Tm_10(),
            Tm01(),
            Tm02(),
        ],
        **kwargs,
    ) -> tuple:
        self.name = self._spectra.name
        # source = self._boundary.data.source
        time = (
            self._spectra.time(data_array=True)
            .sel(time=slice(start_time, end_time))
            .values
        )
        lon = self._spectra.lon(strict=True)
        lat = self._spectra.lat(strict=True)
        x = self._spectra.x(strict=True)
        y = self._spectra.y(strict=True)

        # efth = self._spectra.spec(data_array=True).sel(time=slice(start_time, end_time), inds=inds)

        data = {}
        for wp in parameters:
            wp = get_wave_parameter(wp)
            # data[wp] = np.swapaxes(wp(self._spectra),0,1)
            data[wp] = wp(self._spectra)
        attrs = self._spectra.ds().attrs
        attrs["integration_range"] = f"{self._freq[0]}-{self._freq[-1]} Hz"
        return time, data, lon, lat, x, y, attrs

    def name(self) -> str:
        if self._spectra is None:
            return "EmptyData"
        return self._spectra.name


class DnoraNc(WaveSeriesReader):
    def __init__(self, files: str) -> None:
        self.files = files

    def get_coordinates(
        self,
        grid,
        start_time,
        source: DataSource,
    ) -> tuple:
        data = xr.open_dataset(self.files[0]).isel(time=[0])
        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(data)
        return lon, lat, x, y

    def __call__(
        self, grid, start_time, end_time, inds, source: DataSource, **kwargs
    ) -> tuple:
        def _crop(ds):
            return ds.sel(time=slice(start_time, end_time))

        msg.info(
            f"Getting boundary spectra from DNORA type netcdf files (e.g. {self.files[0]}) from {start_time} to {end_time}"
        )

        ds = xr.open_mfdataset(self.files, preprocess=_crop, data_vars="minimal")
        ds = ds.sel(inds=inds)
        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds)
        data = {}
        for var in ds.data_vars:
            if var not in ["lon", "lat", "x", "y"]:
                if get_wave_parameter(var) is not None:
                    data[get_wave_parameter(var)] = ds.get(var).values

        return ds.get("time"), data, lon, lat, x, y, ds.attrs


class E39(WaveSeriesReader):
    def __init__(self, loc: str = "D", mode="wave"):
        self._loc = loc  # Given as "D", or "D_Breisundet"
        self._mode = mode  # 'wind' or 'wave'

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

    def _convert_var(self, var: str) -> str:
        var_dict = {
            "Hm0": "hs",
            "tm02": "tm02",
            "tp": "tp",
            "tm01": "tm01",
            "mdir": "dirm",
            "WindSpeed": "ff",
            "WindDirection": "dd",
        }  # , 'thtp': 'dirp'}
        return var_dict.get(var)

    def loc(self) -> list[str]:
        if len(self._loc) > 2:
            return self._loc
        else:
            return self._buoy_dict()[self._loc]

    def get_coordinates(
        self,
        grid,
        start_time,
        source: DataSource,
    ) -> tuple:
        start_time = pd.to_datetime(start_time)
        url = self.get_url(start_time, self.loc(), source)
        ds = xr.open_dataset(url).isel(time=[0])
        return ds.longitude.values, ds.latitude.values, None, None

    def __call__(
        self, grid, start_time, end_time, inds, source: DataSource, **kwargs
    ) -> tuple:
        # loc = np.array(self._buoys())[inds][0]
        months = aux_funcs.month_list(start_time, end_time)

        ds_list = []
        for month in months:
            url = self.get_url(pd.to_datetime(month), self.loc(), source)
            ds_list.append(xr.open_dataset(url))
        ds = xr.concat(ds_list, dim="time").sel(time=slice(start_time, end_time))
        ds["lon"] = np.nanmedian(ds.longitude.values)
        ds["lat"] = np.nanmedian(ds.latitude.values)

        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds)
        data = {}
        for var in ds.data_vars:
            if var not in ["lon", "lat", "longitude", "latitude", "x", "y"]:
                dnora_var = self._convert_var(var)
                if dnora_var is not None:
                    data[get_wave_parameter(dnora_var)] = np.swapaxes(
                        np.expand_dims(ds.get(var).values, axis=1), 0, 1
                    )
        return ds.time.values, data, lon, lat, x, y, ds.attrs

    def get_url(self, month, loc, source) -> str:
        if source == "remote":
            return (
                "https://thredds.met.no/thredds/dodsC/obs/buoy-svv-e39/"
                + month.strftime("%Y")
                + "/"
                + month.strftime("%m")
                + "/"
                + month.strftime("%Y")
                + month.strftime("%m")
                + "_E39_"
                + loc
                + f"_{self._mode}.nc"
            )


class WW3Nc(WaveSeriesReader):
    def __init__(
        self, filename: str = "ww3.%Y%m.nc", folder: str = "", mode: str = "single"
    ):
        """Mode can be 'single' (one file), 'monthly'"""
        self._filename = filename
        self._mode = mode
        self._folder = folder

    def _filenames(self, start_time, end_time, folder):
        filenames = []
        if self._mode == "single":
            filenames.append(f"{folder}/{self._filename}")
        else:
            for file in aux_funcs.month_list(start_time, end_time, fmt=self._filename):
                filenames.append(f"{folder}/{file}")
        return filenames

    def _convert_var(self, var: str) -> str:
        var_dict = {
            "hs": "hs",
            "t02": "tm02",
            "fp": "fp",
            "t0m1": "tm_10",
            "t01": "tm01",
            "spr": "sprm",
            "dir": "dirm",
            "dp": "dirp",
            "wnd": "ff",
            "wnddir": "dd",
        }
        return var_dict.get(var)

    def get_coordinates(
        self,
        grid,
        start_time,
        source: DataSource,
    ) -> tuple:
        """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""
        # day = pd.date_range(start_time, start_time,freq='D')
        filename = self._filenames(start_time, start_time, folder=self._folder)[0]
        aux_funcs.check_if_file(filename, halt=True)
        ds = xr.open_dataset(filename).isel(time=[0])
        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds)

        return lon, lat, x, y

    def __call__(
        self, grid, start_time, end_time, inds, source: DataSource, **kwargs
    ) -> tuple:
        """Reads in all wave data between the given times and at for the given indeces"""

        msg.info(
            f"Getting wave data from WW3 netcdf files from {start_time} to {end_time}"
        )

        # for file in self._filenames(start_time, end_time, self._folder):
        msg.from_file(self._filenames(start_time, end_time, self._folder))

        def _crop(ds):
            """Crop to given time and index."""

            try:  # Unstructured mesh output
                cropped_ds = ds.sel(
                    time=slice(start_time, end_time), node=inds
                ).drop_dims(["noel", "element"])
                if not self._file_identification_message_printed:
                    msg.info("Identified file as WW3 unstructured output")
                    self._file_identification_message_printed = True
                return cropped_ds
            except KeyError:
                try:  # WW3 Spectral output (get e.g. wind, hs from those files)
                    cropped_ds = ds.sel(
                        time=slice(start_time, end_time), station=(inds + 1)
                    ).drop_dims(["frequency", "direction", "string40"])
                    if not self._file_identification_message_printed:
                        msg.info("Identified file as WW3 spectral output")
                        self._file_identification_message_printed = True
                    return cropped_ds
                except:
                    raise KeyError(
                        "Could not identify the file as a known format (WW3 Unstructured or Spectral output)!"
                    )

        import dask

        with dask.config.set(**{"array.slicing.split_large_chunks": True}):
            self._file_identification_message_printed = False
            with xr.open_mfdataset(
                self._filenames(start_time, end_time, self._folder), preprocess=_crop
            ) as ds:
                lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds)

                data = {}
                for var in ds.data_vars:
                    # if var not in ['lon', 'lat', 'longitude', 'latitude', 'x', 'y']:
                    dnora_var = self._convert_var(var)
                    if dnora_var is not None:
                        # if dnora_var == 'fp':
                        #    data[get_wave_parameter('tp')] = np.swapaxes(ds.get(var).values**-1,0,1)
                        # else:
                        msg.info(f"Importing {dnora_var} << {var}")
                        data[get_wave_parameter(dnora_var)] = np.swapaxes(
                            ds.get(var).values, 0, 1
                        )
                    else:
                        msg.info(
                            f"Could not identify {var} as a DNORA WaveSeries variable"
                        )

                return ds.time.values, data, lon, lat, x, y, ds.attrs


class WW3Nc_old(WaveSeriesReader):
    def __init__(self, filename: str, parameters: list[str] = ["hs", "tp"]):
        self._filename = filename

    def _convert_var(self, var: str) -> str:
        var_dict = {
            "hs": "hs",
            "t02": "tm02",
            "fp": "fp",
            "t0m1": "tm_10",
            "t01": "tm01",
            "spr": "sprm",
            "dir": "dirm",
            "dp": "dirp",
        }
        return var_dict.get(var)

    def get_coordinates(
        self,
        grid,
        start_time,
        source: DataSource,
    ) -> tuple:
        start_time = pd.to_datetime(start_time)
        ds = xr.open_dataset(self._filename).isel(time=[0])
        return ds.longitude.values, ds.latitude.values, None, None

    def __call__(
        self, grid, start_time, end_time, inds, source: DataSource, **kwargs
    ) -> tuple:
        ds = xr.open_dataset(self._filename).sel(node=inds)

        ds["lon"] = ds.longitude.values
        ds["lat"] = ds.latitude.values
        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds)

        data = {}
        for var in ds.data_vars:
            if var not in ["lon", "lat", "longitude", "latitude", "x", "y"]:
                dnora_var = self._convert_var(var)
                if dnora_var is not None:
                    if dnora_var == "fp":
                        data[get_wave_parameter("tp")] = np.swapaxes(
                            ds.get(var).values ** -1, 0, 1
                        )
                    else:
                        data[get_wave_parameter(dnora_var)] = np.swapaxes(
                            ds.get(var).values, 0, 1
                        )
        return ds.time.values, data, lon, lat, x, y, ds.attrs
