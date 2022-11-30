import xarray as xr
import numpy as np
from copy import copy
from abc import ABC, abstractmethod
from typing import Tuple
from ..grd.grd_mod import Grid
# Import aux_funcsiliry functions
from .. import file_module
from .. import msg
from .. import aux_funcs
#from .conventions import SpectralConvention
from .wave_parameters import WaveParameter, Hs, Tm01, Tp, Dirm, Sprm
from ..spc import Spectra
from ..bnd.conventions import SpectralConvention



from dnora.wsr import wave_parameters
import inspect

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
                    for n in np.linspace(-10,10,41):
                        wps.append(member[1](n))
                except:
                    wps = None

            if wps is None:
                try:
                    wps = []
                    for n in np.linspace(-10,10,41):
                        for m in np.linspace(-10,10,41):
                            wps.append(member[1](n,m))
                except:
                    wps = None


            if wps is not None and isinstance(wps[0], WaveParameter):
                for wp in wps:
                    dict_of_wps[wp.name()] = wp
    return dict_of_wps

class WaveSeriesReader(ABC):
    """Reads boundary spectra from some source and provide it to the object."""
    def __init__(self):
        pass

    @abstractmethod
    def get_coordinates(self, grid, start_time):
        """Return a list of all the available coordinated in the source.

        These are needed fo the PointPicker object to choose the relevant
        point to actually read in.

        Provide the result as two equally long nump arrays.
        """
        return lon, lat

    @abstractmethod
    def __call__(self, grid, start_time, end_time, inds) -> Tuple:
        """Reads in the spectra from inds between start_time and end_time.

        The variables needed to be returned are:

        time:   Time stamps as numpy.datetime64 array
        data: dict{WaveParameter, np.ndarray [inds, time]}
        lon:    Longitude vector as numpy array (None if Cartesian)
        lat:    Latitude vector as numpy array (None if Cartesian)
        x:    Longitude vector as numpy array (None if Spherical)
        y:    Latitude vector as numpy array (None if Spherical)
        attributes: metadata: dict{key, value} will be set as attributes of the xr.Dataset
        """

        return time, data, lon, lat, x, y, attributes

    def name(self) -> str:
        return type(self).__name__

    def set_source(self, source: str) -> None:
        self._source = source

    def source(self) -> str:
        if hasattr(self, '_source'):
            return self._source
        return 'remote'

    #def __str__(self):
        #return (f"{self.start_time} - {self.end_time}")

class SpectraToWaveSeries(WaveSeriesReader):
    """Integrates spectra to wave series"""
    def __init__(self, spectra: Spectra) -> None:
        self._spectra = copy(spectra)

    def get_coordinates(self, grid, start_time: str) -> Tuple[np.ndarray, np.ndarray]:
        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(self._spectra.ds())
        return lon, lat, x, y

    def __call__(self, grid, start_time, end_time, inds) -> Tuple:
        self.name = self._spectra.name
        #source = self._boundary.data.source
        time = self._spectra.time(data_array=True).sel(time=slice(start_time, end_time)).values
        lon = self._spectra.lon(strict=True)
        lat = self._spectra.lat(strict=True)
        x = self._spectra.x(strict=True)
        y = self._spectra.y(strict=True)


        efth = self._spectra.spec(data_array=True).sel(time=slice(start_time, end_time), inds=inds)

        parameters = [Hs(), Dirm(), Tp(), Sprm()]
        data = {}
        for wp in parameters:
            data[wp] = wp(self._spectra)
        return time, data, lon, lat, x, y, self._spectra.ds().attrs

class DnoraNc(WaveSeriesReader):
    def __init__(self, files: str) -> None:
        self.files = files

    def get_coordinates(self, grid, start_time) -> tuple:
        data = xr.open_dataset(self.files[0]).isel(time = [0])
        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(data)
        return lon, lat, x, y

    def __call__(self, grid, start_time, end_time, inds) -> tuple:
        def _crop(ds):
            return ds.sel(time=slice(start_time, end_time))
        msg.info(f"Getting boundary spectra from DNORA type netcdf files (e.g. {self.files[0]}) from {start_time} to {end_time}")
        ds = xr.open_mfdataset(self.files, preprocess=_crop, data_vars='minimal')
        ds = ds.sel(inds=inds)
        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds)
        dict_of_wps = dict_of_wave_parameters()
        data = {}

        for var in ds.data_vars:
            if var not in ['lon', 'lat', 'x', 'y']:
                data[dict_of_wps[var]] = ds.get(var).values

        return ds.get('time'), data, lon, lat, x, y, ds.attrs
