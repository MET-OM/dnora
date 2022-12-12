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
from .conventions import SpectralConvention
import pandas as pd
class BoundaryReader(ABC):
    """Reads boundary spectra from some source and provide it to the object."""

    @abstractmethod
    def get_coordinates(self, grid, start_time):
        """Return a list of all the available coordinated in the source.

        These are needed fo the PointPicker object to choose the relevant
        point to actually read in.

        Provide the result as two equally long nump arrays.
        """
        return lon, lat, x, y

    @abstractmethod
    def convention(self) -> SpectralConvention:
        """Return the convention of the spectra returned to the object.

        The conventions to choose from are predetermined:

        OCEAN:    Oceanic convention
                    Directional vector monotonically increasing.
                    Direction to. North = 0, East = 90.

        MET:      Meteorological convention
                    Directional vector monotonically increasing.
                    Direction from. North = 0, East = 90.

        MATH:     Mathematical convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        MATHVEC:  Mathematical convention in vector
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        WW3:      WAVEWATCH III output convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 0, East = 90.
        """
        return convention

    @abstractmethod
    def __call__(self, grid, start_time, end_time, inds, **kwargs) -> Tuple:
        """Reads in the spectra from inds between start_time and end_time.

        The variables needed to be returned are:

        time:   Time stamps as numpy.datetime64 array
        freq:   Frequency vector as numpy array
        dirs:   Directional vector as numpy array
        spec:   Boundary spectra [time, station, freq, dirs] as numpy array
        lon:    Longitude vector as numpy array (None if Cartesian)
        lat:    Latitude vector as numpy array (None if Cartesian)
        x:    Longitude vector as numpy array (None if Spherical)
        y:    Latitude vector as numpy array (None if Spherical)
        metadata: metadata: dict{key, value} will be set as attributes of the xr.Dataset
        """

        return time, freq, dirs, spec, lon, lat, x, y, metadata


    def name(self) -> str:
        return type(self).__name__

    def set_convention(self, convention: SpectralConvention) -> None:
        if isinstance(convention, str):
            self._convention = SpectralConvention[convention.upper()]
        else:
            self._convention = convention

    def convention(self) -> SpectralConvention:
        return self._convention

    def set_source(self, source: str) -> None:
        self._source = source

    def source(self) -> str:
        if hasattr(self, '_source'):
            return self._source
        return 'remote'

    def post_processing(self):
        return None
    #def __str__(self):
        #return (f"{self.start_time} - {self.end_time}")

class ConstantBoundary(BoundaryReader):
    def __init__(self, grid: Grid, spec: float=1, cartesian: bool=False, metadata: dict=None, spectral_convention: SpectralConvention=SpectralConvention.OCEAN):
        self.spec = spec
        self.metadata = metadata
        self.cartesian = cartesian
        self.set_convention(spectral_convention)
        self.grid = grid

    def get_coordinates(self, grid, start_time) -> tuple:
        lon, lat, x, y = aux_funcs.get_coordinates_from_grid(self.grid, self.cartesian, list=True)

        return lon, lat, x, y

    def __call__(self, grid, start_time, end_time, inds, **kwargs):
        time = pd.date_range(start=start_time, end=end_time, freq='H').values
        lon, lat, x, y = aux_funcs.get_coordinates_from_grid(self.grid, self.cartesian, list=True)
        freq = np.array(range(1,11))/10.

        if self.convention() in [SpectralConvention.WW3, SpectralConvention.MATHVEC]:
            dirs = np.mod(np.linspace(90.,-255.,24),360)
        else:
            dirs = np.linspace(0.,345.,24)

        if self.convention() in [SpectralConvention.MATH, SpectralConvention.MATHVEC]:
            north = 90
        elif self.convention() in [SpectralConvention.MET]:
            north = 180
        else:
            north = 0
        #dirs = np.array(range(0,360,15))

        spec = np.full((len(time), len(inds), len(freq), len(dirs)), 0)
        north_ind = np.where(dirs==north)[0][0]
        spec[:,:,:,north_ind] = np.full((len(time), len(inds), len(freq)), self.spec)
        metadata = {'metadata': 'this is a constant boundary'}

        return time, freq, dirs, spec, lon, lat, x, y, metadata


class DnoraNc(BoundaryReader):
    def __init__(self, files: str) -> None:
        self.files = files

    def get_coordinates(self, grid, start_time) -> Tuple:
        data = xr.open_dataset(self.files[0]).isel(time = [0])
        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(data)
        return lon, lat, x, y

    def __call__(self, grid, start_time, end_time, inds, **kwargs) -> Tuple:
        def _crop(ds):
            return ds.sel(time=slice(start_time, end_time))
        msg.info(f"Getting boundary spectra from DNORA type netcdf files (e.g. {self.files[0]}) from {start_time} to {end_time}")
        ds = xr.open_mfdataset(self.files, preprocess=_crop, data_vars='minimal')
        ds = ds.sel(inds=inds)
        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds)
        if not hasattr(self, '_convention'):
            self.set_convention(ds.spectral_convention)

        return ds.time.values, ds.freq.values, ds.dirs.values, ds.spec.values, lon, lat, x, y, ds.attrs

class ForceFeed(BoundaryReader):
    def __init__(self, time, freq, dirs, spec, lon, lat, convention: SpectralConvention) -> None:
        self.time = copy(time)
        self.freq = copy(freq)
        self.dirs = copy(dirs)
        self.spec = copy(spec)
        self.lon = copy(lon)
        self.lat = copy(lat)
        self.set_convention(convention)
        return

    def get_coordinates(self, grid, start_time) -> Tuple:
        return copy(self.lon), copy(self.lat)

    def __call__(self, grid, start_time, end_time, inds, **kwargs) -> Tuple:
        #return  copy(self.time), copy(self.freq), copy(self.dirs), np.reshape(self.spec, (len(self.time), len(self.lon), self.spec.shape[0], self.spec.shape[1])), copy(self.lon), copy(self.lat), ''
        return  copy(self.time), copy(self.freq), copy(self.dirs), copy(self.spec), copy(self.lon), copy(self.lat), None, None, {}

class File_WW3Nc(BoundaryReader):
    def __init__(self, folder: str='', filename: str='ww3_T0', dateftm: str='%Y%m%dT%H%M', stride: int=6, hours_per_file: int=73, last_file: str='', lead_time: int=0) -> None:
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)

        if (not folder == '') and (not folder[-1] == '/'):
            self.folder = folder + '/'
        else:
            self.folder = copy(folder)

        self.filestring = copy(filename)
        self.datestring = copy(dateftm)

        return

    def convention(self) -> SpectralConvention:
        return SpectralConvention.WW3

    def get_coordinates(self, grid, start_time) -> Tuple:
        """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""
        #day = pd.date_range(start_time, start_time,freq='D')
        start_times, end_times, file_times = aux_funcs.create_time_stamps(start_time, start_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)
        filename = self.get_filename(file_times[0])

        data = xr.open_dataset(filename).isel(time = [0])

        lon_all = data.longitude.values[0]
        lat_all = data.latitude.values[0]

        return lon_all, lat_all

    def __call__(self, grid, start_time, end_time, inds, **kwargs) -> Tuple:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = aux_funcs.create_time_stamps(start_time, end_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)

        msg.info(f"Getting boundary spectra from NORA3 from {self.start_time} to {self.end_time}")
        bnd_list = []
        for n in range(len(file_times)):
            filename = self.get_filename(file_times[n])
            msg.from_file(filename)
            msg.plain(f"Reading boundary spectra: {start_times[n]}-{end_times[n]}")

            bnd_list.append(xr.open_dataset(filename).sel(time = slice(start_times[n], end_times[n]), station = (inds+1)))

        bnd=xr.concat(bnd_list, dim="time")


        # for x in range(len(inds)):
        #     for t in range(len(bnd.time.values)):
        #         new_spec, new_dirs = WW3ToOcean()(bnd.efth.values[t,x,:,:],bnd.direction.values)
        #         bnd.efth.values[t,x,:,:] = new_spec

        time = bnd.time.values
        freq = bnd.frequency.values
        dirs = bnd.direction.values
        spec = bnd.efth.values
        lon = bnd.longitude.values[0,:]
        lat = bnd.latitude.values[0,:]

        source = f"ww3_ouput_spectra"

        return  time, freq, dirs, spec, lon, lat, None, None, bnd.attrs


    def get_filename(self, time) -> str:
        filename = self.folder + file_module.replace_times(self.filename,
                                                        self.dateformat,
                                                        [time]) + '.nc'
        return filename
