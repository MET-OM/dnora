import xarray as xr
import numpy as np
from copy import copy
from abc import ABC, abstractmethod

# Import aux_funcsiliry functions
from dnora import msg, aux_funcs
from dnora.spectral_conventions import SpectralConvention
import pandas as pd
from dnora.data_sources import DataSource


# class ForceFeed(SpectraReader):
#     def __init__(
#         self, time, freq, dirs, spec, lon, lat, convention: SpectralConvention
#     ) -> None:
#         self.time = copy(time)
#         self.freq = copy(freq)
#         self.dirs = copy(dirs)
#         self.spec = copy(spec)
#         self.lon = copy(lon)
#         self.lat = copy(lat)
#         self.set_convention(convention)
#         return

#     def get_coordinates(self, grid, start_time, source, folder) -> tuple:
#         return copy(self.lon), copy(self.lat)

#     def __call__(
#         self, grid, start_time, end_time, inds, source, folder, **kwargs
#     ) -> tuple:
#         # return  copy(self.time), copy(self.freq), copy(self.dirs), np.reshape(self.spec, (len(self.time), len(self.lon), self.spec.shape[0], self.spec.shape[1])), copy(self.lon), copy(self.lat), ''
#         return (
#             copy(self.time),
#             copy(self.freq),
#             copy(self.dirs),
#             copy(self.spec),
#             copy(self.lon),
#             copy(self.lat),
#             None,
#             None,
#             {},
#         )


# class WW3Nc(SpectraReader):
#     def __init__(
#         self, filename: str = "ww3.%Y%m_spec.nc", folder: str = "", mode: str = "single"
#     ):
#         """Mode can be 'single' (one file), 'monthly'"""
#         self._filename = filename
#         self._mode = mode
#         self._folder = folder

#     def convention(self) -> SpectralConvention:
#         return SpectralConvention.WW3

#     def _filenames(self, start_time, end_time, folder):
#         filenames = []
#         if self._mode == "single":
#             filenames.append(f"{folder}/{self._filename}")
#         else:
#             for file in aux_funcs.month_list(start_time, end_time, fmt=self._filename):
#                 filenames.append(f"{folder}/{file}")
#         return filenames

#     def get_coordinates(self, grid, start_time, source, folder) -> tuple:
#         """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""
#         # day = pd.date_range(start_time, start_time,freq='D')

#         data = xr.open_dataset(
#             self._filenames(start_time, start_time, folder=self._folder)[0]
#         ).isel(time=[0])

#         lon_all = data.longitude.values[0]
#         lat_all = data.latitude.values[0]

#         return lon_all, lat_all, None, None

#     def __call__(
#         self, grid, start_time, end_time, inds, source, folder, **kwargs
#     ) -> tuple:
#         """Reads in all boundary spectra between the given times and at for the given indeces"""

#         msg.info(
#             f"Getting boundary spectra from WW3 netcdf files from {start_time} to {end_time}"
#         )

#         def _crop(ds):
#             """
#             EMODNET tiles overlap by two cells on each boundary.
#             """
#             return ds.sel(time=slice(start_time, end_time), station=(inds + 1))

#         import dask

#         with dask.config.set(**{"array.slicing.split_large_chunks": True}):
#             with xr.open_mfdataset(
#                 self._filenames(start_time, end_time, self._folder), preprocess=_crop
#             ) as ds:
#                 time = ds.time.values
#                 freq = ds.frequency.values
#                 dirs = ds.direction.values
#                 spec = ds.efth.values
#                 lon = ds.longitude.values[0, :]
#                 lat = ds.latitude.values[0, :]

#                 return time, freq, dirs, spec, lon, lat, None, None, ds.attrs


# class File_WW3Nc(BoundaryReader):
#     def __init__(self, folder: str='', filename: str='ww3_T0', dateftm: str='%Y%m%dT%H%M', stride: int=6, hours_per_file: int=73, last_file: str='', lead_time: int=0) -> None:
#         self.stride = copy(stride)
#         self.hours_per_file = copy(hours_per_file)
#         self.lead_time = copy(lead_time)
#         self.last_file = copy(last_file)

#         if (not folder == '') and (not folder[-1] == '/'):
#             self.folder = folder + '/'
#         else:
#             self.folder = copy(folder)

#         self.filestring = copy(filename)
#         self.datestring = copy(dateftm)

#         return

#     def convention(self) -> SpectralConvention:
#         return SpectralConvention.WW3

#     def get_coordinates(self, grid, start_time) -> tuple:
#         """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""
#         #day = pd.date_range(start_time, start_time,freq='D')
#         start_times, end_times, file_times = aux_funcs.create_time_stamps(start_time, start_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)
#         filename = self.get_filename(file_times[0])

#         data = xr.open_dataset(filename).isel(time = [0])

#         lon_all = data.longitude.values[0]
#         lat_all = data.latitude.values[0]

#         return lon_all, lat_all

#     def __call__(self, grid, start_time, end_time, inds, **kwargs) -> tuple:
#         """Reads in all boundary spectra between the given times and at for the given indeces"""
#         self.start_time = start_time
#         self.end_time = end_time

#         start_times, end_times, file_times = aux_funcs.create_time_stamps(start_time, end_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)

#         msg.info(f"Getting boundary spectra from NORA3 from {self.start_time} to {self.end_time}")
#         bnd_list = []
#         for n in range(len(file_times)):
#             filename = self.get_filename(file_times[n])
#             msg.from_file(filename)
#             msg.plain(f"Reading boundary spectra: {start_times[n]}-{end_times[n]}")

#             bnd_list.append(xr.open_dataset(filename).sel(time = slice(start_times[n], end_times[n]), station = (inds+1)))

#         bnd=xr.concat(bnd_list, dim="time")


#         # for x in range(len(inds)):
#         #     for t in range(len(bnd.time.values)):
#         #         new_spec, new_dirs = WW3ToOcean()(bnd.efth.values[t,x,:,:],bnd.direction.values)
#         #         bnd.efth.values[t,x,:,:] = new_spec

#         time = bnd.time.values
#         freq = bnd.frequency.values
#         dirs = bnd.direction.values
#         spec = bnd.efth.values
#         lon = bnd.longitude.values[0,:]
#         lat = bnd.latitude.values[0,:]

#         source = f"ww3_ouput_spectra"

#         return  time, freq, dirs, spec, lon, lat, None, None, bnd.attrs


#     def get_filename(self, time) -> str:
#         filename = self.folder + file_module.replace_times(self.filename,
#                                                         self.dateformat,
#                                                         [time]) + '.nc'
#         return filename
