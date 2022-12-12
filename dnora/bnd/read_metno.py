import os
import xarray as xr
import glob
import numpy as np
from copy import copy
from abc import ABC, abstractmethod
from typing import Tuple
import pandas as pd
# Import abstract classes and needed instances of them
from .read import BoundaryReader
from .process import RemoveEmpty
from .conventions import SpectralConvention
# Import aux_funcsiliry functions
from .. import msg
from ..aux_funcs import create_time_stamps

class WAM4km(BoundaryReader):
    def __init__(self, ignore_nan: bool=True, stride: int=6,
                 hours_per_file: int=73, last_file: str='', lead_time: int=0) -> None:
        self.ignore_nan = copy(ignore_nan)
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)

    def convention(self) -> str:
        return SpectralConvention.OCEAN

    def post_processing(self):
        return RemoveEmpty()

    def get_coordinates(self, grid, start_time) -> Tuple:
        """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""

        start_times, end_times, file_times = create_time_stamps(start_time, start_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)
        url = self.get_url(file_times[0])

        data = xr.open_dataset(url).isel(time = [0])

        lon_all = data.longitude.values[0]
        lat_all = data.latitude.values[0]

        return lon_all, lat_all, None, None

    def __call__(self, grid, start_time, end_time, inds) -> Tuple:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time


        start_times, end_times, file_times = create_time_stamps(start_time, end_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)


        msg.info(f"Getting boundary spectra from WAM4 from {self.start_time} to {self.end_time}")
        bnd_list = []
        for n in range(len(file_times)):


            msg.plain(f"Reading boundary spectra: {start_times[n]}-{end_times[n]}")

            file_time = file_times[n]
            ct = 1
            keep_trying = True
            while keep_trying:

                url = self.get_url(file_time)


                try:
                    with xr.open_dataset(url) as f:
                        this_ds = f.sel(
                            time=slice(start_times[n], end_times[n]),
                            x=inds+1,
                        )[
                            ['SPEC', 'longitude', 'latitude', 'time', 'freq', 'direction']
                        ].copy()
                except OSError:
                        this_ds = None

                file_consistent = self.file_is_consistent(this_ds, bnd_list, url)

                if file_consistent:
                    keep_trying =False

                elif self.data_left_to_try_with(n, ct, file_times, end_times):
                    file_time = file_times[n-ct]
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
        bnd = xr.concat(bnd_list, dim="time").squeeze('y')

        time = bnd.time.values
        freq = bnd.freq.values
        dirs = bnd.direction.values
        spec = bnd.SPEC.values
        #spec = np.moveaxis(spec,1,0) # time, inds ... -> inds, time, ...
        lon = bnd.longitude.values
        lat = bnd.latitude.values

        return  time, freq, dirs, spec, lon, lat, None, None, bnd.attrs

    def file_is_consistent(self, this_ds, bnd_list, url) -> bool:
        if this_ds is None:
            msg.plain(f'SKIPPING, file not found: {url}')
            return False

        if (
            not bnd_list # always trust the first file that is read
            or ( # for the rest, check for consistency with first file
                (this_ds.longitude == bnd_list[0].longitude).all() and
                (this_ds.latitude  == bnd_list[0].latitude ).all()

            )):
            return True

        msg.plain(f'SKIPPING, file inconsistent: {url}')
        return False

    def data_left_to_try_with(self, n, ct, file_times, end_times) -> bool:
        if n-ct<0:
            return False

        if pd.Timestamp(end_times[n])-pd.Timestamp(file_times[n-ct])>pd.Timedelta(self.hours_per_file, 'hours'):
            return False

        return True

    def get_url(self, day):
        if self.source() == 'remote':
            return 'https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/'+day.strftime('%d')+'/MyWave_wam4_SPC_'+day.strftime('%Y%m%d')+'T'+day.strftime('%H')+'Z.nc'
        if self.source() == 'met':
            return '/lustre/storeB/project/fou/om/xxxxxxxxx/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/MyWave_wam4_SPC_'+day.strftime('%Y%m%d')
        else:
            return self.source()+'/MyWave_wam4_SPC_'+day.strftime('%Y%m%d')

class WAM4kmOld(BoundaryReader):
    def __init__(self, ignore_nan: bool=True, stride: int=6,
                 hours_per_file: int=73, last_file: str='', lead_time: int=0) -> None:
        self.ignore_nan = copy(ignore_nan)
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)

    def convention(self) -> str:
        return SpectralConvention.OCEAN


    def get_coordinates(self, grid, start_time) -> Tuple:
        """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""

        start_times, end_times, file_times = create_time_stamps(start_time, start_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)
        url = self.get_url(file_times[0])

        data = xr.open_dataset(url).isel(time = [0])

        lon_all = data.longitude.values[0]
        lat_all = data.latitude.values[0]

        if self.ignore_nan:
            msg.info('ignore_nan = True. The following points are NOT provided to the point_picker')
            mask = [True]*len(lon_all)
            for n in range(len(lon_all)):
                if np.isnan(data.SPEC.values[:,0,n,:,:]).any():
                    msg.plain(f"{lon_all[n]:10.7f}, {lat_all[n]:10.7f}")
                    mask[n] = False

            lon_all = lon_all[mask]
            lat_all = lat_all[mask]
            self.pointers = np.where(mask) # These are used to map back the indeces from the PointPicker to indeces in the original data

        return lon_all, lat_all, None, None

    def __call__(self, grid, start_time, end_time, inds) -> Tuple:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        # If we have removed some spectra (NaN's) we need to remap the indeces
        if hasattr(self, 'pointers'):
            inds = self.pointers[0][inds]

        start_times, end_times, file_times = create_time_stamps(start_time, end_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)


        msg.info(f"Getting boundary spectra from WAM4 from {self.start_time} to {self.end_time}")
        bnd_list = []
        for n in range(len(file_times)):


            msg.plain(f"Reading boundary spectra: {start_times[n]}-{end_times[n]}")

            file_time = file_times[n]
            ct = 1
            keep_trying = True
            while keep_trying:

                url = self.get_url(file_time)


                try:
                    with xr.open_dataset(url) as f:
                        this_ds = f.sel(
                            time=slice(start_times[n], end_times[n]),
                            x=inds+1,
                        )[
                            ['SPEC', 'longitude', 'latitude', 'time', 'freq', 'direction']
                        ].copy()
                except OSError:
                        this_ds = None

                file_consistent = self.file_is_consistent(this_ds, bnd_list, url)

                if file_consistent:
                    keep_trying =False

                elif self.data_left_to_try_with(n, ct, file_times, end_times):
                    file_time = file_times[n-ct]
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
        bnd = xr.concat(bnd_list, dim="time").squeeze('y')

        time = bnd.time.values
        freq = bnd.freq.values
        dirs = bnd.direction.values
        spec = bnd.SPEC.values
        #spec = np.moveaxis(spec,1,0) # time, inds ... -> inds, time, ...
        lon = bnd.longitude.values
        lat = bnd.latitude.values

        return  time, freq, dirs, spec, lon, lat, None, None, bnd.attrs

    def file_is_consistent(self, this_ds, bnd_list, url) -> bool:
        if this_ds is None:
            msg.plain(f'SKIPPING, file not found: {url}')
            return False

        if (
            not bnd_list # always trust the first file that is read
            or ( # for the rest, check for consistency with first file
                (this_ds.longitude == bnd_list[0].longitude).all() and
                (this_ds.latitude  == bnd_list[0].latitude ).all()

            )):
            return True

        msg.plain(f'SKIPPING, file inconsistent: {url}')
        return False

    def data_left_to_try_with(self, n, ct, file_times, end_times) -> bool:
        if n-ct<0:
            return False

        if pd.Timestamp(end_times[n])-pd.Timestamp(file_times[n-ct])>pd.Timedelta(self.hours_per_file, 'hours'):
            return False

        return True

    def get_url(self, day):
        if self.source() == 'remote':
            return 'https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/'+day.strftime('%d')+'/MyWave_wam4_SPC_'+day.strftime('%Y%m%d')+'T'+day.strftime('%H')+'Z.nc'
        if self.source() == 'met':
            return '/lustre/storeB/project/fou/om/xxxxxxxxx/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/MyWave_wam4_SPC_'+day.strftime('%Y%m%d')
        else:
            return self.source()+'/MyWave_wam4_SPC_'+day.strftime('%Y%m%d')



class NORA3(BoundaryReader):
    def __init__(self, stride: int=24, hours_per_file: int=24,
                last_file: str='', lead_time: int=0) -> None:
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)

    def convention(self) -> str:
        return SpectralConvention.OCEAN

    def get_coordinates(self, grid, start_time) -> Tuple:
        """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""
        start_times, end_times, file_times = create_time_stamps(start_time, start_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)
        url = self.get_url(file_times[0])
        data = xr.open_dataset(url).isel(time = [0])

        lon_all = data.longitude.values[0]
        lat_all = data.latitude.values[0]

        return lon_all, lat_all, None, None

    def __call__(self, grid, start_time, end_time, inds) -> Tuple:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(start_time, end_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)

        msg.info(f"Getting boundary spectra from NORA3 from {self.start_time} to {self.end_time}")
        bnd_list = []
        for n in range(len(file_times)):
            url = self.get_url(file_times[n])
            msg.from_file(url)
            msg.plain(f"Reading boundary spectra: {start_times[n]}-{end_times[n]}")
            with xr.open_dataset(url) as f:
                this_ds = f.sel(time = slice(start_times[n], end_times[n]), x = (inds+1))[['SPEC', 'longitude', 'latitude', 'time', 'freq', 'direction']].copy()
            bnd_list.append(this_ds)
        bnd=xr.concat(bnd_list, dim="time").squeeze('y')

        time = bnd.time.values
        freq = bnd.freq.values
        dirs = bnd.direction.values
        spec = bnd.SPEC.values
        #spec = np.moveaxis(spec,1,0) # time, inds ... -> inds, time, ...
        lon = bnd.longitude.values[0,:]
        lat = bnd.latitude.values[0,:]
        bnd.attrs.pop('direction_convention')

        return  time, freq, dirs, spec, lon, lat, None, None, bnd.attrs


    def get_url(self, day) -> str:
        if self.source() == 'remote':
            return 'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/SPC'+day.strftime('%Y%m%d')+'00.nc'
        if self.source() == 'met':
            return '/lustre/storeB/project/fou/om/WINDSURFER/mw3hindcast/spectra/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/SPC'+day.strftime('%Y%m%d')+'00.nc'
        else:
            return self.source()+'/SPC'+day.strftime('%Y%m%d')+'00.nc'
