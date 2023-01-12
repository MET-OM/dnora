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

# Import aux_funcsiliry functions
from .. import msg
from ..aux_funcs import day_list, create_time_stamps

class WAM4km(BoundaryReader):
    def __init__(self, ignore_nan: bool=True, stride: int=6,
                 hours_per_file: int=73, last_file: str='', lead_time: int=0,
                 cache: bool=True, clean_cache: bool=False) -> None:
        self.ignore_nan = copy(ignore_nan)
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)
        self.cache = copy(cache)
        self.clean_cache = copy(clean_cache)
        self.cache_folder = 'dnora_bnd_temp'
        return

    def convention(self) -> str:
        return 'Ocean'


    def get_coordinates(self, start_time) -> Tuple:
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

        return lon_all, lat_all

    def __call__(self, start_time, end_time, inds) -> Tuple:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        # If we have removed some spectra (NaN's) we need to remap the indeces
        if hasattr(self, 'pointers'):
            inds = self.pointers[0][inds]

        if self.clean_cache:
            for f in glob.glob(os.path.join(self.cache_folder, '*')):
                os.remove(f)

        if self.cache:
            if not os.path.isdir(self.cache_folder):
                os.mkdir(self.cache_folder)
                print("Creating cache folder %s..." % self.cache_folder)

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
                url_or_cache, data_from_cache = self.get_filepath_if_cached(url)
                if data_from_cache and ct == 1:
                    this_ds = xr.load_dataset(url_or_cache)
                else:
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
                msg.from_file(url_or_cache)
                bnd_list.append(this_ds)
                if self.cache:
                    self.write_to_cache(this_ds, url)
                this_ds = None

        msg.info("Merging dataset together (this might take a while)...")
        bnd = xr.concat(bnd_list, dim="time").squeeze('y')

        time = bnd.time.values
        freq = bnd.freq.values
        dirs = bnd.direction.values
        spec = bnd.SPEC.values
        lon = bnd.longitude.values
        lat = bnd.latitude.values

        source = f"{bnd.title}, {bnd.institution}"

        return  time, freq, dirs, spec, lon, lat, source

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
        url = 'https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/'+day.strftime('%d')+'/MyWave_wam4_SPC_'+day.strftime('%Y%m%d')+'T'+day.strftime('%H')+'Z.nc'
        return url

    def get_filepath_if_cached(self, url: str) -> Tuple:
        """
        Returns the filepath if the file is cached locally, otherwise
        hands back the URL.
        """
        maybe_cache = self._url_to_filename(url)
        if os.path.exists(maybe_cache):
            return maybe_cache, True
        else:
            return url, False

    def write_to_cache(self, ds, url):
        cache = self._url_to_filename(url)
        msg.plain(f'Caching {url} to {cache}')
        ds.to_netcdf(cache)

    def _url_to_filename(self, url):
        """
        Sanitizes a url to a valid file name.
        """
        fname = "".join(x for x in url if x.isalnum() or x == '.')
        return os.path.join(self.cache_folder, fname)


class NORA3(BoundaryReader):
    def __init__(self, stride: int=24, hours_per_file: int=24,
                last_file: str='', lead_time: int=0,
                source: str='thredds',
                folder: str='input/NORA3') -> None:
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)
        self.source = source
        self.folder = folder
        return

    def convention(self) -> str:
        return 'Ocean'

    def get_coordinates(self, start_time) -> Tuple:
        """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""
        #day = pd.date_range(start_time, start_time,freq='D')
        start_times, end_times, file_times = create_time_stamps(start_time, start_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)
        url = self.get_url(file_times[0], source=self.source)

        data = xr.open_dataset(url).isel(time = [0])

        lon_all = data.longitude.values[0]
        lat_all = data.latitude.values[0]

        return lon_all, lat_all

    def __call__(self, start_time, end_time, inds) -> Tuple:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(start_time, end_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)
        #days = day_list(start_time = self.start_time, end_time = self.end_time)
        msg.info(f"Getting boundary spectra from NORA3 from {self.start_time} to {self.end_time}")
        bnd_list = []
        for n in range(len(file_times)):
            url = self.get_url(file_times[n], source=self.source)
            msg.from_file(url)
            msg.plain(f"Reading boundary spectra: {start_times[n]}-{end_times[n]}")
            with xr.open_dataset(url) as f:
                this_ds = f.sel(time = slice(start_times[n], end_times[n]), x = (inds+1))[['SPEC', 'longitude', 'latitude', 'time', 'freq', 'direction']].copy()
            bnd_list.append(this_ds)
            #bnd_list.append(xr.open_dataset(url).sel(time = slice(start_times[n], end_times[n]), x = (inds+1)))
        bnd=xr.concat(bnd_list, dim="time").squeeze('y')

        time = bnd.time.values
        freq = bnd.freq.values
        dirs = bnd.direction.values
        spec = bnd.SPEC.values
        lon = bnd.longitude.values[0,:]
        lat = bnd.latitude.values[0,:]

        source = f"{bnd.title}, {bnd.institution}"

        return  time, freq, dirs, spec, lon, lat, source


    def get_url(self, day, source) -> str:
        if source == 'thredds':
            return 'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/SPC'+day.strftime('%Y%m%d')+'00.nc'
        if source == 'lustre':
            return '/lustre/storeB/project/fou/om/WINDSURFER/mw3hindcast/spectra/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/SPC'+day.strftime('%Y%m%d')+'00.nc'
        if source == 'local':
            return self.folder + '/SPC'+day.strftime('%Y%m%d')+'00.nc'


class WW3_4km(BoundaryReader):
    def __init__(self, ignore_nan: bool=True, stride: int=6,
                 hours_per_file: int=73, last_file: str='', lead_time: int=0,
                 cache: bool=True, clean_cache: bool=False, tile='POI') -> None:
        self.ignore_nan = copy(ignore_nan)
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)
        self.cache = copy(cache)
        self.clean_cache = copy(clean_cache)
        self.cache_folder = 'dnora_bnd_temp'
        self.tile = tile # 'POI for all the domain', 'C0', ... 'C5' for WAM800 coastal areas
        return

    def convention(self) -> str:
        return 'Ocean'


    def get_coordinates(self, start_time) -> Tuple:
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

        return lon_all, lat_all

    def __call__(self, start_time, end_time, inds) -> Tuple:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        # If we have removed some spectra (NaN's) we need to remap the indeces
        if hasattr(self, 'pointers'):
            inds = self.pointers[0][inds]

        if self.clean_cache:
            for f in glob.glob(os.path.join(self.cache_folder, '*')):
                os.remove(f)

        if self.cache:
            if not os.path.isdir(self.cache_folder):
                os.mkdir(self.cache_folder)
                print("Creating cache folder %s..." % self.cache_folder)

        start_times, end_times, file_times = create_time_stamps(start_time, end_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)


        msg.info(f"Getting boundary spectra from WW3_4km from {self.start_time} to {self.end_time}")
        bnd_list = []
        for n in range(len(file_times)):


            msg.plain(f"Reading boundary spectra: {start_times[n]}-{end_times[n]}")

            file_time = file_times[n]
            ct = 1
            keep_trying = True
            while keep_trying:

                url = self.get_url(file_time)
                url_or_cache, data_from_cache = self.get_filepath_if_cached(url)
                if data_from_cache and ct == 1:
                    this_ds = xr.load_dataset(url_or_cache)
                else:
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
                msg.from_file(url_or_cache)
                bnd_list.append(this_ds)
                if self.cache:
                    self.write_to_cache(this_ds, url)
                this_ds = None

        msg.info("Merging dataset together (this might take a while)...")
        bnd = xr.concat(bnd_list, dim="time").squeeze('y')

        time = bnd.time.values
        freq = bnd.freq.values
        dirs = bnd.direction.values
        spec = bnd.SPEC.values
        lon = bnd.longitude.values
        lat = bnd.latitude.values

        source = f"{bnd.title}, {'Norwegian Meterological Institute'}"

        return  time, freq, dirs, spec, lon, lat, source

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
        url = 'https://thredds.met.no/thredds/dodsC/ww3_4km_archive_files/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/'+day.strftime('%d')+'/ww3_4km_'+self.tile+'_SPC_'+day.strftime('%Y%m%d')+'T'+day.strftime('%H')+'Z.nc'
        return url

    def get_filepath_if_cached(self, url: str) -> Tuple:
        """
        Returns the filepath if the file is cached locally, otherwise
        hands back the URL.
        """
        maybe_cache = self._url_to_filename(url)
        if os.path.exists(maybe_cache):
            return maybe_cache, True
        else:
            return url, False

    def write_to_cache(self, ds, url):
        cache = self._url_to_filename(url)
        msg.plain(f'Caching {url} to {cache}')
        ds.to_netcdf(cache)

    def _url_to_filename(self, url):
        """
        Sanitizes a url to a valid file name.
        """
        fname = "".join(x for x in url if x.isalnum() or x == '.')
        return os.path.join(self.cache_folder, fname)
