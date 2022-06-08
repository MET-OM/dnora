import os
import xarray as xr
import glob
import numpy as np
from copy import copy
from abc import ABC, abstractmethod
from typing import Tuple

# Import abstract classes and needed instances of them
from .read import BoundaryReader

# Import auxiliry functions
from .. import msg
from ..aux import day_list, create_time_stamps, CachedReaderMixin

class CachedBoundaryReader(CachedReaderMixin):
    """
    Not usable on its own.
    Child class needs to supply some methods.
    """
    def __init__(self, stride: int,
                 hours_per_file: int, last_file: str, lead_time: int, *argv,
                 clean_cache: bool=False,  **kwarg) -> None:
        super().__init__(*argv, **kwarg)
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)
        self.cache = copy(cache)
        self.clean_cache = copy(clean_cache)
        self.cache_folder = 'dnora_bnd_temp'
        return

    def __call__(self, start_time, end_time, inds) -> Tuple:
        """Reads in all boundary spectra between the given times and at for the given indices"""
        self.start_time = start_time
        self.end_time = end_time

        # If we have removed some spectra (NaN's) we need to remap the indices
        if hasattr(self, 'pointers'):
            inds = self.pointers[0][inds]

        if self.clean_cache:
            for f in glob.glob(os.path.join(self.cache_folder, '*')):
                os.remove(f)

        if self.cache:
            if not os.path.isdir(self.cache_folder):
                os.mkdir(self.cache_folder)
                print("Creating cache folder %s..." % cache_folder)

        start_times, end_times, file_times = create_time_stamps(start_time, end_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)

        msg.info(f"Getting boundary spectra from {self.backend_name} from {self.start_time} to {self.end_time}")
        bnd_list = []
        urls = [self.get_url(file_time) for file_time in file_times]
        for file_time, start_time, end_time, url in zip(file_times, start_times, end_times, urls):
            url_or_cache, data_from_cache = self.get_filepath_if_cached(url)
            msg.from_file(url_or_cache)
            msg.plain(f"Reading boundary spectra: {start_time}-{end_time}")
            if data_from_cache:
                this_ds = xr.load_dataset(url_or_cache)
                bnd_list.append(this_ds)
            else:
                try:
                    this_ds = self.load_data_from_url(url, start_time, end_time, inds)
                    if self.is_consistent(this_ds, bnd_list):
                        bnd_list.append(this_ds)
                        cache = self._url_to_filename(url)
                        self.write_to_cache(this_ds, url, cache)
                    else:
                        msg.plain(f'SKIPPING, file inconsistent: {url}')
                except OSError:
                    msg.plain(f'SKIPPING, file not found: {url}')

        msg.info("Merging dataset together (this might take a while)...")
        ds = xr.concat(bnd_list, dim='time')
        return self.extract_data(ds)


class CachedMetnoBoundaryReader(CachedBoundaryReader):
    def load_data_from_url(self, url, start_time, end_time, inds):
        with xr.open_dataset(url) as f:
            this_ds = f.sel(
                time=slice(start_time, end_time),
                x=inds+1,
            )[
                ['SPEC', 'longitude', 'latitude', 'time', 'freq', 'direction']
            ].copy()
        return this_ds

    def convention(self) -> str:
        return 'Ocean'

    def get_coordinates(self, time) -> Tuple:
        """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""
        start_times, end_times, file_times = create_time_stamps(time, time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)
        url = self.get_url(file_times[0])

        data = xr.open_dataset(url).isel(time=[0])

        lon_all = data.longitude.values[0]
        lat_all = data.latitude.values[0]
        return lon_all, lat_all

class WAM4km(CachedMetnoBoundaryReader, CachedReaderMixin):
    """Read WAM4km boundary spectra"""
    def __init__(self, ignore_nan=True, stride=6,
                 hours_per_file=73, last_file='', lead_time=0,
                 **kwarg) -> None:
        super().__init__(
            *argv,
            stride=stride,
            hours_per_file=hours_per_file,
            last_file=last_file,
            lead_time=lead_time,
            backend_name='WAM4km'
            **kwarg,
        )
        self.ignore_nan = copy(ignore_nan)
        return

    def get_coordinates(self, time) -> Tuple:
        lon_all, lat_all = super().get_coordinates(time)
        if self.ignore_nan:
            msg.info('ignore_nan = True. The following points are NOT provided to the point_picker')
            mask = [True]*len(lon_all)
            for n in range(len(lon_all)):
                if np.isnan(data.SPEC.values[:,0,n,:,:]).any():
                    msg.plain(f"{lon_all[n]:10.7f}, {lat_all[n]:10.7f}")
                    mask[n] = False

            lon_all = lon_all[mask]
            lat_all = lat_all[mask]
            self.pointers = np.where(mask) # These are used to map back the indices from the PointPicker to indices in the original data

        return lon_all, lat_all

    def get_url(self, day):
        url = 'https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/'+day.strftime('%d')+'/MyWave_wam4_SPC_'+day.strftime('%Y%m%d')+'T'+day.strftime('%H')+'Z.nc'
        return url

    @staticmethod
    def extract_data(ds):
        ds = ds.squeeze('y')

        time = ds.time.values
        freq = ds.freq.values
        dirs = ds.direction.values
        spec = ds.SPEC.values
        lon = ds.longitude.values
        lat = ds.latitude.values

        source = f"{ds.title}, {ds.institution}"

        return time, freq, dirs, spec, lon, lat, source


class NORA3(CachedMetnoBoundaryReader, CachedReaderMixin):
    """Read NORA3 boundary spectra"""
    def __init__(self, *argv,
                 stride: int=24, hours_per_file: int=24, last_file: str='',
                 lead_time: int=0,
                 **kwarg) -> None:
        super().__init__(
            *argv,
            stride=stride,
            hours_per_file=hours_per_file,
            last_file=last_file,
            lead_time=lead_time,
            backend_name='NORA3',
            **kwarg,
        )
        return

    def get_url(self, day) -> str:
        url = 'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/SPC'+day.strftime('%Y%m%d')+'00.nc'
        return url

    @staticmethod
    def extract_data(ds):
        ds = ds.squeeze('y')

        time = ds.time.values
        freq = ds.freq.values
        dirs = ds.direction.values
        spec = ds.SPEC.values
        lon = ds.longitude.values[0,:]
        lat = ds.latitude.values[0,:]

        source = f"{ds.title}, {ds.institution}"

        return time, freq, dirs, spec, lon, lat, source
