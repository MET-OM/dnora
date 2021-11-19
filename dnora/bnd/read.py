import xarray as xr
import numpy as np
from ..aux import day_list, create_time_stamps, create_filename_time
from copy import copy
from .. import msg

from .bnd_mod import BoundaryReader # Abstract class
from .process import WW3ToOcean

class ForceFeed(BoundaryReader):
    def __init__(self, time, freq, dirs, spec, lon, lat):
        self.time = copy(time)
        self.freq = copy(freq)
        self.dirs = copy(dirs)
        self.spec = copy(spec)
        self.lon = copy(lon)
        self.lat = copy(lat)
        return

    def get_coordinates(self, start_time):
        return copy(self.lon), copy(self.lat)

    def __call__(self, start_time, end_time, inds):
        return  copy(self.time), copy(self).freq, copy(self).dirs, np.reshape(self.spec, (len(self.time), len(self.lon), self.spec.shape[0], self.spec.shape[1])), copy(self.lon), copy(self.lat), ''


class MetNo_WAM4km(BoundaryReader):
    def __init__(self, ignore_nan: bool=True, stride: int=6, hours_per_file: int=73, last_file: str='', lead_time: int=0):
        self.ignore_nan = copy(ignore_nan)
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)
        return

    def get_coordinates(self, start_time):
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

    def __call__(self, start_time, end_time, inds):
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
            url = self.get_url(file_times[n])
            msg.info(url)
            msg.plain(f"Reading boundary spectra: {start_times[n]}-{end_times[n]}")

            bnd_list.append(xr.open_dataset(url).sel(time = slice(start_times[n], end_times[n]), x = (inds+1)))

        msg.info("Merging dataset together (this might take a while)...")
        bnd=xr.concat(bnd_list, dim="time").squeeze('y')

        time = bnd.time.values
        freq = bnd.freq.values
        dirs = bnd.direction.values
        spec = bnd.SPEC.values
        lon = bnd.longitude.values
        lat = bnd.latitude.values

        source = f"{bnd.title}, {bnd.institution}"

        return  time, freq, dirs, spec, lon, lat, source

    def get_url(self, day):
        url = 'https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/'+day.strftime('%d')+'/MyWave_wam4_SPC_'+day.strftime('%Y%m%d')+'T'+day.strftime('%H')+'Z.nc'
        return url


class MetNo_NORA3(BoundaryReader):
    def __init__(self, stride: int=24, hours_per_file: int=24, last_file: str='', lead_time: int=0):
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)

    def get_coordinates(self, start_time):
        """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""
        #day = pd.date_range(start_time, start_time,freq='D')
        start_times, end_times, file_times = create_time_stamps(start_time, start_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)
        url = self.get_url(file_times[0])

        data = xr.open_dataset(url).isel(time = [0])

        lon_all = data.longitude.values[0]
        lat_all = data.latitude.values[0]

        return lon_all, lat_all

    def __call__(self, start_time, end_time, inds):
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(start_time, end_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)
        #days = day_list(start_time = self.start_time, end_time = self.end_time)

        msg.info(f"Getting boundary spectra from NORA3 from {self.start_time} to {self.end_time}")
        bnd_list = []
        for n in range(len(file_times)):
            url = self.get_url(file_times[n])
            msg.info(url)
            msg.plain(f"Reading boundary spectra: {start_times[n]}-{end_times[n]}")
            bnd_list.append(xr.open_dataset(url).sel(time = slice(start_times[n], end_times[n]), x = (inds+1)))

        bnd=xr.concat(bnd_list, dim="time").squeeze('y')

        time = bnd.time.values
        freq = bnd.freq.values
        dirs = bnd.direction.values
        spec = bnd.SPEC.values
        lon = bnd.longitude.values[0,:]
        lat = bnd.latitude.values[0,:]

        source = f"{bnd.title}, {bnd.institution}"

        return  time, freq, dirs, spec, lon, lat, source


    def get_url(self, day):
        url = 'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/SPC'+day.strftime('%Y%m%d')+'00.nc'
        return url

class File_WW3Nc(BoundaryReader):
    def __init__(self, folder: str='', filestring: str='ww3_T0', datestring: str='%Y%m%dT%H%M', stride: int=6, hours_per_file: int=73, last_file: str='', lead_time: int=0):
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)

        if (not folder == '') and (not folder[-1] == '/'):
            self.folder = folder + '/'
        else:
            self.folder = copy(folder)

        self.filestring = copy(filestring)
        self.datestring = copy(datestring)

    def get_coordinates(self, start_time):
        """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""
        #day = pd.date_range(start_time, start_time,freq='D')
        start_times, end_times, file_times = create_time_stamps(start_time, start_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)
        filename = self.get_filename(file_times[0])

        data = xr.open_dataset(filename).isel(time = [0])

        lon_all = data.longitude.values[0]
        lat_all = data.latitude.values[0]

        return lon_all, lat_all

    def __call__(self, start_time, end_time, inds):
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(start_time, end_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)

        msg.info(f"Getting boundary spectra from NORA3 from {self.start_time} to {self.end_time}")
        bnd_list = []
        for n in range(len(file_times)):
            filename = self.get_filename(file_times[n])
            msg.info(filename)
            msg.plain(f"Reading boundary spectra: {start_times[n]}-{end_times[n]}")

            bnd_list.append(xr.open_dataset(filename).sel(time = slice(start_times[n], end_times[n]), station = (inds+1)))

        bnd=xr.concat(bnd_list, dim="time")


        for x in range(len(inds)):
            for t in range(len(bnd.time.values)):
                new_spec, new_dirs = WW3ToOcean()(bnd.efth.values[t,x,:,:],bnd.direction.values)
                bnd.efth.values[t,x,:,:] = new_spec

        time = bnd.time.values
        freq = bnd.frequency.values
        dirs = new_dirs
        spec = bnd.efth.values
        lon = bnd.longitude.values[0,:]
        lat = bnd.latitude.values[0,:]

        source = f"ww3_ouput_spectra"

        return  time, freq, dirs, spec, lon, lat, source


    def get_filename(self, day):
        filename = self.folder + create_filename_time(self.filestring, [day], self.datestring) + '.nc'
        return filename
