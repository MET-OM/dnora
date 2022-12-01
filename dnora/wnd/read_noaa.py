from abc import ABC,  abstractmethod
from copy import copy
import numpy as np
import xarray as xr
from subprocess import call
import os, glob
import time
# Import objects
from ..grd.grd_mod import Grid

# Import abstract classes
from .read import ForcingReader

# Import aux_funcsiliry functions
from .. import msg
from ..aux_funcs import create_time_stamps, u_v_from_dir, expand_area, lon_in_km


class GFS(ForcingReader):
    """Reads wind data of the GFS global forecast
    """

    def __init__(self, stride: int=6, hours_per_file: int=121, last_file: str='',
                lead_time: int=0):
        """The data is currently in hourly files. Do not change the default
        setting unless you have a good reason to do so.
        """

        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)


    def __call__(self, grid: Grid, start_time: str, end_time: str, expansion_factor: float):
        """Reads wind data between given times and given area around
        the Grid object."""

        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(
            start_time, end_time, self.stride, self.hours_per_file, self.last_file, self.lead_time)

        msg.info(
            f"Getting wind forcing from GFS from {self.start_time} to {self.end_time}")

        # Define area to search in
        lon, lat = expand_area(grid.edges('lon'), grid.edges('lat'), expansion_factor)

        wnd_list = []
        for n in range(len(file_times)):
            url = self.get_url(file_times[n], start_times[n], first_ind=self.lead_time)

            msg.from_file(url)
            msg.plain(f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")

            import dask
            with dask.config.set(**{'array.slicing.split_large_chunks': True}):
                with xr.open_dataset(url) as ds:
                    ds['time'] = ds.time.dt.round('H')
                    ds = ds.sel(time = slice(start_times[n], end_times[n]), lon=slice(lon[0], lon[1]), lat=slice(lat[0], lat[1]))[['lon', 'lat', 'time', 'ugrd10m', 'vgrd10m']]

            wnd_list.append(ds)

        wind_forcing = xr.concat(wnd_list, dim="time")

        u = wind_forcing.ugrd10m.values
        v = wind_forcing.vgrd10m.values
        u = np.moveaxis(u,0,2)
        v = np.moveaxis(v,0,2)
        # u = u.fillna(0)
        # v = v.fillna(0)

        time = wind_forcing.time.values
        lon = wind_forcing.lon.values
        lat = wind_forcing.lat.values
        x = None
        y = None

        return time, u, v, lon, lat, x, y, wind_forcing.attrs

    def get_url(self, time_stamp_file, time_stamp, first_ind) -> str:
        h0 = int(time_stamp_file.hour)
        folder = 'gfs' + time_stamp_file.strftime('%Y%m%d')
        filename = f'gfs_0p25_1hr_{h0:02.0f}z'

        return 'http://nomads.ncep.noaa.gov:80/dods/gfs_0p25_1hr/'+folder + '/' + filename


class MyWave3km(ForcingReader):
    """Reads wind data from the MyWave 3km hindcast directly from MET Norways
    servers. You should probably use MetNo_NORA3 because:

    The wind data is from NORA3 (see the MetNo_NORA3 reader), is taken
    from the wave model output. This means that model land points have no data.
    """

    def __init__(self, stride: int=24, hours_per_file: int=24, last_file: str='', lead_time: int=0):
        """The data is currently in daily files. Do not change the default
        setting unless you have a good reason to do so.
        """

        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)
        return

    def __call__(self, grid: Grid, start_time: str, end_time: str, expansion_factor: float):
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(
            start_time, end_time, self.stride, self.hours_per_file, self.last_file, self.lead_time)

        msg.info(
            f"Getting wind forcing from MEPS from {self.start_time} to {self.end_time}")

        temp_folder = 'dnora_wnd_temp'
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        msg.plain("Removing old files from temporary folder...")
        for f in glob.glob("dnora_wnd_temp/*MetNo_MyWave3km.nc"):
            os.remove(f)


        # Define area to search in
        lon, lat = expand_area(grid.edges('lon'), grid.edges('lat'), expansion_factor)

        # Setting resolution to roughly 3 km
        dlat = 3/111
        mean_lon_in_km = (lon_in_km(grid.edges('lat')[0])+lon_in_km(grid.edges('lat')[-1]))*0.5
        dlon = 3/mean_lon_in_km

        wnd_list = []
        for n in range(len(file_times)):
            url = self.get_url(file_times[n])

            msg.from_file(url)
            msg.plain(
                f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")

            nc_fimex = f'dnora_wnd_temp/wind_{n:04.0f}_MetNo_MyWave3km.nc'

            fimex_command = ['fimex', '--input.file='+url,
                             '--interpolate.method=bilinear',
                             '--interpolate.projString=+proj=latlong +ellps=sphere +a=6371000 +e=0',
                             '--interpolate.xAxisValues='+ str(lon[0])+','+str(lon[0]+dlon)+ ',...,'+str(lon[1])+'',
                             '--interpolate.yAxisValues='           + str(lat[0])+','+str(lat[0]+dlat)+ ',...,'+str(lat[1])+'',
                             '--interpolate.xAxisUnit=degree', '--interpolate.yAxisUnit=degree',
                             '--process.rotateVector.all',
                             '--extract.selectVariables=ff', '--extract.selectVariables=dd',
                             '--extract.reduceTime.start=' + \
                             start_times[n].strftime('%Y-%m-%dT%H:%M:%S'),
                             '--extract.reduceTime.end=' + \
                             end_times[n].strftime('%Y-%m-%dT%H:%M:%S'),
                             '--process.rotateVector.direction=latlon',
                             '--output.file='+nc_fimex]

            call(fimex_command)
            wnd_list.append(xr.open_dataset(nc_fimex).squeeze())

        wind_forcing = xr.concat(wnd_list, dim="time")


        # Go to u and v components
        u, v = u_v_from_dir(wind_forcing.ff, wind_forcing.dd)  # factor 1000

        u = u.fillna(0)
        v = v.fillna(0)
        u = np.moveaxis(u.values,0,2)
        v = np.moveaxis(v.values,0,2)

        time = wind_forcing.time.values
        lon = wind_forcing.rlon.values
        lat = wind_forcing.rlat.values
        x = None
        y = None

        return time, u, v, lon, lat, x, y, wind_forcing.attrs

    def get_url(self, time_stamp):
        filename = time_stamp.strftime('%Y')+time_stamp.strftime('%m')+time_stamp.strftime('%d')+'_MyWam3km_hindcast.nc'

        if self.source() == 'remote':
            return 'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/' + time_stamp.strftime('%Y')+'/' + time_stamp.strftime('%m')+'/'+filename
        else:
            return self.source() + '/'+filename


class MEPS(ForcingReader):
    """Reads wind data from MET Norways MEPS forecast.

    The data is from a 2.5 km AROME model.
    """

    def __init__(self, stride: int = 6, hours_per_file: int = 67, last_file: str = '', lead_time: int = 0):
        """The data is currently in 6 hourly files. Do not change the default
        setting unless you have a good reason to do so.
        """

        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)
        return

    def __call__(self, grid: Grid, start_time: str, end_time: str, expansion_factor: float):
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(
            start_time, end_time, self.stride, self.hours_per_file, self.last_file, self.lead_time)

        msg.info(
            f"Getting wind forcing from MEPS from {self.start_time} to {self.end_time}")

        temp_folder = 'dnora_wnd_temp'
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            msg.info("Creating folder %s..." % temp_folder)

        msg.plain("Removing old files from temporary folder...")
        for f in glob.glob("dnora_wnd_temp/*MetNo_MEPS.nc"):
            os.remove(f)

        # Check weather to use 'det' or 'subset' files
        url = self.get_url(file_times[0], 'det')
        try:
            xr.open_dataset(url)
            prefix = 'det'
        except:
            prefix = 'subset'

        # Set resolution to about 2.5 km
        dlat = 2.5/111
        mean_lon_in_km = (lon_in_km(grid.edges('lat')[0])+lon_in_km(grid.edges('lat')[-1]))*0.5
        dlon = 2.5/mean_lon_in_km

        # Define area to search in
        lon, lat = expand_area(grid.edges('lon'), grid.edges('lat'), expansion_factor)

        wnd_list = []
        for n in range(len(file_times)):
            msg.plain(
                f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")

            nc_fimex = f'dnora_wnd_temp/wind_{n:04.0f}_MetNo_MEPS.nc'
            url = self.get_url(file_times[n], prefix)
            msg.from_file(url)

            fimex_command = ['fimex', '--input.file='+url,
                             '--interpolate.method=bilinear',
                             '--interpolate.projString=+proj=latlong +ellps=sphere +a=6371000 +e=0',
                             '--interpolate.xAxisValues='
                             + str(lon[0])+','+str(lon[0]+dlon)
                             + ',...,'+str(lon[1])+'',
                             '--interpolate.yAxisValues='
                             + str(lat[0])+','+str(lat[0]+dlat)
                             + ',...,'+str(lat[1])+'',
                             '--interpolate.xAxisUnit=degree', '--interpolate.yAxisUnit=degree',
                             '--process.rotateVector.all',
                             '--extract.selectVariables=x_wind_10m', '--extract.selectVariables=y_wind_10m',
                             '--extract.selectVariables=latitude', '--extract.selectVariables=longitude',
                             '--extract.reduceTime.start=' + \
                             start_times[n].strftime('%Y-%m-%dT%H:%M:%S'),
                             '--extract.reduceTime.end=' + \
                             end_times[n].strftime('%Y-%m-%dT%H:%M:%S'),
                             '--process.rotateVector.direction=latlon',
                             '--output.file='+nc_fimex]

            if prefix == 'subset':
                fimex_command.insert(-2,
                                     '--extract.reduceDimension.name=ensemble_member')
                fimex_command.insert(-2, '--extract.reduceDimension.start=1')
                fimex_command.insert(-2, '--extract.reduceDimension.end=1')

            call(fimex_command)

            wnd_list.append(xr.open_dataset(nc_fimex).squeeze())

        wind_forcing = xr.concat(wnd_list, dim="time")

        u = wind_forcing.x_wind_10m.values
        v = wind_forcing.y_wind_10m.values
        u = np.moveaxis(u,0,2)
        v = np.moveaxis(v,0,2)

        time = wind_forcing.time.values
        lon = wind_forcing.x.values
        lat = wind_forcing.y.values
        x = None
        y = None

        return time, u, v, lon, lat, x, y, wind_forcing.attrs


    def get_url(self, time_stamp, prefix):
        filename = 'meps_'+prefix+'_2_5km_'+time_stamp.strftime('%Y')+time_stamp.strftime('%m')+time_stamp.strftime('%d')+'T'+time_stamp.strftime('%H')+'Z.nc'
        if self.source() == 'remote':
            return 'https://thredds.met.no/thredds/dodsC/meps25epsarchive/'+time_stamp.strftime('%Y')+'/'+time_stamp.strftime('%m')+'/'+time_stamp.strftime('%d')+'/' + filename
        else:
            return self.source() + '/' + filename
