from abc import ABC,  abstractmethod
from copy import copy
import numpy as np
import xarray as xr
from subprocess import call
from .. import msg
from ..aux import create_time_stamps, u_v_from_dir, expand_area, lon_in_km

import os

from .wnd_mod import ForcingReader # Abstract class
from .wnd_mod import Forcing # Forcing object

from ..grd.grd_mod import Grid # Grid object

class MetNo_NORA3(ForcingReader):
    def __init__(self, stride: int = 1, hours_per_file: int = 1, last_file: str = '', lead_time: int = 4):
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

        wnd_list = []

        temp_folder = 'dnora_wnd_temp'
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        msg.info(
            f"Getting wind forcing from Arome 2.5 km from {self.start_time} to {self.end_time}")
        for n in range(len(file_times)):

            url = self.get_url(file_times[n], start_times[n], first_ind=self.lead_time)

            msg.info(url)
            msg.plain(f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")

            nc_fimex = f'dnora_wnd_temp/wind_{n:04.0f}.nc'
            #nc_fimex = 'dnora_wnd_temp.nc'

            # Define area to search in
            lon_min, lon_max, lat_min, lat_max = expand_area(grid.lon()[0], grid.lon()[-1], grid.lat()[0], grid.lat()[-1], expansion_factor)

            # Set resolution to about 3 km
            dlat = 3/111
            mean_lon_in_km = (lon_in_km(grid.lat[0])+lon_in_km(grid.lat[-1]))*0.5
            dlon = 3/mean_lon_in_km

            fimex_command = ['fimex', '--input.file='+url,
                             '--interpolate.method=bilinear',
                             '--interpolate.projString=+proj=latlong +ellps=sphere +a=6371000 +e=0',
                             '--interpolate.xAxisValues='
                             + str(lon_min)+','+str(lon_min+dlon)
                             + ',...,'+str(lon_max)+'',
                             '--interpolate.yAxisValues='
                             + str(lat_min)+','+str(lat_min+dlat)
                             + ',...,'+str(lat_max)+'',
                             '--interpolate.xAxisUnit=degree', '--interpolate.yAxisUnit=degree',
                             '--process.rotateVector.all',
                             '--extract.selectVariables=wind_speed', '--extract.selectVariables=wind_direction',
                             '--extract.reduceTime.start=' + \
                             start_times[n].strftime('%Y-%m-%dT%H:%M:%S'),
                             '--extract.reduceTime.end=' + \
                             end_times[n].strftime('%Y-%m-%dT%H:%M:%S'),
                             '--process.rotateVector.direction=latlon',
                             '--output.file='+nc_fimex]


            call(fimex_command)
            wnd_list.append(xr.open_dataset(nc_fimex).squeeze())

        wind_forcing = xr.concat(wnd_list, dim="time")

        # Fimex has already rotated the longitudes and latitudes, so calling them rlon/rlat is now incorrect
        wind_forcing = wind_forcing.rename_dims({'y': 'lat', 'x': 'lon'})
        wind_forcing = wind_forcing.rename_vars({'y': 'lat', 'x': 'lon'})

        # Go to u and v components
        u, v = u_v_from_dir(wind_forcing.wind_speed,
                            wind_forcing.wind_direction)
        u = u.fillna(0)
        v = v.fillna(0)

        # Remove speed and dir and add components to dataset
        wind_forcing = wind_forcing.drop_vars(['wind_speed', 'wind_direction'])
        wind_forcing["u"] = (['time', 'lat', 'lon'],  u)
        wind_forcing["v"] = (['time', 'lat', 'lon'],  v)
        return wind_forcing

    def get_url(self, time_stamp_file, time_stamp, first_ind=4):
        h0 = int(time_stamp_file.hour) % 6
        folder = time_stamp_file.strftime('%Y')+'/'+time_stamp_file.strftime('%m')+'/'+time_stamp_file.strftime('%d')+'/'+(time_stamp_file - np.timedelta64(h0, 'h')).strftime('%H')
        ind = int((time_stamp.hour-first_ind) % 6) + first_ind
        filename = 'fc' + time_stamp_file.strftime('%Y')+time_stamp_file.strftime('%m')+time_stamp_file.strftime('%d')+(time_stamp_file - np.timedelta64(h0, 'h')).strftime('%H')+'_' + f"{ind:03d}" + '_fp.nc'
        url = 'https://thredds.met.no/thredds/dodsC/nora3/'+folder + '/' + filename
        return url


class MetNo_MyWave3km(ForcingReader):
    def __init__(self, stride: int = 24, hours_per_file: int = 24, last_file: str = '', lead_time: int = 0):
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

        wnd_list = []

        temp_folder = 'dnora_wnd_temp'
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        msg.info(
            f"Getting wind forcing from MEPS from {self.start_time} to {self.end_time}")
        for n in range(len(file_times)):
            url = self.get_url(file_times[n])

            msg.info(url)
            msg.plain(
                f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")

            nc_fimex = f'dnora_wnd_temp/wind_{n:04.0f}.nc'

            # Define area to search in
            lon_min, lon_max, lat_min, lat_max = expand_area(grid.lon()[0], grid.lon()[-1], grid.lat()[0], grid.lat()[-1], expansion_factor)

            dlat = 3/111
            mean_lon_in_km = (lon_in_km(grid.lat[0])+lon_in_km(grid.lat[-1]))*0.5
            dlon = 3/mean_lon_in_km

            fimex_command = ['fimex', '--input.file='+url,
                             '--interpolate.method=bilinear',
                             '--interpolate.projString=+proj=latlong +ellps=sphere +a=6371000 +e=0',
                             '--interpolate.xAxisValues='+ str(lon_min)+','+str(lon_min+dlon)+ ',...,'+str(lon_max)+'',
                             '--interpolate.yAxisValues='           + str(lat_min)+','+str(lat_min+dlat)+ ',...,'+str(lat_max)+'',
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

        # Fimex has already rotated the longitudes and latitudes, so calling them rlon/rlat is now incorrect
        wind_forcing = wind_forcing.rename_dims({'rlat': 'lat', 'rlon': 'lon'})
        wind_forcing = wind_forcing.rename_vars({'rlat': 'lat', 'rlon': 'lon'})

        # Go to u and v components
        u, v = u_v_from_dir(wind_forcing.ff, wind_forcing.dd)  # factor 1000
        u = u.fillna(0)
        v = v.fillna(0)

        # Remove speed and dir and add components to dataset
        wind_forcing = wind_forcing.drop_vars(['ff', 'dd'])
        wind_forcing["u"] = (['time', 'lat', 'lon'],  u)
        wind_forcing["v"] = (['time', 'lat', 'lon'],  v)
        return wind_forcing

    def get_url(self, time_stamp):
        filename = time_stamp.strftime(
            '%Y')+time_stamp.strftime('%m')+time_stamp.strftime('%d')+'_MyWam3km_hindcast.nc'
        url = 'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/' + \
            time_stamp.strftime('%Y')+'/' + \
            time_stamp.strftime('%m')+'/'+filename
        return url


class MetNo_MEPS(ForcingReader):
    def __init__(self, prefix: str = 'subset', stride: int = 24, hours_per_file: int = 24, last_file: str = '', lead_time: int = 0):
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.prefix = copy(prefix)
        self.last_file = copy(last_file)
        return

    def __call__(self, grid: Grid, start_time: str, end_time: str, expansion_factor: float):
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(
            start_time, end_time, self.stride, self.hours_per_file, self.last_file, self.lead_time)

        wnd_list = []

        temp_folder = 'dnora_wnd_temp'
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        msg.info(
            f"Getting wind forcing from MEPS from {self.start_time} to {self.end_time}")
        for n in range(len(file_times)):
            url = self.get_url(file_times[n], self.prefix)

            msg.info(url)
            msg.plain(
                f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")

            nc_fimex = f'dnora_wnd_temp/wind_{n:04.0f}.nc'

            # Define area to search in
            lon_min, lon_max, lat_min, lat_max = expand_area(grid.lon()[0], grid.lon()[-1], grid.lat()[0], grid.lat()[-1], expansion_factor)

            # Set resolution to about 2.5 km
            dlat = 2.5/111
            mean_lon_in_km = (lon_in_km(grid.lat[0])+lon_in_km(grid.lat[-1]))*0.5
            dlon = 2.5/mean_lon_in_km

            fimex_command = ['fimex', '--input.file='+url,
                             '--interpolate.method=bilinear',
                             '--interpolate.projString=+proj=latlong +ellps=sphere +a=6371000 +e=0',
                             '--interpolate.xAxisValues='
                             + str(lon_min)+','+str(lon_min+dlon)
                             + ',...,'+str(lon_max)+'',
                             '--interpolate.yAxisValues='
                             + str(lat_min)+','+str(lat_min+dlat)
                             + ',...,'+str(lat_max)+'',
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

            if self.prefix == 'subset':
                fimex_command.insert(-2,
                                     '--extract.reduceDimension.name=ensemble_member')
                fimex_command.insert(-2, '--extract.reduceDimension.start=1')
                fimex_command.insert(-2, '--extract.reduceDimension.end=1')

            call(fimex_command)
            wnd_list.append(xr.open_dataset(nc_fimex).squeeze())

        wind_forcing = xr.concat(wnd_list, dim="time")

        # Fimex has already rotated the longitudes and latitudes, so calling them rlon/rlat is now incorrect
        wind_forcing = wind_forcing.rename_dims({'y': 'lat', 'x': 'lon'})
        wind_forcing = wind_forcing.rename_vars({'y': 'lat', 'x': 'lon'})

        wind_forcing = wind_forcing.rename_vars(
            {'x_wind_10m': 'u', 'y_wind_10m': 'v'})

        #wind_forcing.to_netcdf('test.nc')
        return wind_forcing

    def get_url(self, time_stamp, prefix):
        filename = 'meps_'+prefix+'_2_5km_'+time_stamp.strftime('%Y')+time_stamp.strftime('%m')+time_stamp.strftime('%d')+'T'+time_stamp.strftime('%H')+'Z.nc'
        url = 'https://thredds.met.no/thredds/dodsC/meps25epsarchive/'+time_stamp.strftime('%Y')+'/'+time_stamp.strftime('%m')+'/'+time_stamp.strftime('%d')+'/' + filename
        return url
