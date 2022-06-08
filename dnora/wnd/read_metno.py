from abc import ABC,  abstractmethod
from copy import copy
import numpy as np
import xarray as xr
from subprocess import call
import os, glob

# Import objects
from ..grd.grd_mod import Grid

# Import abstract classes
from .read import ForcingReader

# Import auxiliry functions
from .. import msg
from ..aux import create_time_stamps, u_v_from_dir, expand_area, lon_in_km


class NORA3(ForcingReader):
    """Reads wind data of the NORA3 hindcast directly from MET Norways servers.

    The NORA3 HARMONIE-AROME high-resolution (ca 3 km) hindcast for the
    North Sea, the Norwegian Sea, and the Barents Sea.

    Haakenstad, H., Breivik, Ø., Furevik, B. R., Reistad, M., Bohlinger, P., &
    Aarnes, O. J. (2021). NORA3: A Nonhydrostatic High-Resolution Hindcast of
    the North Sea, the Norwegian Sea, and the Barents Sea,
    Journal of Applied Meteorology and Climatology, 60(10), 1443-1464,
    DOI: 10.1175/JAMC-D-21-0029.1
    """

    def __init__(self, stride: int=1, hours_per_file: int=1, last_file: str='', lead_time: int=4):
        """The data is currently in hourly files. Do not change the default
        setting unless you have a good reason to do so.
        """

        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)
        return

    def __call__(self, grid: Grid, start_time: str, end_time: str, expansion_factor: float):
        """Reads boundary spectra between given times and given area around
        the Grid object."""

        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(
            start_time, end_time, self.stride, self.hours_per_file, self.last_file, self.lead_time)

        msg.info(
            f"Getting wind forcing from NORA3 from {self.start_time} to {self.end_time}")


        temp_folder = 'dnora_wnd_temp'
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        msg.plain("Removing old files from temporary folder...")
        for f in glob.glob("dnora_wnd_temp/*MetNo_NORA3.nc"):
            os.remove(f)

        # Define area to search in
        lon_min, lon_max, lat_min, lat_max = expand_area(min(grid.lon()), max(grid.lon()), min(grid.lat()), max(grid.lat()), expansion_factor)

        # Set resolution to about 3 km
        dlat = 3/111
        mean_lon_in_km = (lon_in_km(grid.lat()[0])+lon_in_km(grid.lat()[-1]))*0.5
        dlon = 3/mean_lon_in_km

        urls = [
            self.get_url(file_time, start_time, first_ind=self.lead_time)
            for file_time, start_time in zip(file_times, start_times)
        ]

        msg.plain('Reading from these URLs:')
        for url in urls:
            msg.plain(url)

        ds = xr.open_mfdataset(urls)
        ds0 = xr.open_dataset(urls[0])
        # subsets
        c = (
            (ds.longitude>=lon_min-dlon) & (ds.longitude<=lon_max+dlon)
            & (ds.latitude>=lat_min-dlat) & (ds.latitude<=lat_max+dlat)
        ).compute()
        cx = c.any('y')
        cy = c.any('x')
        ct = ((ds.time>=start_times[0]) & (ds.time<=end_times[-1])).compute()
        ds = ds.isel(x=cx, y=cy, time=ct)[[
            'wind_direction', 'wind_speed',
            'projection_lambert', 'forecast_reference_time',
        ]]

        nc_fimex_in = 'dnora_fimex_in.nc'
        nc_fimex_out = 'dnora_fimex_out.nc'
        ds.to_netcdf(nc_fimex_in)

        fimex_command = ['fimex', '--input.file='+nc_fimex_in,
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
                         # '--extract.reduceTime.start=' + \
                         # start_times[0].strftime('%Y-%m-%dT%H:%M:%S'),
                         # '--extract.reduceTime.end=' + \
                         # end_times[-1].strftime('%Y-%m-%dT%H:%M:%S'),
                         '--process.rotateVector.direction=latlon',
                         '--debug',
                         '--output.file='+nc_fimex_out]

        call(fimex_command)

        wind_forcing = xr.open_dataset(nc_fimex_out).squeeze()

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
        wind_forcing['u'] = u #(['time', 'lat', 'lon'],  u)
        wind_forcing['v'] = v # (['time', 'lat', 'lon'],  v)

        return wind_forcing

    def get_url(self, time_stamp_file, time_stamp, first_ind) -> str:
        h0 = int(time_stamp_file.hour) % 6
        folder = time_stamp_file.strftime('%Y')+'/'+time_stamp_file.strftime('%m')+'/'+time_stamp_file.strftime('%d')+'/'+(time_stamp_file - np.timedelta64(h0, 'h')).strftime('%H')
        ind = int((time_stamp.hour-first_ind) % 6) + first_ind
        filename = 'fc' + time_stamp_file.strftime('%Y')+time_stamp_file.strftime('%m')+time_stamp_file.strftime('%d')+(time_stamp_file - np.timedelta64(h0, 'h')).strftime('%H')+'_' + f"{ind:03d}" + '_fp.nc'
        url = 'https://thredds.met.no/thredds/dodsC/nora3/'+folder + '/' + filename
        return url


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
        lon_min, lon_max, lat_min, lat_max = expand_area(min(grid.lon()), max(grid.lon()), min(grid.lat()), max(grid.lat()), expansion_factor)

        # Setting resolution to roughly 3 km
        dlat = 3/111
        mean_lon_in_km = (lon_in_km(grid.lat()[0])+lon_in_km(grid.lat()[-1]))*0.5
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
            print('No')
            prefix = 'subset'

        # Set resolution to about 2.5 km
        dlat = 2.5/111
        mean_lon_in_km = (lon_in_km(grid.lat()[0])+lon_in_km(grid.lat()[-1]))*0.5
        dlon = 2.5/mean_lon_in_km

        # Define area to search in
        lon_min, lon_max, lat_min, lat_max = expand_area(min(grid.lon()), max(grid.lon()), min(grid.lat()), max(grid.lat()), expansion_factor)

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

            if prefix == 'subset':
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
