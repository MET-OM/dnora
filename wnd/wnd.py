from abc import ABC,  abstractmethod
from dnora2 import msg
from dnora2.aux import create_time_stamps
from copy import copy
import numpy as np
import pandas as pd
from subprocess import call
import xarray as xr
import os


def u_v_from_dir(ws, wdir):
    # see http://tornado.sfsu.edu/geosciences/classes/m430/Wind/WindDirection.html
    u = -ws * (np.sin(np.deg2rad(wdir)))
    v = -ws * (np.cos(np.deg2rad(wdir)))
    return u, v

# =============================================================================
# FORCING FETCHER CLASSES RESPONSIBLE FOR ACTUALLY READING THE SPECTRA
# =============================================================================


class ForcingFetcher(ABC):
    def __init__(self):
        pass

# =============================================================================
#     @abstractmethod
#     def get_coordinates(self, start_time):
#         pass
# =============================================================================

    @abstractmethod
    def __call__(self, start_time, end_time, inds):
        pass

    def get_time_limits_day(self, ind):
        """Determines star and end time for the day. First and last day doesn't start at 00:00 or end at 23:59"""

        days = bnd.day_list(start_time=self.start_time, end_time=self.end_time)

        if ind == 0:
            t0 = self.start_time
            t1 = days[0].strftime('%Y-%m-%d') + "T23:59:59"
        elif ind == (len(days)-1):
            t0 = days[-1].strftime('%Y-%m-%d') + "T00:00:00"
            t1 = self.end_time
        else:
            t0 = days[ind].strftime('%Y-%m-%d') + "T00:00:00"
            t1 = days[ind].strftime('%Y-%m-%d') + "T23:59:59"
        return t0, t1

    def __str__(self):
        return (f"{self.start_time} - {self.end_time}")


class ForcingArome25thredds(ForcingFetcher):
    def __init__(self, stride=1, hours_per_file=1, last_file=None, lead_time=4):
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)
        return

    def __call__(self, grid, start_time, end_time, expansion_factor):
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        #days = bnd.day_list(start_time = self.start_time, end_time = self.end_time)
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
            #print(time_stamp)
            #print(days[n].strftime('%Y-%m-%d'))
            url = self.get_url(
                file_times[n], start_times[n], first_ind=self.lead_time)

            msg.info(url)
            msg.plain(
                f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")

            nc_fimex = f'dnora_wnd_temp/wind_{n:04.0f}.nc'
            #nc_fimex = 'dnora_wnd_temp.nc'

            # Define area to search in
            expand_lon = (grid.lon()[-1] - grid.lon()
                          [0])*(expansion_factor-1)*0.5
            expand_lat = (grid.lat()[-1] - grid.lat()
                          [0])*(expansion_factor-1)*0.5

            lon_min = grid.lon()[0] - expand_lon
            lon_max = grid.lon()[-1] + expand_lon

            lat_min = grid.lat()[0] - expand_lat
            lat_max = grid.lat()[-1] + expand_lat

            # Temporary hack: set resolution to about 2.5 km
            dlat = 2.5/111
            dlon = dlat*2

            # This doesn't work if we havent set a grid resolution, for example when only having the outer boundaries of an unstructured grid...
            #dlon = grid.data.dlon*5
            #dlat = grid.data.dlat*5

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
                             #'--extract.reduceToBoundingBox.south=' + str(lat_min),
                             #'--extract.reduceToBoundingBox.north=' + str(lat_max),
                             #'--extract.reduceToBoundingBox.west=' + str(lon_min),
                             #'--extract.reduceToBoundingBox.east=' + str(lon_max),
                             '--process.rotateVector.all',
                             '--extract.selectVariables=wind_speed', '--extract.selectVariables=wind_direction',
                             #'--extract.selectVariables=latitude','--extract.selectVariables=longitude',
                             '--extract.reduceTime.start=' + \
                             start_times[n].strftime('%Y-%m-%dT%H:%M:%S'),
                             '--extract.reduceTime.end=' + \
                             end_times[n].strftime('%Y-%m-%dT%H:%M:%S'),
                             '--process.rotateVector.direction=latlon',
                             '--output.file='+nc_fimex]

            #if self.prefix == 'subset':
            #    fimex_command.insert(-2,'--extract.reduceDimension.name=ensemble_member')
            #    fimex_command.insert(-2,'--extract.reduceDimension.start=1')
            #    fimex_command.insert(-2,'--extract.reduceDimension.end=1')

            call(fimex_command)
            wnd_list.append(xr.open_dataset(nc_fimex).squeeze())

        wind_forcing = xr.concat(wnd_list, dim="time")

        # Fimex has already rotated the longitudes and latitudes, so calling them rlon/rlat is now incorrect
        wind_forcing = wind_forcing.rename_dims({'y': 'lat', 'x': 'lon'})
        wind_forcing = wind_forcing.rename_vars({'y': 'lat', 'x': 'lon'})

        # Go to u and v components
        u, v = u_v_from_dir(wind_forcing.wind_speed,
                            wind_forcing.wind_direction)  # factor 1000
        u = u.fillna(0)
        v = v.fillna(0)

        # Remove speed and dir and add components to dataset
        wind_forcing = wind_forcing.drop_vars(['wind_speed', 'wind_direction'])
        wind_forcing["u"] = (['time', 'lat', 'lon'],  u)
        wind_forcing["v"] = (['time', 'lat', 'lon'],  v)
        #wind_forcing.to_netcdf('test.nc')
        return wind_forcing

    def get_url(self, time_stamp_file, time_stamp, first_ind=4):
        h0 = int(time_stamp_file.hour) % 6
        folder = time_stamp_file.strftime('%Y')+'/'+time_stamp_file.strftime('%m')+'/'+time_stamp_file.strftime(
            '%d')+'/'+(time_stamp_file - np.timedelta64(h0, 'h')).strftime('%H')
        ind = int((time_stamp.hour-first_ind) % 6) + first_ind
        filename = 'fc' + time_stamp_file.strftime('%Y')+time_stamp_file.strftime('%m')+time_stamp_file.strftime(
            '%d')+(time_stamp_file - np.timedelta64(h0, 'h')).strftime('%H')+'_' + f"{ind:03d}" + '_fp.nc'
        url = 'https://thredds.met.no/thredds/dodsC/nora3/'+folder + '/' + filename
        return url


class ForcingMyWave3km(ForcingFetcher):
    def __init__(self, stride=24, hours_per_file=24, last_file=None, lead_time=0):
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)
        return

    def __call__(self, grid, start_time, end_time, expansion_factor):
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        #days = bnd.day_list(start_time = self.start_time, end_time = self.end_time)
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
            #print(time_stamp)
            #print(days[n].strftime('%Y-%m-%d'))
            url = self.get_url(file_times[n])

            msg.info(url)
            msg.plain(
                f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")

            nc_fimex = f'dnora_wnd_temp/wind_{n:04.0f}.nc'
            #nc_fimex = 'dnora_wnd_temp.nc'

            # Define area to search in
            expand_lon = (grid.lon()[-1] - grid.lon()
                          [0])*(expansion_factor-1)*0.5
            expand_lat = (grid.lat()[-1] - grid.lat()
                          [0])*(expansion_factor-1)*0.5

            lon_min = grid.lon()[0] - expand_lon
            lon_max = grid.lon()[-1] + expand_lon

            lat_min = grid.lat()[0] - expand_lat
            lat_max = grid.lat()[-1] + expand_lat

            dlon = grid.data.dlon*5
            dlat = grid.data.dlat*5

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
                             #'--extract.reduceToBoundingBox.south=' + str(lat_min),
                             #'--extract.reduceToBoundingBox.north=' + str(lat_max),
                             #'--extract.reduceToBoundingBox.west=' + str(lon_min),
                             #'--extract.reduceToBoundingBox.east=' + str(lon_max),
                             '--process.rotateVector.all',
                             '--extract.selectVariables=ff', '--extract.selectVariables=dd',
                             #'--extract.selectVariables=latitude','--extract.selectVariables=longitude',
                             '--extract.reduceTime.start=' + \
                             start_times[n].strftime('%Y-%m-%dT%H:%M:%S'),
                             '--extract.reduceTime.end=' + \
                             end_times[n].strftime('%Y-%m-%dT%H:%M:%S'),
                             '--process.rotateVector.direction=latlon',
                             '--output.file='+nc_fimex]

            #if self.prefix == 'subset':
            #    fimex_command.insert(-2,'--extract.reduceDimension.name=ensemble_member')
            #    fimex_command.insert(-2,'--extract.reduceDimension.start=1')
            #    fimex_command.insert(-2,'--extract.reduceDimension.end=1')

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
        #wind_forcing.to_netcdf('test.nc')
        return wind_forcing

    def get_url(self, time_stamp):
        filename = time_stamp.strftime(
            '%Y')+time_stamp.strftime('%m')+time_stamp.strftime('%d')+'_MyWam3km_hindcast.nc'
        url = 'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/' + \
            time_stamp.strftime('%Y')+'/' + \
            time_stamp.strftime('%m')+'/'+filename
        return url


class ForcingMEPS(ForcingFetcher):
    def __init__(self, prefix='subset', stride=24, hours_per_file=24, last_file=None, lead_time=0):
        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.prefix = copy(prefix)
        self.last_file = copy(last_file)
        return

    def __call__(self, grid, start_time, end_time, expansion_factor):
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time

        #days = bnd.day_list(start_time = self.start_time, end_time = self.end_time)
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
            #print(time_stamp)
            #print(days[n].strftime('%Y-%m-%d'))
            url = self.get_url(file_times[n], self.prefix)

            msg.info(url)
            msg.plain(
                f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")

            nc_fimex = f'dnora_wnd_temp/wind_{n:04.0f}.nc'
            #nc_fimex = 'dnora_wnd_temp.nc'

            # Define area to search in
            expand_lon = (grid.lon()[-1] - grid.lon()
                          [0])*(expansion_factor-1)*0.5
            expand_lat = (grid.lat()[-1] - grid.lat()
                          [0])*(expansion_factor-1)*0.5

            lon_min = grid.lon()[0] - expand_lon
            lon_max = grid.lon()[-1] + expand_lon

            lat_min = grid.lat()[0] - expand_lat
            lat_max = grid.lat()[-1] + expand_lat

            # Temporary hack: set resolution to about 2.5 km
            dlat = 2.5/111
            dlon = dlat*2

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
                             #'--extract.reduceToBoundingBox.south=' + str(lat_min),
                             #'--extract.reduceToBoundingBox.north=' + str(lat_max),
                             #'--extract.reduceToBoundingBox.west=' + str(lon_min),
                             #'--extract.reduceToBoundingBox.east=' + str(lon_max),
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
        filename = 'meps_'+prefix+'_2_5km_'+time_stamp.strftime('%Y')+time_stamp.strftime(
            '%m')+time_stamp.strftime('%d')+'T'+time_stamp.strftime('%H')+'Z.nc'
        url = 'https://thredds.met.no/thredds/dodsC/meps25epsarchive/'+time_stamp.strftime(
            '%Y')+'/'+time_stamp.strftime('%m')+'/'+time_stamp.strftime('%d')+'/' + filename
        #url = 'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/'+days[i].strftime('%Y')+'/'+days[i].strftime('%m')+'/'+days[i].strftime('%Y%m%d')+'_MyWam3km_hindcast.nc'
        #url = 'https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/'+day.strftime('%d')+'/MyWave_wam4_SPC_'+day.strftime('%Y%m%d')+'T00Z.nc'
        return url




# =============================================================================
# OUTPUT MODEL CLASSES RESPONSIBLE FOR WRITING SPECTRA IN CORRECT FORMAT
# =============================================================================


class OutputModel(ABC):
    @abstractmethod
    def __call__(self, forcing_out):
        pass


class DumpToNc(OutputModel):
    def __init__(self):
        pass

    def __call__(self, forcing_out: Forcing):
        msg.header(f"Writing output with {type(self).__name__}")
        output_file = f"wind_{forcing_out.name}_{forcing_out.grid.name()}_{str(forcing_out.time()[0])[0:10]}_{str(forcing_out.time()[-1])[0:10]}.nc"
        msg.to_file(output_file)
        forcing_out.data.to_netcdf(output_file)

        return


class OutputSWANascii(OutputModel):
    def __init__(self):
        pass

    def __call__(self, forcing_out: Forcing):
        msg.header(
            f'{type(self).__name__}: writing wind forcing from {forcing_out.name}')

        days = forcing_out.days()
        output_file = f"{forcing_out.grid.name()}_wind{days[0].strftime('%Y%m%d')}_{days[-1].strftime('%Y%m%d')}.asc"
        msg.info(f'Writing wind forcing to: {output_file}')
        #print(output_file)
        with open(output_file, 'w') as file_out:
            for day in days:
                msg.plain(day.strftime('%Y-%m-%d'))
                times = forcing_out.times_in_day(day)
                for n in range(len(times)):
                    time_stamp = pd.to_datetime(
                        times[n]).strftime('%Y%m%d.%H%M%S')+'\n'
                    file_out.write(time_stamp)
                    np.savetxt(file_out, forcing_out.u()
                               [n, :, :]*1000, fmt='%i')
                    file_out.write(time_stamp)
                    np.savetxt(file_out, forcing_out.v()
                               [n, :, :]*1000, fmt='%i')


#                     file_out.write(time_stamp)
#                     file_out.write(str(data.time.time[time_step].values).split('-')[0]+str(data.time.time[time_step].values).split('-')[1]+\
#                                    str(data.time.time[time_step].values).split('-')[2][:2]+'.'+str(data.time.time[time_step].values).split('-')[2][3:5]+'0000\n')
#                     np.savetxt(file_out,u[time_step]*1000, fmt='%i') #
#                     file_out.write(str(data.time.time[time_step].values).split('-')[0]+str(data.time.time[time_step].values).split('-')[1]+\
#                                    str(data.time.time[time_step].values).split('-')[2][:2]+'.'+str(data.time.time[time_step].values).split('-')[2][3:5]+'0000\n')
#                     np.savetxt(file_out,v[time_step]*1000, fmt='%i') #
# =============================================================================

        return
