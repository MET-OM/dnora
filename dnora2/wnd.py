from abc import ABC,  abstractmethod
from dnora2.bnd import day_list
from dnora2 import msg
from copy import copy
import numpy as np
import pandas as pd
from subprocess import call
import xarray as xr
import os


def create_time_stamps(start_time, end_time, stride, hours_per_file = None, last_file = None, lead_time = 0):
    """Create time stamps to read in blocks of wind forcing from files"""
    if hours_per_file is None:
        hours_per_file = stride

    t0 = np.datetime64(start_time)
    if last_file is not None:
        t1 = np.datetime64(last_file)    
    else:
        t1 = np.datetime64(end_time)    
    
    file_times = pd.date_range(start = t0, end = t1, freq=f'{stride}H') - np.timedelta64(lead_time,'h')
    start_times = file_times
    end_times = start_times + np.timedelta64(stride-1,'h')
    
    # Make sure we don't exceed the user requested end time
    end_times.values[-1] = min([end_times.values[-1], np.datetime64(end_time)])
    
    # This is the hard upper limit for how much data we can get considering the existing last file and hours in one file
    hard_end_time = min([np.datetime64(end_time), np.datetime64(last_file) + np.timedelta64(hours_per_file,'h')])
    
    end_times.values[-1] = max([end_times.values[-1], hard_end_time])

    return start_times, end_times, file_times

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
        
        days = day_list(start_time = self.start_time, end_time = self.end_time)
        
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







class ForcingMEPS(ForcingFetcher):
    def __init__(self, prefix = 'subset', stride = 24, hours_per_file = 24, last_file = None, lead_time = 0):
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
              
        
        #days = day_list(start_time = self.start_time, end_time = self.end_time)
        start_times, end_times, file_times = create_time_stamps(start_time, end_time, self.stride, self.hours_per_file, self.last_file, self.lead_time)

        wnd_list =[]
        
        temp_folder = 'dnora_wnd_temp'
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print ("Creating folder %s..." % temp_folder)
        
        
        for n in range(len(file_times)):
            #print(time_stamp)
            #print(days[n].strftime('%Y-%m-%d')) 
            url = self.get_url(file_times[n], self.prefix)
            
            msg.info(url)
            msg.plain(f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")
            
            
            
            
            nc_fimex = f'dnora_wnd_temp/wind_{n:04.0f}.nc'
            #nc_fimex = 'dnora_wnd_temp.nc'
        
            
            # Define area to search in
            expand_lon = (grid.lon()[-1] - grid.lon()[0])*(expansion_factor-1)*0.5
            expand_lat = (grid.lat()[-1] - grid.lat()[0])*(expansion_factor-1)*0.5
        
            lon_min = grid.lon()[0] - expand_lon
            lon_max = grid.lon()[-1] + expand_lon
        
            lat_min = grid.lat()[0] - expand_lat
            lat_max = grid.lat()[-1] + expand_lat
            
            dlon = grid.data.dlon*5
            dlat = grid.data.dlat*5
            
            call(['fimex-1.6', '--input.file='+url,
                  '--interpolate.method=bilinear',
                  '--interpolate.projString=+proj=latlong +ellps=sphere +a=6371000 +e=0',
                  '--interpolate.xAxisValues='+str(lon_min)+','+str(lon_min+dlon)+',...,'+str(lon_max)+'',
                  '--interpolate.yAxisValues='+str(lat_min)+','+str(lat_min+dlat)+',...,'+str(lat_max)+'',
                  '--interpolate.xAxisUnit=degree', '--interpolate.yAxisUnit=degree',
                  #'--extract.reduceToBoundingBox.south=' + str(lat_min),
                  #'--extract.reduceToBoundingBox.north=' + str(lat_max),
                  #'--extract.reduceToBoundingBox.west=' + str(lon_min),
                  #'--extract.reduceToBoundingBox.east=' + str(lon_max),
                  '--process.rotateVector.all',
                  '--extract.selectVariables=x_wind_10m','--extract.selectVariables=y_wind_10m',
                  '--extract.selectVariables=latitude','--extract.selectVariables=longitude',
                  '--extract.reduceTime.start='+start_times[n].strftime('%Y-%m-%dT%H:%M:%S'),'--extract.reduceTime.end='+end_times[n].strftime('%Y-%m-%dT%H:%M:%S'),
                  #'--extract.reduceDimension.name=ensemble_member',
                  #'--extract.reduceDimension.start=1',
                  #'--extract.reduceDimension.end=1',
                  '--process.rotateVector.direction=latlon',
                  '--output.file='+nc_fimex])
            
            wnd_list.append(xr.open_dataset(nc_fimex).squeeze())
        
        wind_forcing = xr.concat(wnd_list, dim="time")
        
        #wind_forcing.to_netcdf('test.nc')
        return wind_forcing
       
        
       
    def get_url(self, time_stamp, prefix):
        filename = 'meps_'+prefix+'_2_5km_'+time_stamp.strftime('%Y')+time_stamp.strftime('%m')+time_stamp.strftime('%d')+'T'+time_stamp.strftime('%H')+'Z.nc'
        url = 'https://thredds.met.no/thredds/dodsC/meps25epsarchive/'+time_stamp.strftime('%Y')+'/'+time_stamp.strftime('%m')+'/'+time_stamp.strftime('%d')+'/' + filename
        #url = 'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/'+days[i].strftime('%Y')+'/'+days[i].strftime('%m')+'/'+days[i].strftime('%Y%m%d')+'_MyWam3km_hindcast.nc' 
        #url = 'https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/'+day.strftime('%d')+'/MyWave_wam4_SPC_'+day.strftime('%Y%m%d')+'T00Z.nc'
        return url

       
class Forcing:
    def __init__(self, grid, name = 'AnonymousForcing'):
        self.grid = copy(grid)
        self.name = name
    def import_forcing(self, start_time: str, end_time: str, forcing_fetcher: ForcingFetcher, expansion_factor = 1.2):
        self.start_time = copy(start_time)
        self.end_time = copy(end_time)
        
        msg.header(f"{type(forcing_fetcher).__name__}: Loading wind forcing...")
        self.data = forcing_fetcher(self.grid, start_time, end_time, expansion_factor)
        
        return

    def days(self):
        """Determins a Pandas data range of all the days in the time span."""
        days = day_list(start_time = self.start_time, end_time = self.end_time)
        return days
    
    def time(self):
        return copy(self.data.time.values)
    
    def u(self):
        return copy(self.data.x_wind_10m.values)
    def v(self):
        return copy(self.data.y_wind_10m.values)
    
    
    def slice_data(self, start_time: str = '', end_time: str = ''):
        if not start_time:
            # This is not a string, but slicing works also with this input
            start_time = self.time()[0] 
        
        if not end_time:
            # This is not a string, but slicing works also with this input
            end_time = self.time()[-1]
            
        sliced_data = self.data.sel(time=slice(start_time, end_time))
            
        return sliced_data
    
    
    def times_in_day(self, day):
        """Determines time stamps of one given day."""
        t0 = day.strftime('%Y-%m-%d') + "T00:00:00"	
        t1 = day.strftime('%Y-%m-%d') + "T23:59:59"
        
        times = self.slice_data(start_time = t0, end_time = t1).time.values
        return times  
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
        output_file = f"wind_{forcing_out.name}_{forcing_out.grid.name()}.nc"
        msg.to_file(output_file)
        forcing_out.data.to_netcdf(output_file)
            
        return
    
class OutputSWANascii(OutputModel):
    def __init__(self):
        pass
    
    def __call__(self, forcing_out: Forcing):
        msg.header(f'{type(self).__name__}: writing wind forcing from {forcing_out.name}')
        
        days = forcing_out.days()
        output_file = f"{forcing_out.grid.name()}_wind{days[0].strftime('%Y%m%d')}_{days[-1].strftime('%Y%m%d')}.asc"
        msg.info(f'Writing wind forcing to: {output_file}')
        #print(output_file)
        with open(output_file, 'w') as file_out:
            for day in days:
                msg.plain(day.strftime('%Y-%m-%d'))
                times = forcing_out.times_in_day(day)
                for n in range(len(times)):
                    time_stamp = pd.to_datetime(times[n]).strftime('%Y%m%d.%H%M%S')+'\n'
                    file_out.write(time_stamp)
                    np.savetxt(file_out,forcing_out.u()[n,:,:]*1000, fmt='%i') #
                    file_out.write(time_stamp)
                    np.savetxt(file_out,forcing_out.v()[n,:,:]*1000, fmt='%i') #
                    


#                     file_out.write(time_stamp)
#                     file_out.write(str(data.time.time[time_step].values).split('-')[0]+str(data.time.time[time_step].values).split('-')[1]+\
#                                    str(data.time.time[time_step].values).split('-')[2][:2]+'.'+str(data.time.time[time_step].values).split('-')[2][3:5]+'0000\n')
#                     np.savetxt(file_out,u[time_step]*1000, fmt='%i') #
#                     file_out.write(str(data.time.time[time_step].values).split('-')[0]+str(data.time.time[time_step].values).split('-')[1]+\
#                                    str(data.time.time[time_step].values).split('-')[2][:2]+'.'+str(data.time.time[time_step].values).split('-')[2][3:5]+'0000\n')
#                     np.savetxt(file_out,v[time_step]*1000, fmt='%i') #
# =============================================================================
            
            
            
            

                    
            
        return
