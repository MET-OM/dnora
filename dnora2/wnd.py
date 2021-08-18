from abc import ABC,  abstractmethod
from dnora2.bnd import day_list
from dnora2 import msg
from copy import copy
import numpy as np
import pandas as pd
from subprocess import call
import xarray as xr
import os

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
    def __init__(self, prefix = 'subset', hours_per_file = 24, lead_time = 0):
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.prefix = copy(prefix)
        return
    
    def create_time_stamps(self, start_time, end_time, hours_per_file):
        
        h0 = int(start_time.split('T')[1][0:2])
        h0 = int(np.ceil(h0/hours_per_file)*hours_per_file)
        
        h1 = int(end_time.split('T')[1][0:2])
        h1 = int(np.ceil(h1/hours_per_file)*hours_per_file)
        
        t0 = np.datetime64(start_time.split('T')[0] + 'T00:00:00')+np.timedelta64(h0,'h')
        t1 = np.datetime64(end_time.split('T')[0] + 'T00:00:00')+np.timedelta64(h1,'h')
        
        time_stamps = pd.date_range(start = t0, end = t1, freq=f'{hours_per_file}H')
        return time_stamps
    
    def __call__(self, grid, start_time, end_time, expansion_factor):
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        self.start_time = start_time
        self.end_time = end_time
              
        
        #days = day_list(start_time = self.start_time, end_time = self.end_time)
        time_stamps = self.create_time_stamps(start_time, end_time, self.hours_per_file)
        
        file_stamps = time_stamps - np.timedelta64(self.lead_time,'h')
        wnd_list =[]
        
        temp_folder = 'dnora_wnd_temp'
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print ("Creating folder %s..." % temp_folder)
        
        
        for n in range(len(file_stamps)):
            #print(time_stamp)
            #print(days[n].strftime('%Y-%m-%d')) 
            url = self.get_url(file_stamps[n], self.prefix)
            
            start_date_fimex = time_stamps[n].strftime('%Y-%m-%dT%H:%M:%S')
            if n == (len(time_stamps)-1):
                end_date_fimex = end_time
                #end_date_fimex = (time_stamps[n] + np.timedelta64(0, 'h')).strftime('%Y-%m-%dT%H:%M:%S')
            else:
                end_date_fimex = (time_stamps[n] + np.timedelta64(self.hours_per_file-1, 'h')).strftime('%Y-%m-%dT%H:%M:%S')
                
            msg.info(url)
            msg.plain(f"Reading wind forcing data: {start_date_fimex}-{end_date_fimex}")
            
            
            
            
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
                  '--extract.reduceTime.start='+start_date_fimex,'--extract.reduceTime.end='+end_date_fimex,
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
        self.data = forcing_fetcher(self.grid, start_time, end_time, expansion_factor)
        return
    
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