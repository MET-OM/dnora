import numpy as np
import pandas as pd
import xarray as xr
from abc import ABC, abstractmethod
import netCDF4
from dnora2 import msg
from dnora2 import grd
from dnora2 import bnd
from dnora2 import spec
from copy import copy
from dnora2.bnd import day_list, month_list
import re
import matplotlib.pyplot as plt

# =============================================================================
# PARAMETER FETCHER CLASSES RESPONSIBLE FOR ACTUALLY READING THE SPECTRA
# =============================================================================
class ParameterFetcher(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def get_coordinates(self, start_time):
        pass

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

    def get_time_limits_month(self, ind):
        """Determines star and end time for the day. First and last day doesn't start at 00:00 or end at 23:59"""
        
        months = month_list(start_time = self.start_time, end_time = self.end_time)
        
        if ind == 0:
            t0 = self.start_time
            t1 = months[0].strftime('%Y-%m-%d') + "T23:59:59"
        elif ind == (len(months)-1):
            t0 = months[-1].strftime('%Y-%m') + "-01T00:00:00"				
            t1 = self.end_time
        else:
            t0 = months[ind].strftime('%Y-%m') + "-01T00:00:00"	
            t1 = months[ind].strftime('%Y-%m-%d') + "T23:59:59"
        return t0, t1


    def __str__(self):
        return (f"{self.start_time} - {self.end_time}")


class ParameterFromThredds(ParameterFetcher):
    def __init__(self, location = '', param_list = {}, dump_rest = True):
        self.location = copy(location)
        self.param_list = copy(param_list)
        self.dump_rest = copy(dump_rest)
        return
      
    def get_coordinates(self, start_time):
        day = day_list(start_time, start_time)
        url = self.get_url(self.location, day[0])
        
        data = xr.open_dataset(url).sel(time = slice(start_time, start_time))
        lon= data.longitude.values
        lat = data.latitude.values
        return lon, lat
    
    
    def __call__(self, start_time, end_time, inds):
        self.start_time = start_time
        self.end_time = end_time
        months = month_list(start_time, end_time)
        
        msg.info(f"Getting wave buoy data from thredds from {start_time} to {end_time}")
        bnd_list = []    
        for n in range(len(months)):
            url = self.get_url(self.location, months[n])
            msg.plain(url)
            t0, t1 = self.get_time_limits_month(n)
            bnd_list.append(xr.open_dataset(url).sel(time = slice(t0, t1)))
            
        data=xr.concat(bnd_list, dim="time")#.squeeze('y')
        N = len(data.time.values)
    
        time = data.time.values
        lon= np.reshape(data.longitude.values, (N,1))
        lat = np.reshape(data.latitude.values, (N,1))
        source = re.sub(r' ', '_', data.station_name)

        wavedata = {}
        list_of_keys = list(data.variables.keys())
        wanted_old = list(self.param_list.keys())
        for key in list_of_keys:
            if key not in ['time', 'longitude', 'latitude']:
                if key in wanted_old:
                    wavedata[self.param_list[key]] = np.reshape(data[key].values, (N,1))
                elif self.dump_rest:
                    wavedata[key] = np.reshape(data[key].values, (N,1))
    
        return time, wavedata, lon, lat, source
    
    def get_url(self, location, month):
        url = 'https://thredds.met.no/thredds/dodsC/obs/buoy-svv-e39/'+month.strftime('%Y')+'/'+month.strftime('%m')+'/'+month.strftime('%Y')+month.strftime('%m')+'_E39_'+location+'_wave.nc'
        #url = 'https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/'+day.strftime('%d')+'/MyWave_wam4_SPC_'+day.strftime('%Y%m%d')+'T00Z.nc'
        return url


class ParameterFromNc(ParameterFetcher):
    def __init__(self, filename = '', param_list = {}, dump_rest = True):
        self.filename = copy(filename)
        self.param_list = copy(param_list)
        self.dump_rest = copy(dump_rest)
        return
    
    def get_coordinates(self, start_time):
        data = xr.open_dataset(self.filename).sel(time = slice(start_time, start_time))
        lon= data.longitude.values
        lat = data.latitude.values
        return lon, lat

    def __call__(self, start_time, end_time, inds):
        data = xr.open_dataset(self.filename).sel(time = slice(start_time, end_time))
        N = len(data.time.values)
    
        time = data.time.values
        lon= np.reshape(data.longitude.values, (N,1))
        lat = np.reshape(data.latitude.values, (N,1))
        source = re.sub(r' ', '_', data.station_name)

        wavedata = {}
        list_of_keys = list(data.variables.keys())
        wanted_old = list(self.param_list.keys())
        for key in list_of_keys:
            if key not in ['time', 'longitude', 'latitude']:
                if key in wanted_old:
                    wavedata[self.param_list[key]] = np.reshape(data[key].values, (N,1))
                elif self.dump_rest:
                    wavedata[key] = np.reshape(data[key].values, (N,1))
    
        return time, wavedata, lon, lat, source

class ParameterForceFeed(ParameterFetcher):
    def __init__(self, time, wavedata, lon, lat):
        self.time = copy(time)
        self.wavedata = copy(wavedata)
        self.lon = copy(lon)
        self.lat = copy(lat)
        return 
    
    def get_coordinates(self, start_time):
        return copy(self.lon), copy(self.lat)
    
    def __call__(self, start_time, end_time, inds):
        return  copy(self.time), copy(self.wavedata), copy(self.lon), copy(self.lat), ''



# =============================================================================
#  PARAMETER OBJECT CONTAINING THE ACTUAL DATA
# =============================================================================

class Parameter:
    def __init__(self, grid = None, name = "AnonymousBoundary"):
        if grid is None:
            self.grid=grd.Grid()
        else:
            self.grid = copy(grid)
        self.name = name
        return


    def import_parameter(self, start_time: str, end_time: str, parameter_fetcher: ParameterFetcher,  point_picker = bnd.TrivialPicker()):
        self.start_time = copy(start_time)
        self.end_time = copy(end_time)
        
        msg.header(f"{type(parameter_fetcher).__name__}: Reading coordinats of spectra...")
        lon_all, lat_all = parameter_fetcher.get_coordinates(self.start_time)
        
        
        msg.header(f"Choosing spectra with {type(point_picker).__name__}")
        inds = point_picker(self.grid, lon_all, lat_all)
        
        msg.header(f"{type(parameter_fetcher).__name__}: Loading boundary spectra...")
        time, wavedata, lon, lat, source = parameter_fetcher(self.start_time, end_time, inds)

        self.data = self.compile_to_xr(time, wavedata, lon, lat, source)
        #self.mask = [True]*len(self.x())
        
        return
    
    
    def compile_to_xr(self, time, wavedata, lon, lat, source):
        x = np.array(range(lon.shape[1]))
        
        # Create data_vars dict
        data_vars={}
        keys=list(wavedata.keys())
        vals=list(wavedata.values())
        
        for n in range(len(wavedata)):
            data_vars[keys[n]] = {"dims": ["time", "x"], "data": vals[n]}
        
        d = {"coords": {"time": {"dims": "time", "data": time}, "x": {"dims": "x", "data": x}, "lon": {"dims": ["time", "x"], "data": lon}, "lat": {"dims": ["time", "x"], "data": lat}},
        "dims": ["time", "x"],
        "data_vars": data_vars,
        }
        data = xr.Dataset.from_dict(d)
        return data  
    
    def plot(self, param_list):
        for p in param_list:
            plt.figure()
            plt.plot(self.data.time.values, self.data[p].values)
        
        pass
    