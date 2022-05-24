#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 12:33:03 2022

@author: emiliebyermoen
"""
import xarray as xr
from abc import ABC, abstractmethod


def is_in_area(lon, lat, lonq, latq):
    """ 
    Lon = (lon_min,lon:max), lat = (lat_min,lat_max), lonq=value, latq=value
    If input values lon and lat is Nona, the function returns True. 
    """
    if lon is not None: 
        if lonq<lon[0] or lon[1]<lonq:
            return False
    if lat is not None: 
        if latq<lat[0] or lat[1]<latq:
            return False
    return True


class TimeSeriesReader(ABC): 
    """ Reads and returns time series withing given time interval
    """
    @abstractmethod
    def is_available(self, time_interval, elements, lon: tuple=None, lat: tuple=None):
        """
        time_interval:  tuple=(start, end) on 'YYYY/MM/DDTHH:mm' or numpy.datetime64 format
        elements:       list of elements to look for in data (str)
        
        Returns dictionary available with 'time', 'wave_elements' and 'wind_elements'
        """
        return available 
    
    @abstractmethod
    def fetch_data(self, time_interval, wave_elements:list, wind_elements:list):
        """
        time_interval:      array of numpy.datetime64 or tuple=(start, end)
        wave_elements:      list/array (str)
        wind_elements:      list/array (str)
        
        Returns data on xarray.Dataset or xarray.DataArray format. 
        """
        return data


class E39BuoyFetcher(TimeSeriesReader):
    """
    A class that searches for and returns avaiable data 
    within given boundaries for the buoy data
    from the E39 project made availiable in Thredds. 
    
    Initiate with one of the following buoys 
    'G':'G_Halsafjorden'
    'G1':'G1_Halsafjorden'
    'G2':'G2_Halsafjorden'
    'F':'F_Vartdalsfjorden'
    'D':'D_Breisundet'
    'C':'C_Sulafjorden'
    'C1':'C1_Sulafjorden'
    'B':'B_Sulafjorden'
    'B1':'B1_Sulafjorden'
    'A':'A_Sulafjorden'
    """
    
    def __init__(self,location: str):
        self.location = location
        
    def get_url(self, year, month, data):
        """
        Returns url to file on Thredds for the given bouy, year and month.
        Data should be 'wave' or 'wind'. 
        """
        locations = {'G':'G_Halsafjorden',
                     'G1':'G1_Halsafjorden',
                     'G2':'G2_Halsafjorden',
                     'F':'F_Vartdalsfjorden',
                     'D':'D_Breisundet',
                     'C':'C_Sulafjorden',
                     'C1':'C1_Sulafjorden',
                     'B':'B_Sulafjorden',
                     'B1':'B1_Sulafjorden',
                     'A':'A_Sulafjorden'
                     }
        return f'https://thredds.met.no/thredds/dodsC/obs/buoy-svv-e39/{year}/{month}/{year}{month}_E39_{locations[self.location]}_{data}.nc'
    
    def get_aggregated_url(self):
        """
        Returns url to file with aggregated data information on Thredds.
        """
        locations = {'G':'G_Halsafjorden',
                     'G1':'G1_Halsafjorden',
                     'G2':'G2_Halsafjorden',
                     'F':'F_Vartdalsfjorden',
                     'D':'D_Breisundet',
                     'C':'C_Sulafjorden',
                     'C1':'C1_Sulafjorden',
                     'B':'B_Sulafjorden',
                     'B1':'B1_Sulafjorden',
                     'A':'A_Sulafjorden'
                     }
        return f'https://thredds.met.no/thredds/dodsC/obs/buoy-svv-e39/Aggregated_BUOY_observations/{locations[self.location]}_wave.ncml'
    
    def position(self):
        """
        Returns position of bouy given in the aggregated dataset on Thredds.
        
        Noe: I haven't used this, but I thought it might be helpful
        if we wanted to make the program choose which bouy to take. 
        """
        ds = xr.open_dataset(self.get_aggregated_url)
        lon = (ds.geospatial_lon_min+ds.geospatial_lon_max)/2 
        lat = (ds.geospatial_lat_min+ds.geospatial_lat_max)/2 
        return lon, lat
    
    def is_available(self, time_interval, elements, lon: tuple=None, lat: tuple=None):
        """ 
        Returns dictionary with (1) array of points in time within 
        the given time interval (for given bouy contained in Thredds, see catalog:
        https://thredds.met.no/thredds/catalog/obs/buoy-svv-e39/Aggregated_BUOY_observations/catalog.html)
        and (2) list of given elements that exist in the Thredds data for given bouy.
        
        Returns None if there is no data matching the reqirements.
        
        Time input should be given on the fomat:
            'YYYY/MM/DDTHH:mm'
        
        Elements in the Thredds E39 project bouy data: 
        'latitude'   Degree north (-90,90)
        'longitude'  Degree east (-180,180)
        
        Wave elements:
        'Hm0'        Significant wave height estimate from spectrum
        'hm0a'       Low frequency band Significant wave height estimated from spectrum
        'hm0b'       High frequency band Significant wave height estimated from spectrum
        'hmax'       Maximum individual wave height from zero crossing analysis
        'tm02'       Mean wave period estimated from 0th and 2nd moment of spectrum
        'tm02a'      Low frequency mean wave period estimated from 0th and 2nd moment of spectrum
        'tm02b'      High frequency mean wave period estimated from 0th and 2nd moment of spectrum
        'tp'         Spectral peak wave period
        'thmax'      Period of highest individual wave from zero crossing analysis
        'tm01'       Mean wave period estimated from 0th and 1st moment of spectrum
        'thhf'       Wave direction of wind sea (from)
        'mdir'       Mean wave direction (from) calculated from directional spectrum
        'mdira'      Low frequency mean wave direction (from) calculated from directional spectrum
        'mdirb'      High frequency mean wave direction (from) calculated from directional spectrum
        'thtp'       Mean wave direction at the peak of the heave variance spectrum
        'sprtp'      Wave spreading at the peak of the heave variance spectrum
        
        Wind elements: 
        'WindSpeed'      Wind speed (m/s)
        'WindGust'       3 second Wind Gust (m/s)
        'WindDirection'  Wind from direction [0,360]
        """

        ds = xr.open_dataset(self.get_aggregated_url())
        
        # Checks if the requested elements are in the Thredds data
        wave_elements = []
        wind_elements = []
        for element in elements:
                if element in ds.data_vars:
                    wave_elements.append(element) 
                elif element in ['WindSpeed','WindGust','WindDirection']: 
                    wind_elements.append(element)
            
        if len(wave_elements)>0 or len(wind_elements)>0:
            # Checks if lat, lon of the bouy is within the given lat, lon boundaries
            # The is_in_area function returns True if boundaries are given as None 
            if is_in_area(lon, lat, ds.geospatial_lon_min, ds.geospatial_lat_min):
            # The following alternative to that over would make the code slower and is not nessicary here I think:
            # if is_in_area(lon, lat, self.position) 
                # Collects array of all the points in time that are contained within the time boundary given
                available_time = ds.sel(time=slice(time_interval[0], time_interval[-1]))['time'].values
            else: 
                available_time = None

        return {'time': available_time, 'wave_elements': wave_elements, 
                'wind_elements': wind_elements}

    def fetch_data(self, time_interval, wave_elements: list=[], wind_elements: list=[]):
        """
        Returns xarray.Dataset for given time interval with given elemnts
        and given bouy (from the E39 bouy stations).
        
        Input date array time_interval should be on np.datetime64 format.
        """
        fetch_wave = len(wave_elements)>0
        fetch_wind = len(wind_elements)>0
        
        data = {}
        yearmonths = []
        for timestep in time_interval: 
            year, split, month = str(timestep.astype('datetime64[M]')).partition('-')
            if (year, month) not in yearmonths:
                if fetch_wave:
                    ds_wave = xr.open_dataset(self.get_url(year, month, 'wave'))
                    for element in wave_elements:
                        data[f'{element}_{year}/{month}'] = ds_wave[element].sel(time=slice(time_interval[0],time_interval[-1]))
                if fetch_wind: 
                    ds_wind = xr.open_dataset(self.get_url(year, month, 'wind'))
                    for element in wind_elements: 
                        data[f'{element}_{year}/{month}'] = ds_wind[element].sel(time=slice(time_interval[0],time_interval[-1]))
                yearmonths.append((year, month))
                
        combined = xr.merge(data[file] for file in data.keys())
        return combined
        
class KystverketBuoyFetcher(TimeSeriesReader):
    """
    Searches for and returns available data within given
    boundaries from Kystverket buoys available at Thredds.
        
    See catalog: https://thredds.met.no/thredds/catalog/obs/kystverketbuoy/catalog.html
        
    Initiate with one of the following bouys:
    'F':    Fauskane
    'V':    Vestfjorden
    """
    
    def __init__(self, location: str):
        self.location = location
    
    def get_url(self, year, month, data):
        if data == 'wave':
            sensor = 'AanderaaMotusSensor'
        elif data == 'wind':
            sensor = 'GillWindSensor'
    
        locations = {'F':'Fauskane',
            'V':'Vestfjorden'}
        
        return f'https://thredds.met.no/thredds/dodsC/obs/kystverketbuoy/{year}/{month}/{year}{month}_Kystverket-Smartbuoy-{locations[self.location]}_{sensor}.nc'

    def position(self):
        """
        Returns position of bouy given in the aggregated dataset on Thredds.
        
        Noe: I haven't used this, but I thought it might be helpful
        if we wanted to make the program choose which bouy to take.
        """
        ds = xr.open_dataset(self.get_url,'2020','01','wave')
        lon = (ds.geospatial_lon_min+ds.geospatial_lon_max)/2
        lat = (ds.geospatial_lat_min+ds.geospatial_lat_max)/2
        return lon, lat
        
        
    def is_available(self, time_interval, elements, lon: tuple=None, lat: tuple=None):
        """
            Returns dictionary with
            1) tuple of time (start,stop) with None if the start time or end time does not exist
            2) list of wave elements avaiable (out of those given)
            3) list of wind elements avaiable (out of those given)
            
            Time input should be given on the fomat:
            'YYYY/MM/DDTHH:mm'
            
            Takes lat=(min, max), lon=(min,max), date=(first, last), and list of elements from lists below:
            
            Elements:
            'longitude'                    longitude
            'latitude'                     latitude
            'netTime'                      Time of measurement
            'latitude_qc'                  Latitude Quality
            'longitude_qc'                 Longitude Quality
            'RecordNumber'                 Record Number
            'Input_Voltage_Logger'         Input voltage on datalogger
            'ActualInterval'               Actual interval on datalogger in second
            
            Wind elements:
            'Gust_Direction'               Wind direction gust (3 second)
            'Gust_Speed'                   Wind speed gust (3 second)
            'Average_Wind_Direction'       Wind direction (average 10 minute)
            'Average_Wind_Speed'           Wind speed (average 10 minute)
            
            Wave elements:
            'Input_Current'                Input Current
            'Input_Voltage'                Input Voltage
            'StDev_Roll'                   Standard deviation roll
            'StDev_Pitch'                  Standard deviation pitch
            'StDev_Heading'                Standard deviation heading
            'Roll'                         Roll
            'Pitch'                        Pitch
            'Heading'                      Heading
            'Long_Crestedness_Parameters'  Long Crestedness Parameters
            'First_Order_Spread'           First Order Spread
            'Mean_Spreading_Angle'         Mean Spreading Angle
            'Wave_Period_Tz'               Sea surface wave zero upcrossing period
            'Wave_Period_Tmax'             Wave Period Tmax (The corresponding period of the wave that is identified as wave height max)
            'Wave_Height_Trough'           Wave Height Trough
            'Wave_Height_Crest'            Wave Height Crest
            'Wave_Height_Hmax'             Maximum individual wave height from zero crossing analysis
            'Wave_Height_Wind_Hm0'         Wave Height Wind Hm0
            'Wave_Height_Swell_Hm0'        Wave Height Swell Hm0
            'Wave_Peak_Period_Wind'        Wave Peak Period Wind
            'Wave_Peak_Period_Swell'       Wave Peak Period Swell
            'Wave_Peak_Period'             Spectral peak wave period
            'Wave_Mean_Period_Tm02'        Mean wave period estimated from 0th and 2nd moment of spectrum
            'Wave_Peak_Direction_Wind'     Wave Peak Direction Wind
            'Wave_Mean_Direction'          Wave Mean Direction
            'Wave_Peak_Direction_Swell'    Wave Peak Direction Swell
            'Wave_Peak_Direction'          Wave Peak Direction
            'Significant_Wave_Height_Hm0'  Significant Wave Height Hm0 estimate from spectrum
            """
        # Initiates start_time and end_time = None, so if the functios does not find
        # those times in the datasets, it will return None
        start_time = None
        end_time = None
        start_year = time_interval[0][0:4]
        start_month = time_interval[0][5:7]
        if len(start_month)==0:
            start_month = '01'
        end_year = time_interval[-1][0:4]
        end_month = time_interval[-1][5:7]
        if len(end_month)==0:
            end_month = '12'

        # Checks if first point of time interval do exist as file in Thredds
        start_ds = xr.open_dataset(self.get_url(start_year, start_month, 'wave'))
        start_time = start_ds.sel(time=slice(time_interval[0], time_interval[-1]))['time'].values[0]
        
        # Checks if the requested elements are in the Thredds data
        wave_elements = []
        wind_elements = []
        for element in elements:
            if element in start_ds.data_vars:
                wave_elements.append(element)
            elif element in ['Gust_Direction','Gust_Speed','Average_Wind_Direction','Average_Wind_Speed']:
                wind_elements.append(element)
        
        if len(wave_elements)>0 or len(wind_elements)>0:
            # Checks if lat, lon of the bouy is within the given lat, lon boundaries
            if is_in_area(lon, lat, start_ds.geospatial_lon_min, start_ds.geospatial_lat_min):
                # The following alternative to that over would make the code slower and is not nessicary here I think:
                # if is_in_area(lon, lat, self.position)
                # checks if the last given timestep is included in Thredds
                if start_year == end_year and start_month == end_month:
                    end_ds = start_ds
                else:
                    end_ds = xr.open_dataset(self.get_url(end_year, end_month, 'wave'))
                end_time = end_ds.sel(time=slice(time_interval[0],time_interval[-1]))['time'].values[-1]

        return {'time': (start_time, end_time), 'wave_elements': wave_elements, 'wind_elements': wind_elements}
    
    def fetch_data(self, time_interval, wave_elements: list=[], wind_elements: list=[]):
        """
        Returns xarray.Dataset for given time interval with given elemnts
        and given bouy (from the E39 bouy stations).
        Input time_interval should be on format (start, end) should be on np.datetime64 format. (start,end)
        """
        start_year = str(time_interval[0])[0:4]
        start_month = str(time_interval[0])[5:7]
        end_year = str(time_interval[-1])[0:4]
        end_month = str(time_interval[-1])[5:7]
        
        yearmonths = []
        if start_year == end_year:
            if start_month == end_month:
                yearmonths.append((start_year,start_month))
            else:
                for month in range(int(start_month),int(end_month)+1):
                    if month<10:
                        yearmonths.append((start_year,f'0{month}')) #need a 0 before the number to get right syntax for creating url later
                    else:
                        yearmonths.append((start_year,f'{month}'))
        else:
            for month in range(int(start_month),13):
                if month<10:
                    yearmonths.append((start_year, f'0{month}'))
                else:
                    yearmonths.append((start_year, f'{month}'))
            if int(end_year)-int(start_year)>1:
                for year in range(int(start_year)+1,int(end_year)):
                    for month in range(1,10):
                        yearmonths.append((f'{year}',f'0{month}'))
                    for month in range(10,13):
                        yearmonths.append((f'{year}',f'{month}'))
            for month in range(1,int(end_month)+1):
                if month<10:
                    yearmonths.append((end_year,f'0{month}'))
                else:
                    yearmonths.append((end_year,f'{month}'))
        
        fetch_wave = len(wave_elements)>0
        fetch_wind = len(wind_elements)>0
        data = {}
        for (year, month) in yearmonths:
            if fetch_wave:
                ds_wave = xr.open_dataset(self.get_url(year, month, 'wave'))
                for element in wave_elements:
                    data[f'{element}_{year}/{month}'] = ds_wave[element].sel(time=slice(time_interval[0],time_interval[-1]))
            if fetch_wind:
                ds_wind = xr.open_dataset(self.get_url(year, month, 'wind'))
                for element in wind_elements:
                    data[f'{element}_{year}/{month}'] = ds_wind[element].sel(time=slice(time_interval[0],time_interval[-1]))

        combined = xr.merge(data[file] for file in data.keys())
        
        return combined

    
    
    
    
