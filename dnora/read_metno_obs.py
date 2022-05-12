#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 12:33:03 2022

@author: emiliebyermoen
"""

import xarray as xr


def is_in_area(lon, lat, lonq, latq):
    """ Lon = (lon_min,lon:max), lat = (lat_min,lat_max), lonq=value, latq=value
    """
    if lon is not None: 
        if lonq<lon[0] or lon[1]<lonq:
            return False
    if lat is not None: 
        if latq<lat[0] or lat[1]<latq:
            return False
    return True


class ThreddsFetcher:
    """
    A class that searches for and returns avaiable data 
    within given boundaries for the bouy data 
    from the E39 project made availiable in Thredds. 
    
    Initiate with one of the following bouys 
    'G':'G_Halsafjorden'
    'F':'F_Vartdalsfjorden'
    'D':'D_Breisundet'
    'C':'C_Sulafjorden'
    'B':'B_Sulafjorden'
    """
    
    def __init__(self,location: str):
        self.location = location
        
    def get_url(self, year=None, month=None):
        """
        Returns url to file on Thredds for the given bouy, year and month.
        """
        locations = {'G':'G_Halsafjorden',
                     'F':'F_Vartdalsfjorden',
                     'D':'D_Breisundet',
                     'C':'C_Sulafjorden',
                     'B':'B_Sulafjorden'
                     }
        return f'https://thredds.met.no/thredds/dodsC/obs/buoy-svv-e39/{year}/{month}/{year}{month}_E39_{locations[self.location]}_wave.nc'
    
    def get_aggregated_url(self):
        """
        Returns url to file with aggregated data information on Thredds.
        """
        locations = {'G':'G_Halsafjorden',
                     'F':'F_Vartdalsfjorden',
                     'D':'D_Breisundet',
                     'C':'C_Sulafjorden',
                     'B':'B_Sulafjorden'
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
    
    def is_available(self,lon: tuple=None, lat: tuple=None, date: tuple=None,
                     elements: list=None):
        """ 
        Returns dictionary with (1) array of points in time within 
        the given time interval (for given bouy contained in Thredds, see catalog:
        https://thredds.met.no/thredds/catalog/obs/buoy-svv-e39/Aggregated_BUOY_observations/catalog.html)
        and (2) list of given elements that exist in the Thredds data for given bouy.
        
        Returns None if there is no data matching the reqirements.
        
        Time input should be given on the fomat:
            'YYYY/MM/DDTHH:mm'
        
        Wave elements in the Thredds E39 project bouy data: 
        'latitude'   Degree north (-90,90)
        'longitude'  Degree east (-180,180)
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
        """

        ds = xr.open_dataset(self.get_aggregated_url())
        
        # Checks if the requested elements are in the Thredds data
        available_elements = []
        for element in elements:
                if element in ds.data_vars:
                    available_elements.append(element) 
    
        if len(available_elements) == 0: 
            available_elements = None
        else: 
            # Checks if lat, lon of the bouy is within the given lat, lon boundaries
            if is_in_area(lon, lat, ds.geospatial_lon_min, ds.geospatial_lat_min):
            # The following alternative to that over would make the code slower and is not nessicary here I think:
            # if is_in_area(lon, lat, self.position) 
                # Collects array of all the points in time that are contained within the time boundary given
                available_time = ds.sel(time=slice(date[0], date[1]))['time'].values
            else: 
                available_time = None

        return {'time': available_time, 'elements': available_elements}

    def fetch_data(self,lon: tuple=None, lat: tuple=None, date=None,
                     elements: list=None):
        """
        Returns xarray.Dataset for given time interval with given elemnts
        and given bouy (from the E39 bouy stations).
        
        Input date array should be on np.datetime64 format.
        """
        
        data = {}
        yearmonths = []
        for timestep in date: 
            year, split, month = str(timestep.astype('datetime64[M]')).partition('-')
            if (year, month) not in yearmonths:
                ds = xr.open_dataset(self.get_url(year=year, month=month))
                #ds = xr.open_dataset(f'https://thredds.met.no/thredds/dodsC/obs/buoy-svv-e39/{year}/{month}/{year}{month}_E39_D_Breisundet_wave.nc')
                yearmonths.append((year, month))
                for element in elements:
                    data[f'{element}_{year}/{month}'] = ds[element].sel(time=slice(date[0],date[-1]))
                    #data[f'{locations[self.location]}_{element}_{year}/{month}'] = ds[element].sel(time=slice(date[0],date[-1]))
        combined = xr.merge(data[file] for file in data.keys())
        return combined