#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 10:14:57 2022

@author: emiliebyermoen
"""
import utm 
import pandas as pd

class KartverketNo50m(TopoReader):
    """Reads data from Kartverket bathymetry.

    High resolution bathymetry dataset for the whole Norwegian Coast.
    Can be found at:  
    https://kartkatalog.geonorge.no/metadata/dybdedata-terrengmodeller-50-meters-grid-landsdekkende/bbd687d0-d34f-4d95-9e60-27e330e0f76e

    """

    def __init__(self, expansion_factor: float=1.2, utmzone: int=33,
                 tile: str='B1008', folder: str='/lustre/storeB/project/fou/om/WW3/bathy/kartverket_50m_x_50m') -> Tuple:
        self.source=f'{folder}/{tile}_grid50_utm33.xyz'
        self.expansion_factor = expansion_factor
        return

    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float):
        # Area is expanded a bit to not get in trouble in the meshing stage
        # when we interpoolate or filter
        lon0, lon1, lat0, lat1 = expand_area(lon_min, lon_max, lat_min, lat_max, self.expansion_factor)
        
        df = pd.read_csv(self.source, sep= ' ', header=None)
        df.columns = ['x','y','z']
        x = np.array(df['x'].astype(float))
        y = np.array(df['y'].astype(float))
        z = np.array(df['z'].astype(float))
        
        # Converting from utm to latitude and longitude 
        lat, lon = utm.to_latlon(df['x'], df['y'], utmzone, northern=True, strict=False)
        
        mask_lat = np.logical_and(lat0 < lat, lat < lat1)
        mask_lon = np.logical_and(lon0 < lon, lon < lon1)
        mask = np.logical_and(mask_lat, mask_lon)
    
        topo_lat = lat[mask]
        topo_lon = lon[mask]
        topo = z[mask]

        return topo, topo_lon, topo_lat

    def __str__(self):
        return(f"Reading Kartverket topography from {self.source}.")
