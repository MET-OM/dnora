import os
import xarray as xr
import glob
import numpy as np
from copy import copy
from abc import ABC, abstractmethod
from typing import Tuple
import pandas as pd
# Import abstract classes and needed instances of them
from .read import BoundaryReader
import cdsapi
# Import aux_funcsiliry functions
from .. import msg
from ..aux_funcs import create_time_stamps, expand_area, day_list


def renormalize_era5_spec(bnd_spec):
    bnd_spec = bnd_spec.assign_coords(direction=np.arange(7.5, 352.5 + 15, 15))
    bnd_spec = bnd_spec.assign_coords(frequency=np.full(30, 0.03453) * (1.1 ** np.arange(0, 30)))
    bnd_spec = 10 ** bnd_spec
    bnd_spec = bnd_spec.fillna(0)
    return bnd_spec

def reshape_bnd_spec(bnd_spec):
    pass


def download_era5_from_cds(start_time, end_time, lon, lat, dlon, dlat, folder='dnora_bnd_temp') -> str:
    """Downloads ERA5 spectral data from the Copernicus Climate Data Store for a
    given area and time period"""
    start_time = pd.Timestamp(start_time)
    end_time = pd.Timestamp(end_time)
    c = cdsapi.Client()

    days = day_list(start_time, end_time)
    #tt = pd.date_range(start=start_time,end=end_time)

    filename = f'{folder}/EC_ERA5.nc'

    # Create string for dates
    #dates = [days[0].strftime('%Y-%m-%d'), days[-1].strftime('%Y-%m-%d')]
    #dates = '/'.join(dates)
    dates = f'{str(start_time)[0:10]}/to/{str(end_time)[0:10]}'


    cds_command ={
        'class': 'ea',
        'date': dates,
        'direction': '1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24',
        'domain': 'g',
        'area': f'{lat[1]}/{lon[0]}/{lat[0]}/{lon[1]}', # north, west, south, east
        'grid': f'{dlat}/{dlon}',
        'expver': '1',
        'frequency': '1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30',
        'param': '251.140',
        'stream': 'wave',
        'time': '00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00',
        'type': 'an',
        'format': 'netcdf',
    }
    #print(cds_command)

    # cds_command ={
    #     'class': 'ea',
    #     'date': dates,
    #     'direction': '/'.join([f'{n+1:01.0f}' for n in range(24)]),
    #     'domain': 'g',
    #     'expver': '1',
    #     'frequency': '/'.join([f'{n+1:01.0f}' for n in range(30)]),
    #     'param': '251.140',
    #     'stream': 'wave',
    #     'time': '00:00:00/03:00:00/06:00:00/09:00:00/12:00:00/15:00:00/18:00:00/21:00:00',
    #     'area': f'{lat[1]}/{lon[0]}/{lat[0]}/{lon[1]}', # north, west, south, east
    #     'grid': f'{dlon}/{dlat}',
    #     'type': 'an',
    #     'format': 'netcdf',
    #     }

    c.retrieve('reanalysis-era5-complete', cds_command, filename)
    return filename
class ERA5(BoundaryReader):

    def __init__(self):
        self.dlon = 0.5
        self.dlat = 0.5

    def convention(self) -> str:
        return 'Ocean'

    def get_coordinates(self, start_time) -> Tuple:
        """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""
        restricted_area = self.get_restricted_area()
        restricted_area.set_spacing(dlon=self.dlon, dlat=self.dlat)
        point_list = restricted_area._point_list()
        lon_all = point_list[:,0]
        lat_all = point_list[:,1]

        return lon_all, lat_all




    def __call__(self, start_time, end_time, inds) -> Tuple:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        msg.info(
            f"Getting ERA5 boundary spectra from {start_time} to {end_time}")


        temp_folder = 'dnora_bnd_temp'
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        local_read = False

        if not local_read:
            msg.plain("Removing old files from temporary folder...")
            for f in glob.glob(f"{temp_folder}/EC_ERA5.nc"):
                os.remove(f)

        restricted_area = self.get_restricted_area()
        lon = np.floor(np.array(restricted_area.lon_edges())/self.dlon)*self.dlon
        lat = np.floor(np.array(restricted_area.lat_edges())/self.dlat)*self.dlat



        if local_read:
            nc_file = f'{temp_folder}/EC_ERA5.nc'
        else:
            nc_file = download_era5_from_cds(start_time, end_time,
                                            lon=(lon[0], lon[1]),
                                            lat=(lat[0], lat[1]),
                                            dlon=self.dlon,
                                            dlat=self.dlat,
                                            folder=temp_folder)

        bnd_spec = xr.open_dataset(nc_file)
        bnd_spec = bnd_spec.sortby("time")
        bnd_spec = renormalize_era5_spec(bnd_spec)

        lon, lat = np.meshgrid(bnd_spec.longitude.values, bnd_spec.latitude.values[::-1])
        lon = lon.ravel()
        lat = lat.ravel()

        # This spec is time, freq, dir, lat, lon
        spec = bnd_spec.d2fd.values
        # Latitude was flipped to be ascending, so flip that dimension
        spec = np.flip(spec, 3)

        # This is time, freq, dir, station
        spec = np.reshape(spec, (len(bnd_spec.time),len(bnd_spec.frequency),len(bnd_spec.direction),len(lon)))
        # This is time, station, freq, dir (as we want it)
        spec = np.moveaxis(spec,3,1)

        freq = bnd_spec.frequency.values
        dirs = bnd_spec.direction.values
        time = bnd_spec.time.values

        source = 'ECMWF-ERA5 from Copernicus Climate Data Store'

        # Inds given by point picker
        lon = lon[inds]
        lat = lat[inds]
        spec = spec[:,inds,:,:]

        return  time, freq, dirs, spec, lon, lat, source
