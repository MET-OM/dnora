from abc import ABC,  abstractmethod
from copy import copy
import numpy as np
import xarray as xr
from subprocess import call
import os, glob
import time
import cdsapi
# Import objects
from ..grd.grd_mod import Grid

# Import abstract classes
from .read import ForcingReader

# Import aux_funcsiliry functions
from .. import msg
from ..aux_funcs import create_time_stamps, u_v_from_dir, expand_area, lon_in_km, create_monthly_time_stamps, day_list
import pandas as pd



def download_era5_from_cds(start_time, end_time, lon, lat, folder='dnora_wnd_temp') -> str:
    """Downloads ERA5 10 m wind data from the Copernicus Climate Data Store for a
    given area and time period"""
    start_time = pd.Timestamp(start_time)
    end_time = pd.Timestamp(end_time)
    c = cdsapi.Client()

    filename = f'{folder}/EC_ERA5.nc'
    # cds_command_test = {
    #     'product_type': 'reanalysis',
    #     'format': 'netcdf',
    #     'variable': [
    #         '10m_u_component_of_wind', '10m_v_component_of_wind',
    #     ],
    #     'year': '2008',
    #     'month': '01',
    #     'day': [
    #         '01', '02',
    #     ],
    #     'time': [
    #         '00:00', '01:00', '02:00',
    #         '03:00', '04:00', '05:00',
    #         '06:00', '07:00', '08:00',
    #         '09:00', '10:00', '11:00',
    #         '12:00', '13:00', '14:00',
    #         '15:00', '16:00', '17:00',
    #         '18:00', '19:00', '20:00',
    #         '21:00', '22:00', '23:00',
    #     ],
    #     'area': [
    #         61.25, 4, 60.53,
    #         5.73,
    #     ],
    # }

    # years = [f'{y:4.0f}' for y in int_list_of_years(start_time, end_time)]
    # if len(years) == 1:
    #     years = years[0]
    # months = [f'{m:02.0f}' for m in int_list_of_months(start_time, end_time)]
    # if len(months) == 1:
    #     months = months[0]
    # days = [f'{d:02.0f}' for d in int_list_of_days(start_time, end_time)]
    # if len(days) == 1:
    #     days = days[0]

    days = day_list(start_time, end_time)
    # Create string for dates
    dates = [days[0].strftime('%Y-%m-%d'), days[-1].strftime('%Y-%m-%d')]
    dates = '/'.join(dates)

    cds_command = {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': [
            '10m_u_component_of_wind', '10m_v_component_of_wind',
        ],
        'date': dates,
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
        'area': [
            lat[1], lon[0], lat[0],
            #lat[1], 4, lat[0],
            lon[1],
        ],
    }

    c.retrieve('reanalysis-era5-single-levels', cds_command, filename)
    return filename


class ERA5(ForcingReader):
    """Reads ERA5 wind data
    """

    def __call__(self, grid: Grid, start_time: str, end_time: str, expansion_factor: float):
        """Reads boundary spectra between given times and given area around
        the Grid object."""

        msg.info(
            f"Getting ERA5 wind forcing from {start_time} to {end_time}")


        temp_folder = 'dnora_wnd_temp'
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        msg.plain("Removing old files from temporary folder...")
        for f in glob.glob("dnora_wnd_temp/EC_ERA5.nc"):
            os.remove(f)

        # Define area to search in
        lon_min, lon_max, lat_min, lat_max = expand_area(min(grid.lon()), max(grid.lon()), min(grid.lat()), max(grid.lat()), expansion_factor)


        # start_times, end_times = create_monthly_time_stamps(start_time, end_time)
        # wnd_list = []
        # for t0, t1 in zip(start_times, end_times):
        #     msg.plain(f"Reading wind forcing data: {t0}-{t1}")
        #     # Creates file dnora_wnd_tmp/EC_ERA5_YYYY_MM.nc

        nc_file = download_era5_from_cds(start_time, end_time, lon=(lon_min, lon_max), lat=(lat_min, lat_max), folder='dnora_wnd_temp')
        wind_forcing = xr.open_dataset(nc_file)
        wind_forcing = wind_forcing.rename_dims({'longitude': 'lon', 'latitude': 'lat'})
        wind_forcing = wind_forcing.rename_vars({'longitude': 'lon', 'latitude': 'lat'})
        wind_forcing = wind_forcing.rename_vars({'u10': 'u', 'v10': 'v'})

        wind_forcing = wind_forcing.isel(lat=slice(None,None,-1)) # ERA5 gives lat as descending

        # Extract relevant timestamps
        wind_forcing = wind_forcing.sel(time=slice(start_time, end_time))

        return wind_forcing
