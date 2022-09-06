from abc import ABC,  abstractmethod
from copy import copy
import numpy as np
import xarray as xr
import subprocess
#from subprocess import call, run
import os, glob
import time
import cdsapi
# Import objects
from ..grd.grd_mod import Grid

# Import abstract classes
from .read import WaterLevelReader

# Import aux_funcsiliry functions
from .. import msg
from ..aux_funcs import create_time_stamps, u_v_from_dir, expand_area, lon_in_km, create_monthly_time_stamps, int_list_of_days, int_list_of_months, int_list_of_years
import pandas as pd
from subprocess import Popen



def download_GTSM_from_cds(start_time, end_time, lon, lat, folder='dnora_wlv_temp') -> str:
    """Downloads GTSM model water level data from the Copernicus Climate Data Store for a
    given area and time period"""
    start_time = pd.Timestamp(start_time)
    end_time = pd.Timestamp(end_time)
    c = cdsapi.Client()

    filename = f'{folder}/EC_GTSM_ERA5.tar.gz'
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

    years = [f'{y:4.0f}' for y in int_list_of_years(start_time, end_time)]
    if len(years) == 1:
        years = years[0]
    months = [f'{m:02.0f}' for m in int_list_of_months(start_time, end_time)]
    if len(months) == 1:
        months = months[0]
    days = [f'{d:02.0f}' for d in int_list_of_days(start_time, end_time)]
    if len(days) == 1:
        days = days[0]


    cds_command = {
        'format': 'tgz',
        'variable': ['total_water_level'],
        'experiment': 'reanalysis',
        'temporal_aggregation': 'hourly',
        'year': years,
        'month': months,
    }


    c.retrieve('sis-water-level-change-timeseries-cmip6', cds_command, filename)
    return filename


class GTSM_ERA5(WaterLevelReader):
    """Reads GTSM_ERA5 waterlevel data
    """

    def __call__(self, grid: Grid, start_time: str, end_time: str, expansion_factor: float):
        """Reads hourly water level from GTSM_ERA5 database"""

        msg.info(
            f"Getting GTSM/ERA5 water level from {start_time} to {end_time}")


        temp_folder = 'dnora_wlv_temp'
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        msg.plain("Removing old files from temporary folder...")
        for f in glob.glob("dnora_wlv_temp/EC_GTSM_ERA5.tar.gz"):
            os.remove(f)

        # Define area to search in
        lon_min, lon_max, lat_min, lat_max = expand_area(min(grid.lon()), max(grid.lon()), min(grid.lat()), max(grid.lat()), expansion_factor)


        # start_times, end_times = create_monthly_time_stamps(start_time, end_time)
        # wnd_list = []
        # for t0, t1 in zip(start_times, end_times):
        #     msg.plain(f"Reading wind forcing data: {t0}-{t1}")
        #     # Creates file dnora_wnd_tmp/EC_ERA5_YYYY_MM.nc

        out_file = download_GTSM_from_cds(start_time, end_time, lon=(lon_min, lon_max), lat=(lat_min, lat_max), folder='dnora_wlv_temp')

        temppath = os.path.dirname(out_file)
        # first unpack the tar.gz file.
        nc_file = subprocess.run(['tar', '-ztf', out_file], stdout=subprocess.PIPE).stdout.decode('utf-8')
        print(nc_file)
        subprocess.run(['tar', '-xzvf', out_file,'--directory',temppath], stdout=subprocess.PIPE) # Extract tar file

        print(os.path.join(temppath,nc_file))
        waterlevel = xr.open_dataset(os.path.join(temppath,nc_file).strip("\n"),engine="netcdf4")
        waterlevel = waterlevel.rename_vars({'station_x_coordinate': 'lon', 'station_y_coordinate': 'lat'})
        df = waterlevel.to_dataframe()


        # todo: add a function which
        #
        # wind_forcing = wind_forcing.isel(lat=slice(None,None,-1)) # ERA5 gives lat as descending

        return waterlevel

    def strip away
