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

from scipy.interpolate import griddata
# Import abstract classes
from .read import WaterLevelReader

# Import aux_funcsiliry functions
from .. import msg
from ..aux_funcs import create_time_stamps, u_v_from_dir, expand_area, lon_in_km, create_monthly_time_stamps, int_list_of_days, int_list_of_months, int_list_of_years
import pandas as pd
import matplotlib.pyplot as plt
from subprocess import Popen



def download_GTSM_from_cds(start_time, end_time, folder='dnora_wlv_temp') -> str:
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


    cds_command = {
        'format': 'tgz',
        'variable': ['total_water_level'],
        'experiment': 'reanalysis',
        'temporal_aggregation': 'hourly',
        'year': years, # 1979-2018
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

        out_file = download_GTSM_from_cds(start_time, end_time, folder='dnora_wlv_temp')

        temppath = os.path.dirname(out_file)
        # first unpack the tar.gz file.
        nc_file = subprocess.run(['tar', '-ztf', out_file], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')[0:-1]
        nc_file = sorted([ff.strip('\r') for ff in nc_file])
        #print(nc_file)
        subprocess.run(['tar', '-xzvf', out_file,'--directory',temppath], stdout=subprocess.PIPE) # Extract tar file

        lon_local = np.arange(lon_min, lon_max, 0.1)
        lat_local = np.arange(lat_min, lat_max, 0.1)
        grid_x, grid_y = np.meshgrid(lon_local, lat_local, indexing='xy')

        print(nc_file)
        grid_tot = []
        time_tot = []
        for ncfile in nc_file:
            #print(os.path.join(temppath,nc_file))
            waterlevel = xr.open_dataset(os.path.join(temppath,ncfile),engine="netcdf4")
            waterlevel = waterlevel.rename_vars({'station_x_coordinate': 'lon', 'station_y_coordinate': 'lat'})
            waterlevel = waterlevel.sel(stations=waterlevel.lon >= lon_min-expansion_factor)
            waterlevel = waterlevel.sel(stations=waterlevel.lon <= lon_max+expansion_factor)
            waterlevel = waterlevel.sel(stations=waterlevel.lat >= lat_min-expansion_factor)
            waterlevel = waterlevel.sel(stations=waterlevel.lat <= lat_max+expansion_factor)

            waterlevel = waterlevel.sel(time=slice(start_time,end_time))
            #grid_x, grid_y = np.mgrid[lon_min:lon_max:100j, lat_min:lat_max:110j]

            points = np.array([waterlevel.lon, waterlevel.lat]).T
            time = waterlevel.time
            grid_z = np.zeros([len(time),len(lat_local), len(lon_local)])
            for i_t, t in enumerate(time):
                values = waterlevel.waterlevel[i_t, :]
                grid_z[i_t,:,:] = griddata(points, values, (grid_x, grid_y), method='cubic', fill_value=0.)
                # plt.imshow(grid_z[i_t,:,:], extent=(lon_min, lon_max, lat_min, lat_max), origin='lower')
                # plt.scatter(waterlevel.lon, waterlevel.lat)
                # plt.colorbar()
                # plt.show()
            grid_tot.append(grid_z)
            time_tot.append(time)

        grid_tot = np.concatenate(grid_tot,axis=0)
        time_tot = np.concatenate(time_tot,axis=0)

        # Finally we put the new gridded data into a dataset
        waterlevel_gridded = xr.Dataset(
            data_vars=dict(
                waterlevel=(["time", "lat", "lon"], grid_tot)
            ),
            coords=dict(
                time=(["time"], time_tot),
                lat=(["lat"], lat_local),
                lon=(["lon"], lon_local),
            ),
            attrs=dict(description="waterlevel"),
        )

        #print(waterlevel_gridded)
        return waterlevel_gridded
