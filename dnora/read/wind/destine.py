# Requirements
# pip install polytope-client


# import os

# from pathlib import Path
# import os
import pandas as pd
from dnora.grid import Grid
from dnora.type_manager.data_sources import DataSource
from dnora.read.abstract_readers import DataReader

from dnora.cacher.caching_strategies import CachingStrategy
from dnora.read.ds_read_functions import setup_temp_dir

from dnora.type_manager.dnora_types import DnoraDataType
from dnora import utils
import xarray as xr
import numpy as np
from dnora import msg
# c = Client(address='polytope.lumi.apps.dte.destination-earth.eu')

# datafolder = 'extremesdt_data'
# basepath = Path.cwd()
# data_path = Path.cwd() / Path(datafolder)

# if not os.path.isdir(data_path):
#         os.mkdir(data_path)

# request_waves = {
#     'class': 'd1',
#     'expver': '0001',
#     'dataset': 'extremes-dt',
#     'stream': 'wave',
#     'type': 'fc',
#     'levtype': 'sfc',
#     'param' : '140221/140229/140230/140231/140232',
#     'time': '00',
#     'step': '0/1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23',
# }



# date_str = '20250511'


def download_ecmwf_from_destine(
    start_time, end_time, lon, lat, folder="dnora_wnd_temp"
) -> str:
    """Downloads ERA5 10 m wind data from the Copernicus Climate Data Store for a
    given area and time period"""
    start_time = pd.Timestamp(start_time)
    end_time = pd.Timestamp(end_time)
    try:
        from polytope.api import Client
    except ImportError as e:
        msg.advice("The Ppolytope package is required to acces these data! Install by e.g. 'python -m pip install polytope-client'")
        raise e
    c = Client(address='polytope.lumi.apps.dte.destination-earth.eu')


    filename = f"{folder}/destine_temp.nc"

    times = pd.date_range(start_time, end_time, freq="1d")
    # years = [f"{y:4.0f}" for y in utils.time.int_list_of_years(start_time, end_time)]
    years = list(set(times.strftime("%Y")))
    years.sort()
    months = list(set(times.strftime("%m")))
    months.sort()
    days = [f"{n:02.0f}" for n in range(1, 32)]

    request_winds = {
        'class': 'd1',
        'expver': '0001',
        'dataset': 'extremes-dt',
        'stream': 'oper',
        'type': 'fc',
        'levtype': 'sfc',
        'param' : '165/166',
        'time': '00',
        'step': '0/1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23',
        #"area":[int(np.ceil(lat[1])), 4, 65, 40]
        "area":[int(np.ceil(lat[1])), int(np.floor(lon[0])), int(np.floor(lat[0])), int(np.ceil(lon[1]))],
        "year": years,
        "month": months,
        "day": days,
    }
    
    #date_str = '20250511'
    #request_winds['date'] = date_str
    breakpoint()
    c.retrieve('destination-earth', request_winds, filename)
    return filename

class ECMWF(DataReader):
    """Reads ECMWF wind data"""

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def caching_strategy(self) -> CachingStrategy:
        return CachingStrategy.SinglePatch

    def __call__(
        self,
        obj_type,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        expansion_factor: float = 1.2,
        **kwargs,
    ):
        """Reads boundary spectra between given times and given area around
        the Grid object."""

        msg.info(f"Getting ERA5 wind forcing from {start_time} to {end_time}")
        temp_folder = setup_temp_dir(DnoraDataType.WIND, self.name())

        # Define area to search in
        lon, lat = utils.grid.expand_area(
            grid.edges("lon"), grid.edges("lat"), expansion_factor)#, dlon=0.25, dlat=0.25)

        nc_file = download_ecmwf_from_destine(
            start_time, end_time, lon=lon, lat=lat, folder=temp_folder
        )
        
        ds = xr.open_dataset(nc_file, engine='cfgrib')

        # wind_forcing = wind_forcing.isel(
        #     latitude=slice(None, None, -1)
        # )  # ERA5 gives lat as descending

        wind_forcing = wind_forcing.sel(valid_time=slice(start_time, end_time))

        time = wind_forcing.valid_time.values

        coord_dict = {
            "lon": wind_forcing.longitude.values,
            "lat": wind_forcing.latitude.values,
            "time": time,
        }
        data_dict = {"u": wind_forcing.u10.data, "v": wind_forcing.v10.data}
        meta_dict = wind_forcing.attrs

        return coord_dict, data_dict, meta_dict


