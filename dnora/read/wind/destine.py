import pandas as pd
from dnora.type_manager.data_sources import DataSource

import xarray as xr
import numpy as np
from dnora import msg
from scipy.interpolate import griddata
from dnora.read.product_readers import ProductReader
from dnora.read.product_configuration import ProductConfiguration
from dnora.read.file_structure import FileStructure
from dnora import utils

def ds_polytope_read(
    start_time: pd.Timestamp,
    end_time: pd.Timestamp,
    url: str,
    lon: np.ndarray,
    lat: np.ndarray,
    dnora_class,
    **kwargs,
):

    grib_file = download_ecmwf_from_destine(
        start_time, end_time, lon=lon, lat=lat, folder='dnora_wind_temp'
    )
    
    ds = xr.open_dataset(grib_file, engine='cfgrib', decode_timedelta=True)
    lons, lats, u10, v10 = ds.u10.longitude.values, ds.u10.latitude.values, ds.u10.values, ds.v10.values
    native_dlon, native_dlat = 1/8, 1/30
    xi = np.arange(min(lons), max(lons), native_dlon)
    yi = np.arange(min(lats), max(lats), native_dlat)
    Xi, Yi = np.meshgrid(xi, yi)
    
    Nt = len(ds.step)
    u10i = np.zeros((Nt, len(yi), len(xi)))
    v10i = np.zeros((Nt, len(yi), len(xi)))
    # If this becomes slow, we need to think about 3D interpolation / resuing weights
    for n in range(Nt):
        u10i[n,:,:] = griddata(list(zip(lons, lats)), u10[n,:], (Xi, Yi), method='nearest')
        v10i[n,:,:] = griddata(list(zip(lons, lats)), v10[n,:], (Xi, Yi), method='nearest')
    
    data = dnora_class(lon=xi, lat=yi, time=ds.time+ds.step)
    data.set_u(u10i)
    data.set_v(v10i)
    lo, la = utils.grid.expand_area(lon, lat, expansion_factor=1, dlon=native_dlon, dlat=native_dlat)
    data = data.sel(lon=slice(*lo), lat=slice(*la))
    
    return data.sel(time=slice(start_time, end_time)).ds()

def download_ecmwf_from_destine(start_time, end_time, lon, lat, folder: str) -> str:
    """Downloads ERA5 10 m wind data from the Copernicus Climate Data Store for a
    given area and time period"""
    start_time = pd.Timestamp(start_time)
    end_time = pd.Timestamp(end_time)
    try:
        from polytope.api import Client
    except ImportError as e:
        msg.advice("The polytope package is required to acces these data! Install by e.g. 'python -m pip install polytope-client' and 'conda install cfgrib eccodes=2.41.0'")
        raise e
    c = Client(address='polytope.lumi.apps.dte.destination-earth.eu')

    filename = f"{folder}/ECMWF_temp.grib" # Switch to this in production. Then the files will be cleaned out
    #filename = f"{folder}/destine_temp.grib"
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
        "area":[int(np.ceil(lat[1])), int(np.floor(lon[0])), int(np.floor(lat[0])), int(np.ceil(lon[1]))],
    }
    date_str = start_time.strftime('%Y%m%d')
    request_winds['date'] = date_str
    
    c.retrieve('destination-earth', request_winds, filename)
    return filename

class ECMWF(ProductReader):
    """Downloads ECMWF data from Desinte using polytope api"""
    
    product_configuration = ProductConfiguration(
        ds_creator_function=ds_polytope_read,
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure(stride=24, hours_per_file=24)
