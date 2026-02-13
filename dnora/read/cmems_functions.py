import os
import pandas as pd
import xarray as xr
from dnora import msg

def ds_cmems_read(
    start_time: pd.Timestamp,
    end_time: pd.Timestamp,
    url: str,
    ## Partial variables from ProductReader
    lon: tuple[float],
    lat: tuple[float],
    ## Partial variables in ProductConfiguration
    **kwargs

):   
    
    cred_file = os.path.expanduser("~/.copernicusmarine/.copernicusmarine-credentials")
    try:
        import copernicusmarine
    except ImportError as e:
        msg.advice("The Copernicus Marine Service Toolbox is required to use ECWMF products! Install by e.g. 'conda install copernicusmarine'")
        raise e
    
    if not os.path.isfile(cred_file):
        msg.advice(
            f"No credentials file {cred_file} was found. Login for the first time to create it."
        )
    
        copernicusmarine.login()
    ds = copernicusmarine.open_dataset(
        dataset_id=kwargs.get('dataset_id'),
        variables=kwargs.get('variables'),
        minimum_longitude=lon[0],
        maximum_longitude=lon[1],
        minimum_latitude=lat[0],
        maximum_latitude=lat[1],
        start_datetime=start_time.strftime("%Y-%m-%dT%H:%M:00"),
        end_datetime=end_time.strftime("%Y-%m-%dT%H:%M:00"),
        minimum_depth=kwargs.get('minimum_depth',0),
        maximum_depth=kwargs.get('maximum_depth',0),

    )
    return ds
