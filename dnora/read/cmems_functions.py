from dnora.read.ds_read_functions import  setup_temp_dir
from dnora.type_manager.dnora_types import DnoraDataType
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
    name: str,
    data_type: DnoraDataType,
    ## Partial variables in ProductConfiguration
    **kwargs

):   
    temp_dir = setup_temp_dir(data_type, name)

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

    copernicusmarine.subset(
        dataset_id=kwargs.get('dataset_id'),
        variables=kwargs.get('variables'),
        minimum_longitude=lon[0],
        maximum_longitude=lon[1],
        minimum_latitude=lat[0],
        maximum_latitude=lat[1],
        start_datetime=start_time.strftime("%Y-%m-%dT%H:%M:00"),
        end_datetime=end_time.strftime("%Y-%m-%dT%H:%M:00"),
        minimum_depth=kwargs.get('minimum_depth'),
        maximum_depth=kwargs.get('maximum_depth'),
        output_directory=temp_dir,
        credentials_file=cred_file,
        force_download=True,
        output_filename=f"{name}_CMEMS_temp.nc",
    )

    ds = xr.open_dataset(f"{temp_dir}/{name}_CMEMS_temp.nc")
    return ds