from abc import ABC,  abstractmethod
import pandas as pd
# Import objects
from ..grd.grd_mod import Grid
from .. import msg
from ..aux import expand_area
import xarray as xr

class ForcingReader(ABC):
    """Reads forcing data from some source and provide it to the object.

    The area is defined from the Grid object that is passed.
    """

    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, grid: Grid, start_time: str, end_time: str, expansion_factor: float):
        pass

class NEMONc(ForcingReader):
    """Reads wind data from NEMO style netcdf files
    """

    def __call__(self, grid: Grid, start_time: str, end_time: str, expansion_factor: float):
        """Reads boundary spectra between given times and given area around
        the Grid object."""


        #start_times, end_times, file_times = monthly_stamps(start_time, end_time)

        msg.info(
            f"Getting wind forcing from {start_time} to {end_time}")


            # Define area to search in
        lon_min, lon_max, lat_min, lat_max = expand_area(min(grid.lon()), max(grid.lon()), min(grid.lat()), max(grid.lat()), expansion_factor)
        breakpoint()
        wind_forcing = xr.open_mfdataset('ERA5/*.nc', drop_variables=['d2m','msl','tcc','t2m','sf','ssr','tp','strd']).sel(time=slice(start_time, end_time), longitude=slice(lon_min, lon_max), latitude=slice(lat_max, lat_min))


        wind_forcing = wind_forcing.rename_dims({'latitude': 'lat', 'longitude': 'lon'})
        wind_forcing = wind_forcing.rename_vars({'u10': 'u', 'v10': 'v'})
        return wind_forcing

    def get_url(self, time_stamp_file) -> str:
        url = f"meteo.ERA5.baltic.1h.{time_stamp_file.strftime('%Y')}.{time_stamp_file.strftime('%m')}"
        return url
