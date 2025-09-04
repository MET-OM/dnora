import xarray as xr

import os, glob

# Import objects
from dnora.grid import Grid


# Import aux_funcsiliry functions
from dnora import msg
from dnora import utils


import pandas as pd

from dnora.type_manager.data_sources import DataSource
from dnora.read.abstract_readers import DataReader

from dnora.cacher.caching_strategies import CachingStrategy
from dnora.read.ds_read_functions import setup_temp_dir

from dnora.type_manager.dnora_types import DnoraDataType

from dnora.read.depreciation_decorator import deprecated_class_call


def download_era5_from_cds(
    start_time, end_time, lon, lat, folder="dnora_wnd_temp"
) -> str:
    """Downloads ERA5 10 m wind data from the Copernicus Climate Data Store for a
    given area and time period"""
    start_time = pd.Timestamp(start_time)
    end_time = pd.Timestamp(end_time)
    try:
        import cdsapi
    except ImportError as e:
        msg.advice(
            "The cdsapi package is required to use ECWMF products! Install by e.g. 'conda install cdsapi'"
        )
        raise e
    c = cdsapi.Client()

    filename = f"{folder}/ERA5_cds_temp.nc"
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
    times = pd.date_range(start_time, end_time, freq="1d")
    # years = [f"{y:4.0f}" for y in utils.time.int_list_of_years(start_time, end_time)]
    years = list(set(times.strftime("%Y")))
    years.sort()
    months = list(set(times.strftime("%m")))
    months.sort()
    days = [f"{n:02.0f}" for n in range(1, 32)]

    cds_command = {
        "product_type": "reanalysis",
        "format": "netcdf",
        "variable": [
            "10m_u_component_of_wind",
            "10m_v_component_of_wind",
        ],
        "year": years,
        "month": months,
        "day": days,
        "time": [
            "00:00",
            "01:00",
            "02:00",
            "03:00",
            "04:00",
            "05:00",
            "06:00",
            "07:00",
            "08:00",
            "09:00",
            "10:00",
            "11:00",
            "12:00",
            "13:00",
            "14:00",
            "15:00",
            "16:00",
            "17:00",
            "18:00",
            "19:00",
            "20:00",
            "21:00",
            "22:00",
            "23:00",
        ],
        "area": [
            lat[1],
            lon[0],
            lat[0],
            # lat[1], 4, lat[0],
            lon[1],
        ],
    }

    c.retrieve("reanalysis-era5-single-levels", cds_command, filename)
    return filename


@deprecated_class_call("ERA5", "era5", "wind")
class ERA5(DataReader):
    """Reads ERA5 wind data"""

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
            grid.edges("lon"), grid.edges("lat"), expansion_factor, dlon=0.25, dlat=0.25
        )

        nc_file = download_era5_from_cds(
            start_time, end_time, lon=lon, lat=lat, folder=temp_folder
        )
        wind_forcing = xr.open_dataset(nc_file)

        wind_forcing = wind_forcing.isel(
            latitude=slice(None, None, -1)
        )  # ERA5 gives lat as descending

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
