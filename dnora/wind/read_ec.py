import xarray as xr

import os, glob
import cdsapi

# Import objects
from dnora.grid import Grid


# Import aux_funcsiliry functions
from dnora import msg
from dnora.aux_funcs import (
    expand_area,
    int_list_of_days,
    int_list_of_months,
    int_list_of_years,
)
import pandas as pd

from dnora.dnora_type_manager.data_sources import DataSource
from dnora.readers.abstract_readers import DataReader

from dnora.cacher.caching_strategies import CachingStrategy


def download_era5_from_cds(
    start_time, end_time, lon, lat, folder="dnora_wnd_temp"
) -> str:
    """Downloads ERA5 10 m wind data from the Copernicus Climate Data Store for a
    given area and time period"""
    start_time = pd.Timestamp(start_time)
    end_time = pd.Timestamp(end_time)
    c = cdsapi.Client()

    filename = f"{folder}/EC_ERA5.nc"
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

    years = [f"{y:4.0f}" for y in int_list_of_years(start_time, end_time)]
    if len(years) == 1:
        years = years[0]
    months = [f"{m:02.0f}" for m in int_list_of_months(start_time, end_time)]
    if len(months) == 1:
        months = months[0]
    days = [f"{d:02.0f}" for d in int_list_of_days(start_time, end_time)]
    if len(days) == 1:
        days = days[0]

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


class ERA5(DataReader):
    """Reads ERA5 wind data"""

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def _caching_strategy(self) -> CachingStrategy:
        return CachingStrategy.SinglePatch

    def __call__(
        self,
        grid: Grid,
        start_time: str,
        end_time: str,
        expansion_factor: float = 1.2,
        **kwargs,
    ):
        """Reads boundary spectra between given times and given area around
        the Grid object."""

        msg.info(f"Getting ERA5 wind forcing from {start_time} to {end_time}")

        temp_folder = "dnora_wnd_temp"
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        msg.plain("Removing old files from temporary folder...")
        for f in glob.glob("dnora_wnd_temp/EC_ERA5.nc"):
            os.remove(f)

        # Define area to search in
        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

        nc_file = download_era5_from_cds(
            start_time, end_time, lon=lon, lat=lat, folder="dnora_wnd_temp"
        )
        wind_forcing = xr.open_dataset(nc_file)
        # wind_forcing = wind_forcing.rename_dims({"longitude": "lon", "latitude": "lat"})
        # wind_forcing = wind_forcing.rename_vars({"longitude": "lon", "latitude": "lat"})
        # wind_forcing = wind_forcing.rename_vars({"u10": "u", "v10": "v"})

        wind_forcing = wind_forcing.isel(
            latitude=slice(None, None, -1)
        )  # ERA5 gives lat as descending

        coord_dict = {
            "lon": wind_forcing.longitude.values,
            "lat": wind_forcing.latitude.values,
            "time": wind_forcing.time.values,
        }
        data_dict = {"u": wind_forcing.u10.values, "v": wind_forcing.v10.values}
        meta_dict = wind_forcing.attrs

        return coord_dict, data_dict, meta_dict
