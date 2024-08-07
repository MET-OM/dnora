import numpy as np
import xarray as xr
import pandas as pd
from datetime import timedelta
import dask

# Import objects
from dnora.grid import Grid

# Import aux_funcsiliry functions
from dnora import msg
from dnora.aux_funcs import create_time_stamps, expand_area, get_url

from dnora.dnora_type_manager.data_sources import DataSource
from dnora.readers.abstract_readers import DataReader


class GFS(DataReader):
    """Reads wind data of the GFS global forecast"""

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def __init__(
        self,
        stride: int = 6,
        hours_per_file: int = 121,
        last_file: str = "",
        lead_time: int = 0,
    ):
        """The data is currently in hourly files. Do not change the default
        setting unless you have a good reason to do so.
        """

        self.stride = stride
        self.hours_per_file = hours_per_file
        self.lead_time = lead_time
        self.last_file = last_file

    def _folder_filename(
        self, source: DataSource, folder: str, filename: str
    ) -> tuple[str]:
        if source == DataSource.REMOTE:
            folder = "http://nomads.ncep.noaa.gov:80/dods/gfs_0p25_1hr/gfs%Y%m%d"
        if filename is None:
            filename = "gfs_0p25_1hr_%Hz"
        return folder, filename

    def __call__(
        self,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        folder: str,
        expansion_factor: float = 1.2,
        **kwargs,
    ):
        """Reads wind data between given times and given area around
        the Grid object."""

        self.start_time = start_time
        self.end_time = end_time

        start_times, end_times, file_times = create_time_stamps(
            start_time,
            end_time,
            self.stride,
            self.hours_per_file,
            self.last_file,
            self.lead_time,
        )

        msg.info(
            f"Getting wind forcing from GFS from {self.start_time} to {self.end_time}"
        )
        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        # Define area to search in
        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

        wnd_list = []
        for n in range(len(file_times)):
            folder, filename = self._folder_filename(source, folder, filename=None)
            url = get_url(folder, filename, file_times[n])
            msg.plain(f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")
            msg.from_file(url)

            # with dask.config.set(**{"array.slicing.split_large_chunks": True}):
            #     with xr.open_dataset(url, decode_times=False) as ds:

            ds = xr.open_dataset(url, decode_times=False)

            t0 = pd.to_datetime(
                ds.time.minimum[-4:] + " " + ds.time.minimum[3:8]
            ) + timedelta(hours=int(ds.time.minimum[:2]))
            t1 = pd.to_datetime(
                ds.time.maximum[-4:] + " " + ds.time.maximum[3:8]
            ) + timedelta(hours=int(ds.time.maximum[:2]))

            # ds["time"] = ds.time.dt.round("H")
            ds["time"] = pd.date_range(t0, t1, freq="1h")

            ds = ds.sel(
                time=slice(start_times[n], end_times[n]),
                lon=slice(lon[0], lon[1]),
                lat=slice(lat[0], lat[1]),
            )[["lon", "lat", "time", "ugrd10m", "vgrd10m"]]
            wnd_list.append(ds)

        wind_forcing = xr.concat(wnd_list, dim="time")

        u = wind_forcing.ugrd10m.data
        v = wind_forcing.vgrd10m.data
        # u = np.moveaxis(u, 0, 2)
        # v = np.moveaxis(v, 0, 2)
        data_dict = {"u": u, "v": v}
        coord_dict = {
            "time": wind_forcing.time.data,
            "lon": wind_forcing.lon.data,
            "lat": wind_forcing.lat.data,
        }
        meta_dict = wind_forcing.attrs

        return coord_dict, data_dict, meta_dict
