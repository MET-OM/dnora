from copy import copy
import numpy as np
import xarray as xr


# Import objects
from dnora.grid import Grid

# Import aux_funcsiliry functions
from dnora import msg
from dnora.aux_funcs import (
    create_time_stamps,
    expand_area,
)

from dnora.data_sources import DataSource
from dnora.readers.abstract_readers import DataReader


class GFS(DataReader):
    """Reads wind data of the GFS global forecast"""

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

        self.stride = copy(stride)
        self.hours_per_file = copy(hours_per_file)
        self.lead_time = copy(lead_time)
        self.last_file = copy(last_file)

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
            url = self.get_url(file_times[n], start_times[n], first_ind=self.lead_time)

            msg.from_file(url)
            msg.plain(f"Reading wind forcing data: {start_times[n]}-{end_times[n]}")

            import dask

            with dask.config.set(**{"array.slicing.split_large_chunks": True}):
                with xr.open_dataset(url) as ds:
                    ds["time"] = ds.time.dt.round("H")
                    ds = ds.sel(
                        time=slice(start_times[n], end_times[n]),
                        lon=slice(lon[0], lon[1]),
                        lat=slice(lat[0], lat[1]),
                    )[["lon", "lat", "time", "ugrd10m", "vgrd10m"]]

            wnd_list.append(ds)

        wind_forcing = xr.concat(wnd_list, dim="time")

        u = wind_forcing.ugrd10m.values
        v = wind_forcing.vgrd10m.values
        u = np.moveaxis(u, 0, 2)
        v = np.moveaxis(v, 0, 2)
        data_dict = {"u": u, "v": v}
        coord_dict = {
            "time": wind_forcing.time.values,
            "lon": wind_forcing.lon.values,
            "lat": wind_forcing.lat.values,
        }
        meta_dict = wind_forcing.attrs
        metaparameter_dict = {}

        return coord_dict, data_dict, meta_dict, metaparameter_dict

    def get_url(self, time_stamp_file, time_stamp, first_ind) -> str:
        h0 = int(time_stamp_file.hour)
        folder = "gfs" + time_stamp_file.strftime("%Y%m%d")
        filename = f"gfs_0p25_1hr_{h0:02.0f}z"

        return (
            "http://nomads.ncep.noaa.gov:80/dods/gfs_0p25_1hr/"
            + folder
            + "/"
            + filename
        )
