from dnora.read.abstract_readers import DataReader
from dnora.type_manager.data_sources import DataSource
from dnora.type_manager.dnora_types import DnoraDataType
from dnora.utils.io import get_url
from dnora.grid import Grid
from dnora import msg
from dnora import utils
import pandas as pd
import numpy as np


class SWAN_Ascii(DataReader):
    _default_data_source = DataSource.LOCAL

    @staticmethod
    def returning_ds():
        return True

    def __init__(self, lon: tuple[float] = None, lat: tuple[float] = None):
        self._lon = lon
        self._lat = lat

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        folder: str,
        filename: str,
        dnora_class,
        expansion_factor: float = 1.2,
        **kwargs,
    ):
        folder = self._folder(folder, source)
        filename = self._filename(filename, source)
        url = get_url(folder, filename)
        msg.info(
            f"Getting {obj_type} from {self.name()} from {start_time} to {end_time}"
        )

        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        lon, lat = utils.grid.expand_area(
            grid.edges("lon"), grid.edges("lat"), expansion_factor
        )
        times = []
        u_list = []
        v_list = []
        with open(url, "r") as file:
            line = True
            reading_u = False
            while line:
                line = file.readline()
                if is_date(line):
                    reading_u = not reading_u
                    times.append(pd.to_datetime("".join(line.split("."))))
                    continue

                if line:
                    if reading_u:
                        u_list.append(np.array(line.split(" ")).astype(float) / 1000)
                    else:
                        v_list.append(np.array(line.split(" ")).astype(float) / 1000)

        times = times[:-1:2]
        nt, ny, nx = len(times), int(len(u_list) / len(times)), len(u_list[0])
        us = np.zeros((nt, ny, nx))
        vs = np.zeros((nt, ny, nx))
        msg.plain(f"Identifying shape: nt={nt}, ny={ny}, nx={nx}")
        for ct in range(len(u_list)):
            ct_lat = np.mod(ct, ny)
            ct_time = np.floor(ct / ny).astype(int)
            # msg.plain(f"ct: {ct} and ct_lat = {ct_lat}, ct_time={ct_time}")
            us[ct_time, ct_lat, :] = u_list[ct]
            vs[ct_time, ct_lat, :] = v_list[ct]

        lon_edges = self._lon or grid.edges("lon")
        lat_edges = self._lat or grid.edges("lat")
        lons = np.linspace(*lon_edges, nx)
        lats = np.linspace(*lat_edges, ny)
        msg.plain(f"Time ({nt} steps): {times[0]}-{times[-1]}")
        msg.plain(f"Lon: {lons[0]}-{lons[-1]}, dlon={(lons[-1]-lons[0])/(nx-1):.7f}")
        msg.plain(f"Lat: {lats[0]}-{lats[-1]}, dlat={(lats[-1]-lats[0])/(ny-1):.7f}")

        data = dnora_class(time=times, lon=lons, lat=lats)
        data.set_u(us)
        data.set_v(vs)

        data = data.sel(time=slice(start_time, end_time))
        data = data.sel(lon=slice(*lon), lat=slice(*lat))
        data.meta.set(
            {
                "source": f"Read from {filename} assuming lon, lat coverage ({lon_edges[0]}, {lon_edges[-1]}), ({lat_edges[0]}, {lat_edges[-1]})"
            }
        )

        return data.ds()


def is_date(line: str):
    """Checks if a line is a new date in a SWAN ascii file"""
    if len(line) != 16:
        return False
    split_line = line.split(".")
    if len(split_line) != 2:
        return False

    day, hour = tuple(split_line)
    if len(day) != 8:
        return False
    if len(hour) != 7:  # 7th is end-line
        return False

    if int(day[:4]) > 2500:
        return False
    if int(day[4:6]) > 12:
        return False
    if int(day[6:8]) > 31:
        return False

    if int(hour[:2]) > 23:
        return False
    if int(hour[2:4]) > 59:
        return False

    return True
