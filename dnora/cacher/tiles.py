from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from dnora.file_module import FileNames

import pandas as pd
import glob
import numpy as np
from dnora import utils

from dnora.grid import Grid


class TileObject:
    def __init__(self, file_object: FileNames):
        self._file_object = file_object
        self._covering_files = []
        self.tile_res = 5

    def existing_files(self) -> list[str]:
        years, months, __ = self.times(self.covering_files())
        times = list(
            set([f"{year:.0f}-{month:02.0f}" for year, month in zip(years, months)])
        )
        times.sort()
        files = []
        for time in times:
            files += glob.glob(
                f"{'_'.join(self._file_object.get_filepath(start_time=time).split('_')[:-3])}*"
            )
        return files

    def create_covering_files(
        self, area: Grid, start_time, end_time, expansion_factor: float = 1.0
    ) -> None:
        lons, lats, times = create_tiles(
            area,
            start_time,
            end_time,
            expansion_factor,
        )
        covering_files = []
        for lon, lat, time in zip(lons, lats, times):
            covering_files.append(
                self._file_object.get_filepath(lon=lon, lat=lat, start_time=time)
            )
        self._covering_files = covering_files

    def covering_files(self) -> list[str]:
        return self._covering_files

    def spatial_extent(self, filenames: list[str]) -> tuple[tuple[float]]:
        lon, lat = self.lonlat(filenames)
        lon_edges = (min([l[0] for l in lon]), max([l[1] for l in lon]))
        lat_edges = (min([l[0] for l in lat]), max([l[1] for l in lat]))
        return lon_edges, lat_edges

    def temporal_extent(self, filenames: list[str]) -> tuple:
        times = pd.to_datetime(self.time_strings(filenames))
        return (min(times), max(times) + pd.Timedelta(hours=23.99))

    def relevant_files(self) -> list[str]:
        """The files that cover the grid and exists in the cache"""
        # if not self._covering_files:
        #     return self.existing_files()

        return list(set(self.existing_files()).intersection(set(self.covering_files())))

    def additional_files(self) -> list[str]:
        """The files that cover the grid but doesn't exist in the cache"""
        # if self.covering_files is None:
        #     return []
        return list(set(self.covering_files()) - set(self.existing_files()))

    @staticmethod
    def times(filenames: list[str]) -> tuple[list[int]]:
        year, month, day = [], [], []
        for filename in filenames:
            year.append(int(filename.split("_")[-1].split("-")[-3][0:4]))
            month.append(int(filename.split("_")[-1].split("-")[-2][0:2]))
            day.append(int(filename.split("_")[-1].split("-")[-1][0:2]))
        return year, month, day

    def time_strings(self, filenames: list[str]) -> list[str]:
        year, month, day = self.times(filenames)
        return [f"{y:04.0f}-{m:02.0f}-{d:02.0f}" for (y, m, d) in zip(year, month, day)]

    def lonlat_strings(self, filenames: list[str]) -> list[str]:
        lonlat_strs = []
        for filename in filenames:

            lonlat = filename.split("_")[-3:-1]
            lonstr = [s for s in lonlat if "E" in s][0]
            latstr = [s for s in lonlat if "N" in s][0]
            lonlat_strs.append(f"{latstr}{lonstr}")

        return lonlat_strs
        #return [
        #    f"{lat[0]:03.0f}N{lon[0]:03.0f}E" for (lat, lon) in self.latlon(filenames)
        #]

    def lonlat(self, filenames: list[str]) -> tuple[list[tuple[int]]]:
        """Creates a list of tuples of lons and lats that correspond to the given filenames"""
        lon, lat = [], []
        for filename in filenames:

            lonlat = filename.split("_")[-3:-1]
            lonstr = [s for s in lonlat if "E" in s][0]
            latstr = [s for s in lonlat if "N" in s][0]
            lon.append((int(lonstr[:-1]), int(lonstr[:-1]) + self.tile_res))
            lat.append((int(latstr[:-1]), int(latstr[:-1]) + self.tile_res))

        return lon, lat


def create_tiles(area, start_time, end_time, expansion_factor) -> tuple:
    """Creates tiles that cover 5x5 degrees in space and one day in time.
    Returns:

    lons/lats that is a list of tuple [(0.0, 5.0), (5.0,10.0), ...]
    times: list of strings ['2020-01-01', '2020-01-01', '2020-01-02', ...]
    """
    tile_res = 5  # degrees
    lon, lat = utils.grid.expand_area(
        area.edges("lon", native=True), area.edges("lat", native=True), expansion_factor
    )
    lon, lat = np.array(lon), np.array(lat)
    days = utils.time.day_list(start_time, end_time)

    lon[0] = np.floor(lon[0] / tile_res) * tile_res
    lon[1] = np.floor(lon[1] / tile_res) * tile_res + tile_res
    lat[0] = np.floor(lat[0] / tile_res) * tile_res
    lat[1] = np.floor(lat[1] / tile_res) * tile_res + tile_res
    lon_vec = np.arange(lon[0], lon[1] - tile_res / 2, tile_res)
    lat_vec = np.arange(lat[0], lat[1] - tile_res / 2, tile_res)
    lons = []
    lats = []
    times = []
    for d in days:
        for la in lat_vec:
            for lo in lon_vec:
                lons.append((lo, lo + tile_res))
                lats.append((la, la + tile_res))
                times.append(d)
    return lons, lats, times


def find_tile(lon, lat, tile_lon, tile_lat) -> int:
    west_of = lon < tile_lon[:, 1]
    east_of = lon >= tile_lon[:, 0]
    north_of = lat >= tile_lat[:, 0]
    south_of = lat < tile_lat[:, 1]
    lon_ok = np.logical_and(east_of, west_of)
    lat_ok = np.logical_and(north_of, south_of)
    if not np.any(np.logical_and(lon_ok, lat_ok)):
        return []

    ind = np.where(np.logical_and(lon_ok, lat_ok))[0][0]
    return ind


def info_from_filename(filenames: list[str]) -> tuple[int]:
    year, month, day = [], [], []
    lon, lat = [], []
    for filename in filenames:
        year.append(int(filename.split("_")[-1].split("-")[-3][0:4]))
        month.append(int(filename.split("_")[-1].split("-")[-2][0:2]))
        day.append(int(filename.split("_")[-1].split("-")[-1][0:2]))

        lonlat = filename.split("_")[-3:-1]
        lonstr = [s for s in lonlat if "E" in s][0]
        latstr = [s for s in lonlat if "N" in s][0]

        lon.append(int(lonstr[0:3]))
        lat.append(int(latstr[0:3]))
    return year, month, day, lon, lat
