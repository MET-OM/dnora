# from ..bnd import Boundary
# from ..wnd import Forcing
# from ..spc import Spectra
# from ..wsr import WaveSeries
# from ..wlv import WaterLevel
# from ..ocr import OceanCurrent
from ..file_module import FileNames
import glob
from copy import copy
from .. import msg
import pandas as pd
import numpy as np
import re
from dnora.aux_funcs import day_list, expand_area
from dnora.dnora_types import DnoraDataType
from dnora.model_formats import ModelFormat
from dnora.readers.generic_readers import DataReader
from inspect import getcallargs
from dnora.grid import Grid

from geo_skeletons import GriddedSkeleton


def create_tiles(area, start_time, end_time, expansion_factor):
    tile_res = 5  # degrees
    lon, lat = expand_area(
        area.edges("lon", native=True), area.edges("lat", native=True), expansion_factor
    )
    lon, lat = np.array(lon), np.array(lat)
    days = day_list(start_time, end_time)

    lon[0] = np.floor(lon[0] / tile_res) * tile_res
    lon[1] = np.ceil(lon[1] / tile_res) * tile_res
    lat[0] = np.floor(lat[0] / tile_res) * tile_res
    lat[1] = np.ceil(lat[1] / tile_res) * tile_res

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


class Area(GriddedSkeleton):
    @classmethod
    def from_grid(cls, grid: Grid) -> "Area":
        return cls(
            lon=grid.edges("lon", strict=True),
            lat=grid.edges("lat", strict=True),
            x=grid.edges("x", strict=True),
            y=grid.edges("y", strict=True),
        )

    def expand(self, expansion_factor) -> "Area":
        lon, lat = expand_area(
            self.edges("lon", native=True),
            self.edges("lat", native=True),
            expansion_factor,
        )

        kwargs = {self.x_str: lon, self.y_str: lat}
        return Area(**kwargs)


class TileObject:
    def __init__(self, file_object: FileNames):
        self._file_object = file_object
        self._covering_files = []
        self.tile_res = 5

    def existing_files(self) -> list[str]:
        return glob.glob(
            f"{'_'.join(self._file_object.get_filepath().split('_')[:-3])}*"
        )

    def create_covering_files(
        self, area: Area, start_time, end_time, expansion_factor: float = 1.0
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

    def extent(self) -> tuple[tuple[float]]:
        lon, lat = self.lonlat(self.covering_files())
        lon_edges = (min([l[0] for l in lon]), max([l[1] for l in lon]))
        lat_edges = (min([l[0] for l in lat]), max([l[1] for l in lat]))
        return lon_edges, lat_edges

    def relevant_files(self) -> list[str]:
        """The files that cover the grid and exists in the cache"""
        if self.covering_files is None:
            return self.existing_files()

        return list(set(self.existing_files()).intersection(set(self.covering_files())))

    def additional_files(self) -> list[str]:
        """The files that cover the grid but doesn't exist in the cache"""
        if self.covering_files is None:
            return []
        return list(set(self.covering_files()) - set(self.existing_files()))

    @staticmethod
    def times(filenames: list[str]) -> tuple[list[int]]:
        year, month, day = [], [], []
        for filename in filenames:
            year.append(int(filename.split("_")[-1].split("-")[-3][0:4]))
            month.append(int(filename.split("_")[-1].split("-")[-2][0:2]))
            day.append(int(filename.split("_")[-1].split("-")[-1][0:2]))
        return year, month, day

    def lonlat(self, filenames: list[str]) -> tuple[list[tuple[int]]]:
        lon, lat = [], []
        for filename in filenames:

            lonlat = filename.split("_")[-3:-1]
            lonstr = [s for s in lonlat if "E" in s][0]
            latstr = [s for s in lonlat if "N" in s][0]

            lon.append((int(lonstr[0:3]), int(lonstr[0:3]) + self.tile_res))
            lat.append((int(latstr[0:3]), int(latstr[0:3]) + self.tile_res))
        return lon, lat

    def determine_patch_period(self, single_patch) -> list[tuple]:
        all_lon, all_lat = self.lonlat(self.covering_files())

        number_of_spatial_tiles = len(
            np.unique(
                [f"{la[0]:03.0f}N{lo[0]:03.0f}E" for lo, la in zip(all_lon, all_lat)]
            )
        )
        year, month, day = self.times(self.relevant_files())
        patch_dates = []
        for y, m, d in zip(year, month, day):
            tiles_for_day = sum(
                np.logical_and(
                    np.logical_and(np.array(year) == y, np.array(month) == m),
                    np.array(day) == d,
                )
            )
            # If even one tile is missing for that day, we will read the day for all area
            if tiles_for_day < number_of_spatial_tiles:
                patch_dates.append(f"{y:04.0f}-{m:02.0f}-{d:02.0f}")

        patch_dates = pd.to_datetime(patch_dates)
        if single_patch:
            return [(min(patch_dates), max(patch_dates) + pd.Timedelta(hours=23.99))]
        else:
            return [(start, start + pd.Timedelta(hours=23.99)) for start in patch_dates]


def cached_reader(obj_type: DnoraDataType, cache_reader: DataReader):
    def outer_wrapper(import_method):  # import_method is e.g. modelrun.import_spectra()
        def wrapper(
            *args,
            cache_name: str = None,
            read_cache: bool = False,
            write_cache: bool = False,
            patch: bool = False,
            **kwargs,
        ):
            def determine_patch_periods():
                """Determines if there is some periods that we need to patch from thredds
                adter reading cached data"""

                # This is not optimal, but seems to work
                times = mrun[obj_type].time()
                dt = times[1] - times[0]
                wanted_times = pd.date_range(
                    start=mrun.start_time(), end=mrun.end_time(), freq=dt
                )
                wanted_times.isin(times)

                if np.all(wanted_times.isin(times)):
                    return [], []
                wt = wanted_times.isin(times)

                was_found = "".join(
                    [str((w * 1)) for w in wt]
                )  # string of '0001110111110000'
                # was_found = '000011111111111100011111111111111111111111100011111111111111011111111111110000' # Testing
                inds = list(range(len(was_found)))
                was_found = re.sub("01", "0.1", was_found)
                was_found = re.sub("10", "1.0", was_found)

                list_of_blocks = was_found.split(".")

                patch_start = []
                patch_end = []
                for block in list_of_blocks:
                    if block[0] == "0":  # These need to be patched
                        ind_subset = inds[0 : len(block)]
                        patch_start.append(wanted_times[ind_subset[0]])
                        patch_end.append(wanted_times[ind_subset[-1]])
                    inds[0 : len(block)] = []

                return patch_start, patch_end

            # If no need for caching services, just execute method and return
            if not (read_cache or write_cache):
                import_method(*args, **kwargs)
                return

            # Get all arguments as keyword arguments
            kwargs = getcallargs(import_method, *args, **kwargs)
            mrun = kwargs.get("self")
            if kwargs.get("kwargs") is not None:
                conflicting_keys = list(
                    set(kwargs.keys()).intersection(set(kwargs.get("kwargs").keys()))
                )
                if conflicting_keys:
                    raise TypeError(
                        f"Got multiple inputs for key: {conflicting_keys[0]}!"
                    )
                kwargs.update(kwargs.get("kwargs"))
                del kwargs["kwargs"]

            # Don't care about caching in dry runs
            if kwargs.get("dry_run", False) or mrun.dry_run():
                import_method(**kwargs)
                return

            # This name will be used in the folders
            # Default use name of the reader, but user can override
            name = kwargs.get("name")
            if name is None:
                given_reader = kwargs.get("reader") or mrun._get_reader(obj_type)
                if given_reader is not None:
                    name = given_reader.name()
            if name is None:
                raise ValueError("Provide a DataReader!")

            ## Define file_obect to get proper file name
            file_object = FileNames(
                format=ModelFormat.CACHE,
                obj_type=obj_type,
                obj_name=name,
                model=mrun,
                edge_object=DnoraDataType.GRID,
                filename=cache_name,
            )

            tiles = TileObject(file_object)
            tiles.create_covering_files(
                Area.from_grid(mrun.grid()),
                mrun.start_time(),
                mrun.end_time(),
                kwargs.get("expansion_factor", 1.0),
            )

            if read_cache and tiles.relevant_files():
                kwargs_cache = copy(kwargs)
                kwargs_cache["reader"] = cache_reader(files=tiles.relevant_files())
                import_method(**kwargs_cache)
            else:
                import_method(**kwargs)

            if tiles.additional_files():
                patch_dates = tiles.determine_patch_period()
                mrun_temp = copy(mrun)
                lon, lat = tiles.extent()
                for start_time, end_time in patch_dates:
                    mrun_temp._time = pd.date_range(start_time, end_time, freq="h")
                    mrun_temp._grid = Grid(lon=lon, lat=lat)
                    mrun_temp._import_data(obj_type, **kwargs)
                    mrun[obj_type]._absorb_object(mrun_temp[obj_type], "time")

            # If write_cache = True, then the data is first read in tiles and saved
            # The data for the specific requested area is then read from these tiles
            if write_cache:
                for (lon, lat) in tiles.lonlat(tiles.covering_files()):




                existing_files = glob.glob(
                    f"{'_'.join(file_object.get_filepath().split('_')[:-3])}*"
                )

                relevant_files = list(
                    set(existing_files).intersection(set(covering_files))
                )

                file_object.create_folder()

                # if read_cache:  # Use existing tiles where possible
                #     files =
                # else:  # Download all needed tiles to write to cache
                #     files = all_files

                # Loop through each tile and cache them
                kwargs_cache = copy(kwargs)
                kwargs_cache["reader"] = given_reader
                kwargs_cache["name"] = name
                del kwargs_cache["self"]
                for n, _ in enumerate(relevant_files):
                    # Import method crops read area based on the Grid

                    # Do daily tiles in time
                    start_time = pd.to_datetime(times[n])
                    end_time = pd.to_datetime(times[n]) + pd.Timedelta(hours=24)

                    # This creates an DNORA object
                    obj = mrun._read_data(
                        obj_type=obj_type,
                        grid=Grid(lon=lons[n], lat=lats[n]),
                        start_time=start_time,
                        end_time=end_time,
                        **kwargs_cache,
                    )

                    # We don't want both edges in the tile so we don't get overlap
                    ind_lon = np.where(obj.lon() < lons[n][1])[0]
                    ind_lat = np.where(obj.lat() < lats[n][1])[0]
                    ind_time = np.where(obj.time() < end_time)[0]
                    obj = obj.isel(lon=ind_lon, lat=ind_lat, time=ind_time)
                    obj.name = name
                    mrun[obj_type] = obj
                    mrun.cache(obj_type)

                # Needs to be set, since we now want to read the data from the cahed tiles
                read_cache = True

            existing_files = glob.glob(
                f"{'_'.join(file_object.get_filepath().split('_')[:-3])}*"
            )
            relevant_files = list(set(existing_files).intersection(set(covering_files)))
            year, month, day, lon, lat = info_from_filename(relevant_files)
            breakpoint()
            if relevant_files and read_cache:
                kwargs_cache = copy(kwargs)
                kwargs_cache["reader"] = cache_reader(files=relevant_files)

                import_method(**kwargs_cache)
            else:
                import_method(**kwargs)

            # If cached files did not cover area, we need to patch
            # patching_times = [time for time in times if ]
            if patching_files:
                kwargs_cache = copy(kwargs)
                kwargs_cache["reader"] = given_reader
                kwargs_cache["name"] = name
                del kwargs_cache["self"]
                for n, _ in enumerate(relevant_files):
                    # Import method crops read area based on the Grid

                    # Do daily tiles in time
                    start_time = pd.to_datetime(times[n])
                    end_time = pd.to_datetime(times[n]) + pd.Timedelta(hours=24)

                    # This creates an DNORA object
                    obj = mrun._read_data(
                        obj_type=obj_type,
                        grid=Grid(lon=lons[n], lat=lats[n]),
                        start_time=start_time,
                        end_time=end_time,
                        **kwargs_cache,
                    )

                breakpoint()
            if files and read_cache:
                patch_start, patch_end = determine_patch_periods()
                if patch_start and patch:
                    msg.info(
                        "Not all data found in cache. Patching from original source..."
                    )
                    temp_args = tuple(list(args)[1:])
                    for t0, t1 in zip(patch_start, patch_end):
                        mrun_temp = copy(mrun)
                        mrun_temp.start_time = t0
                        mrun_temp.end_time = t1
                        exec(
                            f"mrun_temp.import_{obj_type.lower()}(*temp_args, **kwargs)"
                        )
                        mrun[obj_type]._absorb_object(mrun_temp[obj_type], "time")

            if write_cache:
                mrun.cache(obj_type)

        return wrapper

    return outer_wrapper
