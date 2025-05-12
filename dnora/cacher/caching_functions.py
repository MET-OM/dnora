import numpy as np
from dnora.export import Cacher
from dnora.type_manager.data_sources import DataSource
from copy import copy
from dnora.grid import Grid
from .caching_strategies import caching_strategies, CachingStrategy
from dnora import msg
import pandas as pd
from dnora.type_manager.spectral_conventions import SpectralConvention
from dnora.type_manager.dnora_types import DnoraDataType


def dont_proceed_with_caching(read_cache, write_cache, strategy, kwargs):
    """Checks if there is any reason not to proceed with the cahcing process"""
    dont_proceed = False
    dont_proceed = dont_proceed or (not (read_cache or write_cache))
    dont_proceed = dont_proceed or strategy == CachingStrategy.DontCacheMe
    dont_proceed = dont_proceed or (
        kwargs.get("dry_run", False) or kwargs.get("self").dry_run()
    )
    return dont_proceed


def expand_area_to_tiles(tiles, dlon, dlat):
    """Expands in time and space to cover full daily tiles"""
    lon, lat = tiles.spatial_extent(tiles.covering_files())
    grid = Grid(lon=lon, lat=lat)
    grid.set_spacing(dlon=dlon, dlat=dlat)
    # Gives full days
    start_time, end_time = tiles.temporal_extent(tiles.covering_files())

    return grid, start_time, end_time


def read_data_from_cache(mrun_cacher, tiles, cache_reader, kwargs_cache):
    """Read all possible data from cached files"""
    if tiles.relevant_files():
        kwargs_read_cache = copy(kwargs_cache)
        kwargs_read_cache["reader"] = cache_reader(
            files=tiles.relevant_files(), used_for_cache=True
        )
        kwargs_read_cache["source"] = DataSource.LOCAL
        mrun_cacher._import_data(**kwargs_read_cache)
    return mrun_cacher


def patch_cached_data(mrun_cacher, tiles, kwargs_cache, strategy: CachingStrategy):
    """Patch data not found in the cached files from the original source"""

    strategy_func = caching_strategies.get(strategy)
    if strategy_func is None:
        msg.info(
            f"Caching strategy {strategy.name} not implemented! Reverting to SinglePatch."
        )
        strategy_func = caching_strategies.get(CachingStrategy.SinglePatch)
    patch_dates, patch_lon, patch_lat, patch_dimension = strategy_func(tiles)
    obj_type = kwargs_cache.get("obj_type")

    for patch_date, lon, lat in zip(patch_dates, patch_lon, patch_lat):
        grid_lon = (
            max(lon[0], mrun_cacher.grid().edges("lon")[0]),
            min(lon[1], mrun_cacher.grid().edges("lon")[1]),
        )
        grid_lat = (
            max(lat[0], mrun_cacher.grid().edges("lat")[0]),
            min(lat[1], mrun_cacher.grid().edges("lat")[1]),
        )
        grid = Grid(lon=grid_lon, lat=grid_lat)
        grid.set_spacing(nx=mrun_cacher.grid().nx(), ny=mrun_cacher.grid().ny())
        mrun_patch = mrun_cacher.empty_copy(
            grid=grid,
            start_time=patch_date[0],
            end_time=patch_date[1],
        )

        mrun_patch._import_data(**kwargs_cache, coverage_warning=False)

        ## Merge patch together with what was found in the cached
        if patch_dimension == "time":
            times_to_keep = np.logical_and(
                mrun_cacher[obj_type].time() < patch_date[0],
                mrun_cacher[obj_type].time() > patch_date[-1],
            )
        if mrun_cacher.get(obj_type) is None or not np.any(times_to_keep):
            mrun_cacher[obj_type] = mrun_patch[obj_type]
        else:
            times_to_keep = mrun_cacher[obj_type].time()[times_to_keep]
            mrun_cacher[obj_type] = mrun_cacher[obj_type].sel(times=times_to_keep)
            mrun_cacher[obj_type].absorb(mrun_patch[obj_type], patch_dimension)

    return mrun_cacher


def write_data_to_cache(mrun_cacher, tiles, obj_type):
    # Write spatial tile for spatial tile
    lons, lats = tiles.lonlat(tiles.covering_files())
    # Get unique values
    lons = list(set(lons))
    lats = list(set(lats))
    first_round = True
    for lon_tuple in lons:
        for lat_tuple in lats:
            mrun_write_tile = mrun_cacher.empty_copy(
                grid=Grid(lon=lon_tuple, lat=lat_tuple),
                start_time=mrun_cacher.start_time(),
                end_time=mrun_cacher.end_time(),
            )
            cropped_obj = mrun_cacher[obj_type]

            lon_mask = np.logical_and(
                cropped_obj.lon() < lon_tuple[1],
                cropped_obj.lon() >= lon_tuple[0],
            )
            lat_mask = np.logical_and(
                cropped_obj.lat() < lat_tuple[1],
                cropped_obj.lat() >= lat_tuple[0],
            )
            ind_lon = np.where(lon_mask)[0]
            ind_lat = np.where(lat_mask)[0]
            if cropped_obj.is_gridded():
                cropped_obj = cropped_obj.isel(lon=ind_lon, lat=ind_lat)
            else:
                sel_inds = np.array(list(set(ind_lon) & set(ind_lat)))
                if len(sel_inds) == 0:
                    coord_dict = cropped_obj.coord_dict()
                    coord_dict[cropped_obj.core.x_str] = lon_tuple
                    coord_dict[cropped_obj.core.y_str] = lat_tuple

                    cropped_obj = cropped_obj.__class__(**coord_dict)

                else:

                    cropped_obj = cropped_obj.isel(inds=sel_inds)

            if hasattr(mrun_cacher[obj_type], "convention"):
                cropped_obj._mark_convention(
                    mrun_cacher[obj_type].convention(), silent=True
                )
                cropped_obj.set_convention(SpectralConvention.OCEAN)
                if first_round:
                    msg.plain(f"Writing data in convention: {cropped_obj.convention()}")
                first_round = False
            cropped_obj.name = mrun_cacher[obj_type].name

            mrun_write_tile[obj_type] = cropped_obj
            exporter = Cacher(mrun_write_tile)  # Writes daily files
            exporter.export(obj_type, edge_object=DnoraDataType.GRID)


# def write_data_to_cache(mrun_cacher, tiles, obj_type):
#     # Write spatial tile for spatial tile
#     lons, lats = tiles.lonlat(tiles.covering_files())
#     years, months, days = tiles.times(tiles.covering_files())
#     for lon_tuple, lat_tuple, year, month, day in zip(lons, lats, years, months, days):

#         mrun_write_tile = mrun_cacher.empty_copy(
#             grid=Grid(lon=lon_tuple, lat=lat_tuple),
#             start_time=mrun_cacher.start_time(),
#             end_time=mrun_cacher.end_time(),
#         )
#         cropped_obj = mrun_cacher[obj_type]
#         lon_mask = np.logical_and(
#             cropped_obj.lon() < lon_tuple[1],
#             cropped_obj.lon() >= lon_tuple[0],
#         )
#         lat_mask = np.logical_and(
#             cropped_obj.lat() < lat_tuple[1],
#             cropped_obj.lat() >= lat_tuple[0],
#         )
#         ind_lon = np.where(lon_mask)[0]
#         ind_lat = np.where(lat_mask)[0]
#         if cropped_obj.is_gridded():
#             cropped_obj = cropped_obj.isel(lon=ind_lon, lat=ind_lat).sel(
#                 time=f"{year:.0f}-{month:02.0f}-{day:02.0f}"
#             )
#         else:
#             sel_inds = np.array(list(set(ind_lon) & set(ind_lat)))
#             cropped_obj = cropped_obj.isel(inds=sel_inds).sel(
#                 time=f"{year:.0f}-{month:02.0f}-{day:02.0f}"
#             )
#         cropped_obj.name = mrun_cacher[obj_type].name
#         if hasattr(mrun_cacher[obj_type], "convention"):
#             cropped_obj._mark_convention(mrun_cacher[obj_type].convention())
#         mrun_write_tile[obj_type] = cropped_obj
#         exporter = Cacher(mrun_write_tile)  # Writes daily files
#         exporter.export(obj_type)
