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
from dnora.dnora_type_manager.dnora_types import DnoraDataType
from dnora.model_formats import ModelFormat
from dnora.readers.generic_readers import DataReader
from inspect import getcallargs
from dnora.grid import Grid
from dnora.pick import Area
from geo_skeletons import GriddedSkeleton

from .tiles import TileObject
from dnora.importer import DataImporter
from dnora.dnora_type_manager.data_sources import DataSource
from dnora.export import Cacher


def get_kwargs(func, args, kwargs) -> dict:

    kwargs = getcallargs(func, *args, **kwargs)

    if kwargs.get("kwargs") is not None:
        conflicting_keys = list(
            set(kwargs.keys()).intersection(set(kwargs.get("kwargs").keys()))
        )
        if conflicting_keys:
            raise TypeError(f"Got multiple inputs for key: {conflicting_keys[0]}!")
        kwargs.update(kwargs.get("kwargs"))
        del kwargs["kwargs"]

    return kwargs


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
            # If no need for caching services, just execute method and return
            if not (read_cache or write_cache):
                import_method(*args, **kwargs)
                return

            # Get all arguments as keyword arguments
            kwargs = get_kwargs(import_method, args, kwargs)
            mrun = kwargs.get("self")

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

            ## Covering tiles are daily 5degx5deg tiles that cover the (possibly expanded) area we want to read
            tiles = TileObject(file_object)
            tiles.create_covering_files(
                mrun.grid(),
                mrun.start_time(),
                mrun.end_time(),
                kwargs.get("expansion_factor", 1.0),
            )

            kwargs_cache = copy(kwargs)
            del kwargs_cache["self"]
            kwargs_cache["obj_type"] = obj_type
            if write_cache:
                ## We should always write full tiles to cache, so we need to read an area covering the tiles
                lon, lat = tiles.spatial_extent(tiles.covering_files())
                grid = Grid(lon=lon, lat=lat)
                grid.set_spacing(dlon=mrun.grid().dlon(), dlat=mrun.grid().dlat())
                kwargs_cache["point_picker"] = Area()
                kwargs_cache["expansion_factor"] = 1.0
                # Gives full days
                start_time, end_time = tiles.temporal_extent(tiles.covering_files())
            else:
                start_time, end_time = mrun.start_time(), mrun.end_time()
                grid = mrun.grid()
            ## If we have some files cached, read them in now
            mrun_in_cache = mrun.empty_copy(
                grid=grid,
                start_time=start_time,
                end_time=end_time,
            )
            if read_cache and tiles.relevant_files():
                kwargs_read_cache = copy(kwargs_cache)
                kwargs_read_cache["reader"] = cache_reader(files=tiles.relevant_files())
                kwargs_read_cache["source"] = DataSource.LOCAL
                mrun_in_cache._import_data(**kwargs_read_cache)

            ## Read everything regardless of what is cached
            mrun_patch = None
            if not read_cache:
                mrun_patch = mrun.empty_copy(
                    grid=grid,
                    start_time=start_time,
                    end_time=end_time,
                )

                mrun_patch._import_data(**kwargs_cache)

            ## Patch the data that was not in cache
            if tiles.additional_files() and read_cache:
                patch_dates = tiles.determine_patch_period(single_patch=True)
                mrun_patch = mrun.empty_copy(
                    grid=grid,
                    start_time=patch_dates[0][0],
                    end_time=patch_dates[0][1],
                )

                mrun_patch._import_data(**kwargs_cache)

            ## Merge patch together with what was found in the cached
            if mrun_in_cache[obj_type] is None:
                mrun_in_cache[obj_type] = mrun_patch[obj_type]
            elif mrun_patch is not None:
                mrun_in_cache[obj_type].absorb(mrun_patch[obj_type], "time")

            ## Write the data to cache if that is requested
            if write_cache:
                # Write spatial tile for spatial tile
                lons, lats = tiles.lonlat(tiles.covering_files())
                for lon_tuple, lat_tuple in zip(lons, lats):
                    mrun_write_tile = mrun.empty_copy(
                        grid=Grid(lon=lon_tuple, lat=lat_tuple),
                        start_time=start_time,
                        end_time=end_time,
                    )
                    cropped_obj = mrun_in_cache[obj_type]
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
                    cropped_obj = cropped_obj.isel(lon=ind_lon, lat=ind_lat)
                    cropped_obj.name = mrun_in_cache[obj_type].name
                    mrun_write_tile[obj_type] = cropped_obj
                    exporter = Cacher(mrun_write_tile)  # Writes daily files
                    exporter.export(obj_type)

            final_object = mrun_in_cache[obj_type]
            lon, lat = expand_area(
                grid.edges("lon", strict=True),
                grid.edges("lat", strict=True),
                expansion_factor=kwargs.get("expansion_factor", 1.0),
            )
            final_object = final_object.sel(
                lon=slice(*lon),
                lat=slice(*lat),
                time=slice(mrun.start_time(), mrun.end_time()),
            )

            mrun[obj_type] = final_object

        return wrapper

    return outer_wrapper
