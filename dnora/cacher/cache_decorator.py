from dnora.file_module import FileNames

from copy import copy
from dnora.aux_funcs import expand_area
from dnora.dnora_type_manager.dnora_types import DnoraDataType
from dnora.model_formats import ModelFormat
from dnora.readers.generic_readers import DataReader
from inspect import getcallargs

from dnora.pick import Area
from .tiles import TileObject
from .caching_functions import (
    dont_proceed_with_caching,
    expand_area_to_tiles,
    read_data_from_cache,
    patch_cached_data,
    write_data_to_cache,
)
from dnora import msg


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
            read_cache: bool = False,
            write_cache: bool = False,
            **kwargs,
        ):

            # Get all arguments as keyword arguments
            kwargs = get_kwargs(import_method, args, kwargs)
            mrun = kwargs.get("self")

            given_reader = kwargs.get("reader") or mrun._get_reader(obj_type)
            if given_reader is None:
                raise ValueError("Provide a DataReader!")

            if dont_proceed_with_caching(read_cache, write_cache, given_reader, kwargs):
                import_method(**kwargs)
                return

            msg.header(">>> CACHE", "Entering caching mode >>>")
            ### Helper objects
            ### -------------------------------------------------------------
            # This name will be used in the folders
            # Default use name of the reader, but user can override
            name = kwargs.get("name")
            if name is None:
                name = given_reader.name()

            # Define file_obect to get proper file name
            file_object = FileNames(
                format=ModelFormat.CACHE,
                obj_type=obj_type,
                obj_name=name,
                model=mrun,
                edge_object=DnoraDataType.GRID,
            )

            # Covering tiles are daily 5degx5deg tiles
            # Covers the (possibly expanded) area we want to read
            tiles = TileObject(file_object)
            tiles.create_covering_files(
                mrun.grid(),
                mrun.start_time(),
                mrun.end_time(),
                kwargs.get("expansion_factor", 1.0),
            )
            ### -------------------------------------------------------------

            ### Determine new kwargs to be given to caching function calls
            kwargs_cache = copy(kwargs)
            del kwargs_cache["self"]
            kwargs_cache["obj_type"] = obj_type
            start_time, end_time = mrun.start_time(), mrun.end_time()
            grid = mrun.grid()

            ## We should always write full tiles to cache, so we need to read an area covering the tiles
            if write_cache:
                kwargs_cache["point_picker"] = Area()
                kwargs_cache["expansion_factor"] = 1.0
                (
                    grid,
                    start_time,
                    end_time,
                ) = expand_area_to_tiles(tiles, grid.dlon(), grid.dlat())

            ## Reading of the data starts here
            mrun_cacher = mrun.empty_copy(
                grid=grid,
                start_time=start_time,
                end_time=end_time,
            )
            if read_cache:
                mrun_cacher = read_data_from_cache(
                    mrun_cacher, tiles, cache_reader, kwargs_cache
                )
                # Patch from original source
                mrun_cacher = patch_cached_data(mrun_cacher, tiles, kwargs_cache)
            else:
                # Import data from original source only
                mrun_cacher._import_data(**kwargs_cache)

            mrun_cacher[obj_type].name = name

            ## Write data if necessary
            if write_cache:
                msg.header(
                    "Netcdf (DataWriter)",
                    f"Writing {obj_type.name} data from {name}",
                )
                write_data_to_cache(mrun_cacher, tiles, obj_type)
                msg.to_multifile(tiles.covering_files())

            ## Crop final object to the desired area since it might have been exanded to tiles
            final_object = mrun_cacher[obj_type]
            lon, lat = expand_area(
                grid.edges("lon", native=True),
                grid.edges("lat", native=True),
                expansion_factor=kwargs.get("expansion_factor", 1.0),
            )

            slice_dict = {
                "time": slice(mrun.start_time(), mrun.end_time()),
                grid.x_str: slice(*lon),
                grid.y_str: slice(*lat),
            }

            final_object = final_object.sel(**slice_dict)
            final_object.name = name
            mrun[obj_type] = final_object

            msg.header("<<< CACHE", "Exiting caching mode <<<")

        return wrapper

    return outer_wrapper
