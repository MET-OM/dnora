from dnora.file_module import FileNames

from copy import copy
from dnora.type_manager.dnora_types import DnoraDataType
from dnora.type_manager.model_formats import ModelFormat
from dnora.read.abstract_readers import DataReader
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
from dnora import msg, utils
from geo_skeletons import PointSkeleton


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
            given_point_picker = kwargs.get("point_picker") or mrun._get_point_picker()

            if given_reader is None:
                raise ValueError("Provide a DataReader!")
            strategy = given_reader.caching_strategy()

            if dont_proceed_with_caching(read_cache, write_cache, strategy, kwargs):
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
                kwargs["expansion_factor"],
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
                msg.info("Expanding area to 5 degree daily tiles...")
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

            if write_cache:
                # Import data from original source only
                kwargs_cache["post_process"] = False
                mrun_cacher._import_data(**kwargs_cache)
                mrun_cacher[obj_type].name = name

                msg.header(
                    "Netcdf (DataWriter)",
                    f"Writing {obj_type.name} data from {name}",
                )
                write_data_to_cache(mrun_cacher, tiles, obj_type)
                msg.to_multifile(tiles.covering_files())

                mrun_cacher._post_process_object(obj_type, mrun_cacher._post_processing)
                read_cache = True

            if read_cache:
                kwargs_cache["point_mask"] = grid.boundary_mask()
                mrun_cacher = read_data_from_cache(
                    mrun_cacher, tiles, cache_reader, kwargs_cache
                )

                # Patch from original source
                if tiles.additional_files() and not kwargs.get("dont_patch"):
                    msg.blank()
                    msg.process(
                        f"Applying {strategy.name} caching strategy in patching."
                    )
                    mrun_cacher = patch_cached_data(
                        mrun_cacher, tiles, kwargs_cache, strategy
                    )
                mrun_cacher[obj_type].name = name

            ## Crop final object to the desired area since it might have been exanded to tiles
            final_object = mrun_cacher[obj_type]
            lon, lat = utils.grid.expand_area(
                mrun.grid().edges("lon", native=True),
                mrun.grid().edges("lat", native=True),
                expansion_factor=kwargs["expansion_factor"],
            )
            slice_dict = {"time": slice(mrun.start_time(), mrun.end_time())}

            msg.plain("")
            msg.print_line()
            msg.plain("CUTTING CACHED DATA")
            msg.print_line()
            orig_lon = mrun.grid().edges("lon", native=True)
            orig_lat = mrun.grid().edges("lat", native=True)
            msg.plain(
                f"Area: lon: ({orig_lon[0]:.2f}, {orig_lon[1]:.2f}), ({orig_lat[0]:.2f}, {orig_lat[1]:.2f})"
            )
            msg.process(f"Using expansion_factor = {kwargs['expansion_factor']:.2f}")
            msg.plain(
                f"Area: lon: ({lon[0]:.2f}, {lon[1]:.2f}), ({lat[0]:.2f}, {lat[1]:.2f})"
            )
            msg.plain(f"{mrun.start_time()} - {mrun.end_time()}")

            if final_object.is_gridded():

                slice_dict[grid.core.x_str] = slice(*lon)
                slice_dict[grid.core.y_str] = slice(*lat)

            else:
                # Get the wanted points from the exanded area using the original PointPicker
                msg.header(given_point_picker, "Choosing points to import...")
                inds = given_point_picker(
                    grid=mrun.grid(),
                    all_points=PointSkeleton.from_skeleton(mrun_cacher[obj_type]),
                    selected_points=PointSkeleton.from_skeleton(
                        mrun.grid(), mask=mrun.grid().sea_mask()
                    ),
                    expansion_factor=kwargs["expansion_factor"],
                )
                slice_dict["inds"] = inds

            if obj_type in [DnoraDataType.SPECTRA, DnoraDataType.SPECTRA1D]:
                convention = final_object.convention()
                final_object = final_object.sel(**slice_dict)
                final_object._mark_convention(convention)
            else:
                final_object = final_object.sel(**slice_dict)
            final_object.name = name

            mrun[obj_type] = final_object

            msg.header("<<< CACHE", "Exiting caching mode <<<")

        return wrapper

    return outer_wrapper
