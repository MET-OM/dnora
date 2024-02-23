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

from dnora.dnora_types import DnoraDataType
from dnora.model_formats import ModelFormat
from dnora.readers.generic_readers import DataReader


def cached_reader(obj_type: DnoraDataType, reader_function: DataReader):
    def outer_wrapper(import_method):  # import_method is e.g. modelrun.import_spectra()
        def wrapper(
            *args,
            cache_name: str = None,
            read_cache: bool = False,
            write_cache: bool = False,
            patch: bool = False,
            **kwargs,
        ):
            def get_reader(args, kwargs):
                reader = None
                if kwargs.get("reader") is not None:
                    reader = kwargs.get("reader")
                else:
                    # new_args = []
                    for arg in args:
                        if isinstance(arg, DataReader):
                            reader = arg
                reader = reader or args[0]._get_reader(obj_type)
                return reader

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

            if cache_name is not None and "#T0" not in cache_name:
                cache_name = cache_name + "_#T0"
            mrun = args[0]

            name = kwargs.get("name")
            if name is None:
                given_reader = get_reader(args, kwargs)
                if given_reader is not None:
                    name = given_reader.name()
            if name is None:
                raise ValueError("Provide a DataReader!")

            file_object = FileNames(
                format=ModelFormat.CACHE,
                obj_type=obj_type,
                obj_name=name,
                model=mrun,
                edge_object=DnoraDataType.GRID,
                filename=cache_name,
            )

            create_folder = (
                not (kwargs.get("dry_run", False) or mrun.dry_run()) and write_cache
            )

            if create_folder:
                file_object.create_folder()
            files = glob.glob(f"{file_object.get_filepath()[0:-3]}*")
            reader = reader_function(files=files)

            if files and read_cache:
                new_kwargs = copy(kwargs)
                new_kwargs["name"] = name
                if kwargs.get("reader") is not None:
                    new_kwargs["reader"] = reader
                    new_args = args
                else:
                    new_args = []
                    for arg in args:
                        if (
                            not reader.__class__.__bases__[0]
                            == arg.__class__.__bases__[0]
                        ):
                            new_args.append(arg)
                    new_args = tuple(new_args)

                new_kwargs["reader"] = reader
                import_method(*new_args, **new_kwargs)
            else:
                import_method(*args, **kwargs)

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
