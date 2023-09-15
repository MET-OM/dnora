from ..bnd import Boundary
from ..wnd import Forcing
from ..spc import Spectra
from ..wsr import WaveSeries
from ..file_module import FileNames
import glob, os
from copy import copy
from calendar import monthrange
from .. import msg
import pandas as pd
import numpy as np
import re


def cached_reader(dnora_obj, reader_function):

    def outer_wrapper(import_method):
        def wrapper(*args, cache_name: str=None, read_cache: bool=False, write_cache: bool=False, patch: bool=False, **kwargs):
            def get_reader(args, kwargs):
                reader = None
                if kwargs.get(f'{dnora_obj}_reader') is not None:
                    reader = kwargs.get(f'{dnora_obj}_reader')
                else:
                    new_args = []
                    for arg in args:
                        if reader_function.__bases__[0] == arg.__class__.__bases__[0]: # HAve we found the e.g. BoundaryReader
                            reader = arg
                reader = reader or eval(f'args[0]._get_{dnora_obj.lower()}_reader()')
                return reader

            def determine_patch_periods():
                """Determines if there is some periods that we need to patch from thredds
                adter reading cached data"""

                # This is not optimal, but seems to work
                times =  mrun.dict_of_objects().get(dnora_obj).time()
                dt =times[1]-times[0]
                wanted_times = pd.date_range(start=mrun.start_time(), end=mrun.end_time(), freq=dt)
                wanted_times.isin(times)

                if np.all(wanted_times.isin(times)):
                    return [], []
                wt=wanted_times.isin(times)

                was_found = ''.join([str((w*1)) for w in wt]) # string of '0001110111110000'
                #was_found = '000011111111111100011111111111111111111111100011111111111111011111111111110000' # Testing
                inds = list(range(len(was_found)))
                was_found = re.sub('01', '0.1', was_found)
                was_found = re.sub('10', '1.0', was_found)

                list_of_blocks = was_found.split('.')

                patch_start = []
                patch_end = []
                for block in list_of_blocks:
                    if block[0] == '0': # These need to be patched
                        ind_subset = inds[0:len(block)]
                        patch_start.append(wanted_times[ind_subset[0]])
                        patch_end.append(wanted_times[ind_subset[-1]])
                    inds[0:len(block)] = []

                return patch_start, patch_end

            if cache_name is not None and '#T0' not in cache_name:
                cache_name = cache_name + '_#T0'
            mrun = args[0]

            name = kwargs.get('name')
            if name is None:
                given_reader = get_reader(args, kwargs)
                if given_reader is not None:
                    name = given_reader.name()

            # FileObject takes names from the dict of objects, so create one here

            exec(f"mrun._{dnora_obj.lower()} = {dnora_obj}(name=name)")
            file_object = FileNames(format='Cache',
                                    obj_type=dnora_obj,
                                    extension='nc',
                                    dict_of_objects=mrun.dict_of_objects(),
                                    edge_object='Grid',
                                    filename=cache_name
                                    )
            file_object.create_folder()

            files = glob.glob(f'{file_object.get_filepath()[0:-3]}*')

            reader = reader_function(files=files)

            if files and read_cache:
                new_kwargs = copy(kwargs)
                new_kwargs['name'] = name
                if kwargs.get(f'{dnora_obj.lower()}_reader') is not None:
                    new_kwargs[f'{dnora_obj.lower()}_reader'] = reader
                    new_args = args
                else:
                    new_args = []
                    for arg in args:
                        if not reader.__class__.__bases__[0] == arg.__class__.__bases__[0]:
                            new_args.append(arg)
                    new_args = tuple(new_args)


                new_kwargs[f'{dnora_obj.lower()}_reader'] = reader
                import_method(*new_args, **new_kwargs)
            else:
                import_method(*args, **kwargs)


            if files and read_cache:
                patch_start, patch_end = determine_patch_periods()
                if patch_start and patch:
                    msg.info('Not all data found in cache. Patching from original source...')
                    temp_args = tuple(list(args)[1:])
                    for t0, t1 in zip(patch_start, patch_end):
                        mrun_temp = copy(mrun)
                        mrun_temp.start_time = t0
                        mrun_temp.end_time = t1
                        exec(f'mrun_temp.import_{dnora_obj.lower()}(*temp_args, **kwargs)')
                        mrun.dict_of_objects().get(dnora_obj)._absorb_object(mrun_temp.dict_of_objects().get(dnora_obj), 'time')

            if write_cache:
                exec(f'mrun.cache_{dnora_obj.lower()}()')
        return wrapper
    return outer_wrapper
