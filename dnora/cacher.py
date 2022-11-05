from . import file_module
from . import msg
from . import aux_funcs
import glob, os, re
from calendar import monthrange
import pandas as pd
import numpy as np
from dataclasses import dataclass
from .file_module import FileNames
@dataclass
class Cacher:
    dnora_obj: str
    cache_name: str

    def __post_init__(self):
        self.file_object = FileNames(format='Cache',
                                dnora_obj=self.dnora_obj,
                                edges_from_grid=True,
                                extension='nc',
                                _filename=self.cache_name
                                )
        self.file_object.create_folder()

    def empty(self):
        return not glob.glob(f'{self.filepath(extension=False)}*')

    def filename(self, start_time: str=None, end_time: str=None):
        return self.file_object.filename(start_time=start_time, end_time=end_time)

    def folder(self):
        return self.file_object.folder()

    def filepath(self, start_time: str=None, end_time: str=None, extension: bool=True):
        if extension:
            return self.file_object.filepath(start_time=start_time, end_time=end_time)
        return self.file_object.filepath(start_time=start_time, end_time=end_time)[0:-3]

    def determine_patch_periods(self, start_time, end_time):
        """Determines if there is some periods that we need to patch from thredds
        adter reading cached data"""

        # This is not optimal, but seems to work
        dt = self.dnora_obj.time()[1]-self.dnora_obj.time()[0]
        wanted_times = pd.date_range(start=start_time, end=end_time, freq=dt)
        wanted_times.isin(self.dnora_obj.time())

        if np.all(wanted_times.isin(self.dnora_obj.time())):
            return [], []
        wt=wanted_times.isin(self.dnora_obj.time())

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

    def write_cache(self) -> None:
        for month in self.dnora_obj.months():
            t0 = f"{month.strftime('%Y-%m-01')}"
            d1 = monthrange(int(month.strftime('%Y')), int(month.strftime('%m')))[1]
            t1 = f"{month.strftime(f'%Y-%m-{d1}')}"
            outfile = self.filepath(start_time=t0)
            if os.path.exists(outfile):
                os.remove(outfile)
            self.dnora_obj.ds().sel(time=slice(t0, t1)).to_netcdf(outfile)
            msg.to_file(outfile)
