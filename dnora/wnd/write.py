import numpy as np
import pandas as pd
from .. import msg
from copy import copy
from ..aux import check_if_folder, create_filename_obj, create_filename_time, add_file_extension, add_folder_to_filename

from .wnd_mod import ForcingWriter # Abstract class
from .wnd_mod import Forcing # Forcing object

from ..defaults import dflt_frc

class WW3(ForcingWriter):
    def __init__(self, folder: str='', filestring: str=dflt_frc['fs']['WW3'], datestring: str=dflt_frc['ds']['WW3']) -> None:
        self.folder = copy(folder)
        self.filestring = copy(filestring)
        self.datestring = copy(datestring)

        return

    def __call__(self, forcing: Forcing) -> None:
        msg.header(f'{type(self).__name__}: writing wind forcing from {forcing.name()}')

        existed = check_if_folder(folder=self.folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {self.folder}")

        output_file = forcing.filename(filestring=self.filestring, datestring=self.datestring, extension='nc')

        # Add folder
        output_path = add_folder_to_filename(output_file, folder=self.folder)

        msg.to_file(output_path)
        forcing.data.to_netcdf(output_path)

        return output_file, self.folder


class SWAN(ForcingWriter):
    def __init__(self, folder: str='', filestring: str=dflt_frc['fs']['SWAN'], datestring: str=dflt_frc['ds']['SWAN']) -> None:
        self.folder = copy(folder)
        self.filestring = copy(filestring)
        self.datestring = copy(datestring)

        return

    def __call__(self, forcing: Forcing) -> None:
        msg.header(f'{type(self).__name__}: writing wind forcing from {forcing.name()}')

        existed = check_if_folder(folder=self.folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {self.folder}")

        output_file = forcing.filename(filestring=self.filestring, datestring=self.datestring, extension='asc')

        # Add folder
        output_path = add_folder_to_filename(output_file, folder=self.folder)
        msg.to_file(output_path)

        days = forcing.days()
        with open(output_path, 'w') as file_out:
            for day in days:
                msg.plain(day.strftime('%Y-%m-%d'))
                times = forcing.times_in_day(day)
                for n in range(len(times)):
                    time_stamp = pd.to_datetime(
                        times[n]).strftime('%Y%m%d.%H%M%S')+'\n'
                    file_out.write(time_stamp)
                    np.savetxt(file_out, forcing.u()
                               [n, :, :]*1000, fmt='%i')
                    file_out.write(time_stamp)
                    np.savetxt(file_out, forcing.v()
                               [n, :, :]*1000, fmt='%i')


        return output_file, self.folder
