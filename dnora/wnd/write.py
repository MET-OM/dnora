import numpy as np
import pandas as pd
from .. import msg
from copy import copy
from ..aux import check_if_folder, create_filename_obj, create_filename_time, add_file_extension, add_folder_to_filename

from .wnd_mod import ForcingWriter # Abstract class
from .wnd_mod import Forcing # Forcing object

from ..defaults import frc

class WW3(ForcingWriter):
    def __init__(self, folder: str='', filestring: str=frc['fs']['WW3'], datestring: str=frc['ds']['WW3']) -> None:
        self.folder = copy(folder)
        self.filestring = copy(filestring)
        self.datestring = copy(datestring)

        return

    def __call__(self, forcing: Forcing) -> None:
        msg.header(f'{type(self).__name__}: writing wind forcing from {forcing.name()}')

        existed = check_if_folder(folder=self.folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {self.folder}")

        ### Creating the filename based on the provided filestring and datestring
        # Substitute placeholders for objects ($Grid etc.)
        output_file = create_filename_obj(filestring=self.filestring, objects=[forcing, forcing.grid])
        # Substitute placeholders for times ($T0 etc.)
        output_file = create_filename_time(filestring=output_file, times=[forcing.time()[0], forcing.time()[-1]], datestring=self.datestring)
        # Add extension .nc if doesn't exist
        output_file = add_file_extension(output_file, extension='nc')

        # Add folder
        output_path = add_folder_to_filename(output_file, folder=self.folder)

        msg.to_file(output_path)
        forcing.data.to_netcdf(output_path)

        # This is set as info in case an input file needs to be generated
        forcing._written_as = output_file
        forcing._written_to = self.folder

        return


class SWAN(ForcingWriter):
    def __init__(self, folder: str='', filestring: str=frc['fs']['SWAN'], datestring: str=frc['ds']['SWAN']) -> None:
        self.folder = copy(folder)
        self.filestring = copy(filestring)
        self.datestring = copy(datestring)

        return

    def __call__(self, forcing: Forcing) -> None:
        msg.header(f'{type(self).__name__}: writing wind forcing from {forcing.name()}')

        existed = check_if_folder(folder=self.folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {self.folder}")

        ### Creating the filename based on the provided filestring and datestring
        # Substitute placeholders for objects ($Grid etc.)
        output_file = create_filename_obj(filestring=self.filestring, objects=[forcing, forcing.grid])
        # Substitute placeholders for times ($T0 etc.)
        output_file = create_filename_time(filestring=output_file, times=[pd.Timestamp(forcing.time()[0]), pd.Timestamp(forcing.time()[-1])], datestring=self.datestring)
        # Add extension .asc if doesn't exist
        output_file = add_file_extension(output_file, extension='asc')

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

        # This is set as info in case an input file needs to be generated
        forcing._written_as = output_file
        forcing._written_to = self.folder

        return
