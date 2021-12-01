from __future__ import annotations
from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from .. import msg
from copy import copy
from ..aux import check_if_folder, add_folder_to_filename, clean_filename

#from .wnd_mod import ForcingWriter # Abstract class
from typing import TYPE_CHECKING, Tuple
if TYPE_CHECKING:
    from .wnd_mod import Forcing # Boundary object
from ..defaults import dflt_frc

class ForcingWriter(ABC):
    """Writes the forcing data to a certain file format.

    This object is provided to the .export_forcing() method.
    """
    def _preferred_format(self):
        return 'General'

    @abstractmethod
    def __call__(self, forcing: Forcing, filename: str, folder: str) -> Tuple[str, str]:
        """Writed the data from the Forcing object and returns the file and
        folder where data were written."""

        return output_file, output_folder

class WW3(ForcingWriter):
    """Writes wind forcing data to WAVEWATH III netcdf format."""
    def _preferred_format(self):
        return 'WW3'

    def __call__(self, forcing: Forcing, filename: str, folder: str) -> Tuple[str, str]:


        # Add folder
        output_path = add_folder_to_filename(filename, folder=folder)

        msg.to_file(output_path)
        forcing.data.to_netcdf(output_path)

        return filename, folder


class SWAN(ForcingWriter):
    """Writes wind forcing data to SWAN ascii format."""
    def __init__(self, out_format = 'SWAN'):
        self.out_format = out_format
        return

    def _preferred_format(self):
        return self.out_format

    def __call__(self, forcing: Forcing, filename: str, folder: str) -> Tuple[str, str]:

        # Add folder
        output_path = add_folder_to_filename(filename, folder=folder)
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

        return filename, folder
