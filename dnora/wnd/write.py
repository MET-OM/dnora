from __future__ import annotations

from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from copy import copy
from typing import TYPE_CHECKING, Tuple

# Import objects
if TYPE_CHECKING:
    from .wnd_mod import Forcing # Boundary object

# Import default values and auxiliry functions
from ..defaults import dflt_frc, list_of_placeholders
from .. import msg
from ..aux import check_if_folder, add_folder_to_filename, clean_filename

class ForcingWriter(ABC):
    """Writes the forcing data to a certain file format.

    This object is provided to the .export_forcing() method.
    """
    @abstractmethod
    def _extension(self):
        pass

    def _im_silent(self) -> bool:
        """Return False if you want to be responsible for printing out the
        file names."""
        return True

    def _clean_filename(self):
        """If this is set to False, then the ModelRun object does not clean
        the filename, and possible placeholders (e.g. #T0) can still be
        present.
        """
        return True

    @abstractmethod
    def __call__(self, forcing: Forcing, filename: str, folder: str) -> List[str]:
        """Writed the data from the Forcing object and returns the file and
        folder where data were written."""

        return output_file

class WW3(ForcingWriter):
    """Writes wind forcing data to WAVEWATH III netcdf format."""
    def _extension(self):
        return 'nc'

    def __call__(self, forcing: Forcing, filename: str, folder: str) -> List[str]:
        # Add folder
        output_file = add_folder_to_filename(filename, folder=folder)

        forcing.data.to_netcdf(output_path)

        return output_file


class SWAN(ForcingWriter):
    """Writes wind forcing data to SWAN ascii format."""

    def _extension(self):
        return 'asc'

    def __call__(self, forcing: Forcing, filename: str, folder: str) -> List[str]:

        # Add folder
        output_file = add_folder_to_filename(filename, folder=folder)

        days = forcing.days()
        with open(output_file, 'w') as file_out:
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

        return output_file
