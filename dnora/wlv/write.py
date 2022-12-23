from __future__ import annotations

from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from copy import copy
from typing import TYPE_CHECKING, Tuple

# Import objects
if TYPE_CHECKING:
    from .wlv_mod import WaterLevel # Boundary object

# Import default values and aux_funcsiliry functions
from .. import msg

class WaterLevelWriter(ABC):
    """Writes the water level data to a certain file format.

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
    def __call__(self, waterlevel: WaterLevel, filename: str) -> List[str]:
        """Writed the data from the Forcing object and returns the file and
        folder where data were written."""

        return output_file

# class WW3(ForcingWriter):
#     """Writes wind forcing data to WAVEWATH III netcdf format."""
#     def _extension(self):
#         return 'nc'
#
#     def __call__(self, forcing: Forcing, filename: str) -> List[str]:
#         forcing.data.to_netcdf(filename)
#         return filename


class SWAN(WaterLevelWriter):
    """Writes wind forcing data to SWAN ascii format."""

    def _extension(self):
        return 'asc'

    def __call__(self, waterlevel: WaterLevel, filename: str) -> List[str]:

        days = waterlevel.days()
        with open(filename, 'w') as file_out:
            ct = 0
            for day in days:
                msg.plain(day.strftime('%Y-%m-%d'))
                times = waterlevel.times_in_day(day)
                for n in range(len(times)):
                    time_stamp = pd.to_datetime(
                        times[n]).strftime('%Y%m%d.%H%M%S')+'\n'
                    file_out.write(time_stamp)
                    np.savetxt(file_out, waterlevel.waterlevel()
                               [ct, :, :]*1000, fmt='%i')
                    ct += 1

        return filename
