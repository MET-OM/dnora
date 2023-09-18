from __future__ import annotations

from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from typing import TYPE_CHECKING

# Import objects
if TYPE_CHECKING:
    from .wlv_mod import WaterLevel # Boundary object

# Import default values and aux_funcsiliry functions
from .. import msg
from ..aux_funcs import write_monthly_nc_files

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
    def __call__(self, waterlevel: WaterLevel, filename: str) -> list[str]:
        """Writed the data from the Forcing object and returns the file and
        folder where data were written."""

        return output_file

class Null(WaterLevelWriter):
    def _extension(self):
        return 'junk'

    def __call__(self, dict_of_objects: dict, file_object):
        return ''

class DnoraNc(WaterLevelWriter):
    def _extension(self) -> str:
        return 'nc'

    def __call__(self, dict_of_objects: dict, file_object) -> tuple[str, str]:
        output_files = write_monthly_nc_files(dict_of_objects['WaterLevel'], file_object)
        return output_files


class SWAN(WaterLevelWriter):
    """Writes wind forcing data to SWAN ascii format."""

    def _extension(self):
        return 'asc'

    def __call__(self, waterlevel: WaterLevel, filename: str) -> list[str]:

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