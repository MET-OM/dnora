from __future__ import annotations

from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from typing import TYPE_CHECKING

# Import objects
if TYPE_CHECKING:
    from ...mdl.mdl_mod import ModelRun
    from ...file_module import FileNames

# Import default values and aux_funcsiliry functions
from ... import msg

from ...aux_funcs import write_monthly_nc_files

class ForcingWriter(ABC):
    """Writes the forcing data to a certain file format.

    This object is provided to the .export_forcing() method.
    """
    def _im_silent(self) -> bool:
        """Return False if you want to be responsible for printing out the
        file names."""
        return True

    def _clean_filename(self) -> bool:
        """If this is set to False, then the ModelRun object does not clean
        the filename, and possible placeholders (e.g. #T0) can still be
        present.
        """

    @abstractmethod
    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> list[str]:
        """Writed the data from the Forcing object and returns the file and
        folder where data were written."""

        return output_file

class Null(ForcingWriter):
    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs):
        return ''

class DnoraNc(ForcingWriter):
    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> tuple[str, str]:
        output_files = write_monthly_nc_files(model.forcing(), file_object)
        return output_files

class WW3(ForcingWriter):
    """Writes wind forcing data to WAVEWATH III netcdf format."""
    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> list[str]:
        filename = file_object.get_filepath()
        model.forcing().ds().to_netcdf(filename)

        return filename


class SWAN(ForcingWriter):
    """Writes wind forcing data to SWAN ascii format."""
    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> list[str]:
        filename = file_object.get_filepath()
        forcing = model.forcing()

        days = forcing.days(datetime=False)
        with open(filename, 'w') as file_out:
            ct = 0
            for day in days:
                msg.plain(day)
                #times = forcing.times_in_day(day)
                times = forcing.time(time=slice(day, day))
                for n in range(len(times)):
                    time_stamp = pd.to_datetime(
                        times[n]).strftime('%Y%m%d.%H%M%S')+'\n'
                    file_out.write(time_stamp)
                    np.savetxt(file_out, forcing.u()
                               [ct, :, :]*1000, fmt='%i')
                    file_out.write(time_stamp)
                    np.savetxt(file_out, forcing.v()
                               [ct, :, :]*1000, fmt='%i')
                    ct += 1

        return filename
