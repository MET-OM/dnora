from __future__ import annotations

from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from typing import TYPE_CHECKING, Union

# Import objects
if TYPE_CHECKING:
    from ...mdl.mdl_mod import ModelRun
    from ...file_module import FileNames

# Import default values and aux_funcsiliry functions
from ... import msg

class IceWriter(ABC):
    """Writes the water level data to a certain file format.

    This object is provided to the .export_ice() method.
    """

    @abstractmethod
    def __call__(
        self, model: ModelRun, file_object: FileNames, **kwargs
    ) -> Union[str, list[str]]:
        """Write the data from the Spectra object and returns the file and
        folder where data were written."""
        pass

class WW3(IceWriter):
    """Writes  ice data to WAVEWATH III netcdf format."""

    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> list[str]:
        filename = file_object.get_filepath()
        model.ice().ds().to_netcdf(filename)

        return filename


class SWAN(IceWriter):
    """Writes ice data to SWAN ascii format."""

    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> list[str]:
        filename = file_object.get_filepath()
        ice = model.ice()
        breakpoint()

        days = ice.days(datetime=False)
        with open(filename, "w") as file_out:
            ct = 0
            for day in days:
                msg.plain(day)
                # times = ice.times_in_day(day)
                times = ice.time(time=slice(day, day))
                for n in range(len(times)):
                    time_stamp = (
                        pd.to_datetime(times[n]).strftime("%Y%m%d.%H%M%S") + "\n"
                    )
                    file_out.write(time_stamp)
                    np.savetxt(file_out, ice.sic()[ct, :, :] * 1000, fmt="%i")
                    ct += 1

        return filename
