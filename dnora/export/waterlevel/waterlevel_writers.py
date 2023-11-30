from __future__ import annotations

from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from typing import TYPE_CHECKING

# Import objects
if TYPE_CHECKING:
    from ...modelrun.modelrun import ModelRun
    from ...file_module import FileNames

# Import default values and aux_funcsiliry functions
from ... import msg


class WaterLevelWriter(ABC):
    """Writes the water level data to a certain file format.

    This object is provided to the .export_waterlevel() method.
    """

    @abstractmethod
    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> list[str]:
        """Writed the data from the WaterLevel object and returns the file and
        folder where data were written."""
        pass


class SWAN(WaterLevelWriter):
    """Writes wind forcing data to SWAN ascii format."""

    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> list[str]:
        filename = file_object.get_filepath()
        waterlevel = model.waterlevel()

        days = waterlevel.days(datetime=False)
        with open(filename, "w") as file_out:
            ct = 0
            for day in days:
                msg.plain(day)
                times = waterlevel.time(time=slice(day, day))
                for n in range(len(times)):
                    time_stamp = (
                        pd.to_datetime(times[n]).strftime("%Y%m%d.%H%M%S") + "\n"
                    )
                    file_out.write(time_stamp)
                    np.savetxt(
                        file_out, waterlevel.waterlevel()[ct, :, :] * 1000, fmt="%i"
                    )
                    ct += 1

        return filename
