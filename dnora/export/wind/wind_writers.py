from __future__ import annotations

import numpy as np
import pandas as pd
from typing import TYPE_CHECKING

# Import objects
if TYPE_CHECKING:
    from ...modelrun.modelrun import ModelRun
    from ...file_module import FileNames

# Import default values and aux_funcsiliry functions
from ... import msg

from ..generic.generic_writers import DataWriter
from dnora.dnora_types import DnoraDataType


class SWAN(DataWriter):
    """Writes wind forcing data to SWAN ascii format."""

    _datavars = {DnoraDataType.WIND: ["u", "v"], DnoraDataType.WATERLEVEL: ["eta"]}

    def __call__(
        self, model: ModelRun, file_object: FileNames, obj_type: DnoraDataType, **kwargs
    ) -> list[str]:
        filename = file_object.get_filepath()
        data = model[obj_type]
        if data is None:
            msg.warning("Can't find any {obj_type.name} data!! Aborting...")
            return ""

        days = data.days(datetime=False)
        with open(filename, "w") as file_out:
            ct = 0
            for day in days:
                msg.plain(day)
                times = data.time(time=slice(day, day))
                for n in range(len(times)):
                    time_stamp = (
                        pd.to_datetime(times[n]).strftime("%Y%m%d.%H%M%S") + "\n"
                    )

                    if self._datavars.get(obj_type) is None:
                        msg.warning(
                            "Don't know how to write {obj_type.name} in SWAN format! Aborting..."
                        )
                        return ""
                    for var in self._datavars.get(obj_type):
                        file_out.write(time_stamp)
                        np.savetxt(file_out, data.get(var)[ct, :, :] * 1000, fmt="%i")

                    ct += 1

        return filename
