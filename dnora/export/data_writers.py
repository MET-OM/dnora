from __future__ import annotations
from abc import ABC, abstractmethod

# Import abstract classes and needed instances of them
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from dnora.modelrun.modelrun import ModelRun
    from dnora.file_module import FileNames


from calendar import monthrange
from dnora.type_manager.dnora_types import DnoraDataType
from dnora.type_manager.spectral_conventions import SpectralConvention

from dnora import file_module
import os
import numpy as np
import pandas as pd
from dnora import msg


class DataWriter(ABC):
    """General writing function that can be used for several objects"""

    @abstractmethod
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        obj_type: DnoraDataType,
        **kwargs,
    ) -> Union[str, list[str]]:
        pass


class Null(DataWriter):
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        obj_type: DnoraDataType,
        **kwargs,
    ) -> str:
        return ""


class Netcdf(DataWriter):
    def __init__(self, monthly_files: bool = False, daily_files: bool = False):
        self._monthly_files = monthly_files
        self._daily_files = daily_files

    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        obj_type: DnoraDataType,
        monthly_files: bool = False,
        daily_files: bool = False,
        **kwargs,
    ) -> str:
        monthly_files = monthly_files or self._monthly_files
        daily_files = daily_files or self._daily_files
        if monthly_files and daily_files:
            raise ValueError(f"Choose montly_files OR daily_files!")
        if monthly_files:
            output_files = write_monthly_nc_files(model[obj_type], file_object)
        elif daily_files:
            output_files = write_daily_nc_files(model[obj_type], file_object)
        else:
            output_files = file_object.get_filepath()
            model[obj_type].ds().to_netcdf(output_files)

        return output_files


def write_monthly_nc_files(
    dnora_obj: DnoraDataType, file_object: FileNames
) -> list[str]:
    "Writes the data of a DNORA object into montly netcdf-files wh the ames specified by the FileNames instance."
    output_files = []
    for month in dnora_obj.months():
        t0 = f"{month.strftime('%Y-%m-01')}"
        d1 = monthrange(int(month.strftime("%Y")), int(month.strftime("%m")))[1]
        t1 = f"{month.strftime(f'%Y-%m-{d1}')}"

        outfile = file_object.get_filepath(start_time=month)

        outfile = file_module.clean_filename(outfile)
        if os.path.exists(outfile):
            os.remove(outfile)
        dnora_obj.ds().sel(time=slice(t0, t1)).to_netcdf(outfile)

        output_files.append(outfile)
    return output_files


def write_daily_nc_files(dnora_obj: DnoraDataType, file_object: FileNames) -> list[str]:
    """Writes the data of a DNORA object into daily netcdf-files wh the ames specified by the FileNames instance.
    if the dnora_object has no xr.Dataset, then empty files will be created (needed in caching)
    """
    output_files = []
    for day in dnora_obj.days():
        t0 = f"{day.strftime('%Y-%m-%d 00:00')}"
        t1 = f"{day.strftime(f'%Y-%m-%d 23:59')}"

        outfile = file_object.get_filepath(start_time=day)

        outfile = file_module.clean_filename(outfile)
        if os.path.exists(outfile):
            os.remove(outfile)
        file_object.create_folder(start_time=day)
        # if dnora_obj.meta.get().get("dummy") == "True":
        #     with open(outfile, "w") as __:
        #         pass
        # else:

        dnora_obj.ds().sel(time=slice(t0, t1)).to_netcdf(outfile)

        output_files.append(outfile)
    return output_files


class SWAN(DataWriter):
    """Writes wind forcing data to SWAN ascii format."""

    # The names of the variables could in theory be read from the object class itself (obj_type.value._coord_manager.added_vars().keys())
    # We prefer to be explicit here, since then we are 100% sure that the order is correct
    def convention(self):
        return SpectralConvention.MET

    _datavars = {
        DnoraDataType.WIND: ["u", "v"],
        DnoraDataType.WATERLEVEL: ["eta"],
        DnoraDataType.ICE: ["sic"],
        DnoraDataType.CURRENT: ["u", "v"],
    }

    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        obj_type: DnoraDataType,
        data_vars: list[str] = None,
        **kwargs,
    ) -> list[str]:
        filename = file_object.get_filepath()
        data = model[obj_type]
        if data is None:
            msg.warning(f"Can't find any {obj_type.name} data!! Aborting...")
            return ""

        data_vars = data_vars or self._datavars.get(obj_type)
        if not isinstance(data_vars, list):
            data_vars = [data_vars]

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
                            f"Don't know how to write {obj_type.name} in SWAN format! Aborting..."
                        )
                        return ""
                    for var in data_vars:
                        file_out.write(time_stamp)
                        np.savetxt(
                            file_out, data.get(var)[ct, :, :] * 1000, fmt="%4.0f"
                        )

                    ct += 1

        return filename
