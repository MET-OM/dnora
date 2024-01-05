from __future__ import annotations
from abc import ABC, abstractmethod

# Import abstract classes and needed instances of them
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from ...modelrun.modelrun import ModelRun
    from ...file_module import FileNames


from calendar import monthrange
from ...dnora_object_type import DnoraDataType
from ... import file_module
import os


class GenericWriter(ABC):
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


class Null(GenericWriter):
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        obj_type: DnoraDataType,
        **kwargs,
    ) -> str:
        return ""


class DnoraNc(GenericWriter):
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        obj_type: DnoraDataType,
        **kwargs,
    ) -> list[str]:
        output_files = write_monthly_nc_files(model[obj_type], file_object)
        return output_files


class DumpToNc(GenericWriter):
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        obj_type: DnoraDataType,
        **kwargs,
    ) -> str:
        output_file = file_object.get_filepath()
        model[obj_type].ds().to_netcdf(output_file)
        return output_file


def write_monthly_nc_files(
    dnora_obj: DnoraDataType, file_object: FileNames
) -> list[str]:
    "Writes the data of a DNORA object into montly netcdf-files wh the ames specified by the FileNames instance."
    output_files = []
    for month in dnora_obj.months():
        t0 = f"{month.strftime('%Y-%m-01')}"
        d1 = monthrange(int(month.strftime("%Y")), int(month.strftime("%m")))[1]
        t1 = f"{month.strftime(f'%Y-%m-{d1}')}"

        outfile = file_object.get_filepath(
            start_time=month, edge_object=DnoraDataType.Grid
        )

        outfile = file_module.clean_filename(outfile)
        if os.path.exists(outfile):
            os.remove(outfile)
        dnora_obj.ds().sel(time=slice(t0, t1)).to_netcdf(outfile)

        output_files.append(outfile)
    return output_files
