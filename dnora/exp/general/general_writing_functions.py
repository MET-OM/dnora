from __future__ import annotations
from abc import ABC, abstractmethod

# Import abstract classes and needed instances of them
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from ...mdl.mdl_mod import ModelRun
    from ...file_module import FileNames

from ...aux_funcs import write_monthly_nc_files


class GeneralWritingFunction(ABC):
    """General writing function that can be used for several objects"""

    @abstractmethod
    def __call__(
        self, model: ModelRun, file_object: FileNames, obj_type: str, **kwargs
    ) -> Union[str, list[str]]:
        return output_files


class Null(GeneralWritingFunction):
    def __call__(
        self, model: ModelRun, file_object: FileNames, obj_type: str, **kwargs
    ) -> str:
        return ""


class DnoraNc(GeneralWritingFunction):
    def __call__(
        self, model: ModelRun, file_object: FileNames, obj_type: str, **kwargs
    ) -> list[str]:
        output_files = write_monthly_nc_files(model[obj_type], file_object)
        return output_files


class DumpToNc(GeneralWritingFunction):
    def __call__(
        self, model: ModelRun, file_object: FileNames, obj_type: str, **kwargs
    ) -> str:
        output_file = file_object.get_filepath()
        model[obj_type].ds().to_netcdf(output_file)
        return output_file
