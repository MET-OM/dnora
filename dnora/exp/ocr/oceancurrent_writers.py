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
from ...aux_funcs import write_monthly_nc_files


class OceanCurrentWriter(ABC):
    """Writes the ocean current data to a certain file format.

    This object is provided to the .export_oceancurrent() method.
    """

    @abstractmethod
    def __call__(
        self, model: ModelRun, file_object: FileNames, **kwargs
    ) -> Union[str, list[str]]:
        pass