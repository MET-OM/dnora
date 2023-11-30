from __future__ import annotations

from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from typing import TYPE_CHECKING, Union

# Import objects
if TYPE_CHECKING:
    from ...modelrun.modelrun import ModelRun
    from ...file_module import FileNames


class IceWriter(ABC):
    """Writes the water level data to a certain file format.

    This object is provided to the .export_iceforcing() method.
    """

    @abstractmethod
    def __call__(
        self, model: ModelRun, file_object: FileNames, **kwargs
    ) -> Union[str, list[str]]:
        """Write the data from the Spectra object and returns the file and
        folder where data were written."""
        pass
