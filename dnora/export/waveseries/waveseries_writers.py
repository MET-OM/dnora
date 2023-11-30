from __future__ import annotations  # For TYPE_CHECKING
from typing import TYPE_CHECKING
from abc import ABC, abstractmethod

# Import abstract classes and needed instances of them
if TYPE_CHECKING:
    from ...modelrun.modelrun import ModelRun
    from ...file_module import FileNames


class WaveSeriesWriter(ABC):
    """Writes WaveSeries data to a certain file format.

    This object is provided to the .export_waveseries() method.
    """

    @abstractmethod
    def __call__(
        self, model: ModelRun, file_object: FileNames, **kwargs
    ) -> list[str, str]:
        """Write the data from the Spectra object and returns the file and
        folder where data were written."""
        pass
