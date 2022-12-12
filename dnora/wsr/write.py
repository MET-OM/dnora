from __future__ import annotations # For TYPE_CHECKING
from typing import TYPE_CHECKING, Tuple
from abc import ABC, abstractmethod

# Import abstract classes and needed instances of them
if TYPE_CHECKING:
    from .wsr_mod import WaveSeries

from ..aux_funcs import write_monthly_nc_files

class WaveSeriesWriter(ABC):
    """Writes WaveSeries data to a certain file format.

    This object is provided to the .export_waveseries() method.
    """

    def _clean_filename(self):
        """If this is set to False, then the ModelRun object does not clean
        the filename, and possible placeholders (e.g. #T0) can still be
        present.
        """
        return True

    @abstractmethod
    def _extension(self) -> str:
        pass

    def _im_silent(self) -> bool:
        """Return False if you want to be responsible for printing out the
        file names."""
        return True

    @abstractmethod
    def __call__(self, waveseries: WaveSeries, file_object: str) -> tuple[str, str]:
        """Write the data from the Spectra object and returns the file and
        folder where data were written."""

class Null(WaveSeriesWriter):
    def _extension(self):
        return 'junk'

    def __call__(self, dict_of_objects: dict, file_object) -> List[str]:
        return ''

class DnoraNc(WaveSeriesWriter):
    def _extension(self) -> str:
        return 'nc'

    def __call__(self, dict_of_objects: dict, file_object) -> tuple[str, str]:
        output_files = write_monthly_nc_files(dict_of_objects['WaveSeries'], file_object)
        return output_files
