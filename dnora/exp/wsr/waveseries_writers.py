from __future__ import annotations # For TYPE_CHECKING
from typing import TYPE_CHECKING, Tuple
from abc import ABC, abstractmethod

# Import abstract classes and needed instances of them
if TYPE_CHECKING:
    from ...mdl.mdl_mod import ModelRun
    from ...file_module import FileNames

from ...aux_funcs import write_monthly_nc_files

class WaveSeriesWriter(ABC):
    """Writes WaveSeries data to a certain file format.

    This object is provided to the .export_waveseries() method.
    """

    def _im_silent(self) -> bool:
        """Return False if you want to be responsible for printing out the
        file names."""
        return True

    def _clean_filename(self) -> bool:
        """If this is set to False, then the ModelRun object does not clean
        the filename, and possible placeholders (e.g. #T0) can still be
        present.
        """

    @abstractmethod
    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> tuple[str, str]:
        """Write the data from the Spectra object and returns the file and
        folder where data were written."""

class Null(WaveSeriesWriter):
    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> list[str]:
        return ''

class DumpToNc(WaveSeriesWriter):
    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> tuple[str, str]:
        filename = file_object.get_filepath(extension='nc')
        model.waveseries().ds().to_netcdf(filename)
        return filename

class DnoraNc(WaveSeriesWriter):
    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> tuple[str, str]:
        output_files = write_monthly_nc_files(model.waveseries(), file_object)
        return output_files
