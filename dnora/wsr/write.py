from __future__ import annotations # For TYPE_CHECKING
from typing import TYPE_CHECKING, Tuple
from abc import ABC, abstractmethod

# Import abstract classes and needed instances of them
if TYPE_CHECKING:
    from .wsr_mod import WaveSeries

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
    def __call__(self, waveseries: WaveSeries, filename: str) -> tuple[str, str]:
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

    def __call__(self, dict_of_objects: dict, file_object) -> Tuple[str, str]:
        output_files = []
        waveseries = dict_of_objects['WaveSeries']
        for month in waveseries.months():
            t0 = f"{month.strftime('%Y-%m-01')}"
            d1 = monthrange(int(month.strftime('%Y')), int(month.strftime('%m')))[1]
            t1 = f"{month.strftime(f'%Y-%m-{d1}')}"

            outfile = file_object.get_filepath(start_time=month, edge_object='Grid')

            outfile = file_object.clean(outfile)
            if os.path.exists(outfile):
                os.remove(outfile)
            waveseries.ds().sel(time=slice(t0, t1)).to_netcdf(outfile)

            output_files.append(outfile)
        return output_files

# class DumpToNc(WaveSeriesWriter):
#     def _extension(self) -> str:
#         return 'nc'
#
#     def __call__(self, dict_of_objects: dict, file_object) -> List[str]:
#         filename = file_object.get_filepath()
#         wavesereies = dict_of_objects['WaveSeries']
#
#         waveseries.ds().to_netcdf(filename)
#
#         return filename
