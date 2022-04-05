from __future__ import annotations # For TYPE_CHECKING

import numpy as np
from copy import copy
from abc import ABC, abstractmethod
import netCDF4
import re

# Import abstract classes and needed instances of them
from typing import TYPE_CHECKING, Tuple
if TYPE_CHECKING:
    from .spc_mod import Spectra

from .. import msg
from ..aux import add_folder_to_filename, clean_filename
from ..defaults import list_of_placeholders

class SpectralWriter(ABC):
    """Writes omnidirectional spectra spectra to a certain file format.

    This object is provided to the .export_spectra() method.
    """

    def _preferred_format(self) -> str:
        """For the file format using defauts.py, e.g. SWAN or WW3"""
        return 'General'

    def _preferred_extension(self) -> str:
        return 'nc'

    def _im_silent(self) -> bool:
        """Return False if you want to be responsible for printing out the
        file names."""
        return True

    @abstractmethod
    def _convention_in(self) -> str:
        """Defines in which format the incoming spectra should be.

        The conventions to choose from are predetermined:

        'Ocean':    Oceanic convention (direction to)

        'Met':      Meteorological convention (direction from)
        """

    @abstractmethod
    def __call__(self, spectra: Spectra, filename: str, folder: str) -> Tuple[str, str]:
        """Write the data from the Spectra object and returns the file and
        folder where data were written."""

class DumpToNc(SpectralWriter):
    def __init__(self, convention: str='Met', out_format: str='General') -> None:
        self._convention = convention
        self.out_format = out_format
        return

    def _convention_in(self) -> str:
        """Convention of spectra"""
        return self._convention

    def _preferred_format(self) -> str:
        """Preferred format of file name"""
        return self.out_format

    def __call__(self, spectra: Spectra, filename: str, folder: str) -> Tuple[str, str]:

        output_file = clean_filename(filename, list_of_placeholders)

        # Add folder
        output_path = add_folder_to_filename(output_file, folder=folder)

        # Dumping to a netcdf-file
        msg.to_file(output_path)
        spectra.data.to_netcdf(output_path)

        return output_file, folder

class REEF3D(SpectralWriter):
    def __init__(self, convention: str='Met', out_format: str='REEF3D') -> None:
        self._convention = convention
        self.out_format = out_format
        return

    def _convention_in(self) -> str:
        """Convention of spectra"""
        return self._convention

    def _preferred_format(self) -> str:
        """Preferred format of file name"""
        return self.out_format

    def _preferred_extension(self) -> str:
        return 'dat'

    def __call__(self, spectra: Spectra, filename: str, folder: str) -> Tuple[str, str]:

        output_file = clean_filename(filename, list_of_placeholders)

        # Add folder
        output_path = add_folder_to_filename(output_file, folder=folder)

        # Take first spectra and first time step for now
        x = 0
        t = 0
        with open(output_path, 'w') as f:
            spec = spectra.data.spec.values[x,t,:]
            freq = spectra.data.freq.values
            freq = freq*2*np.pi
            spec = spec/2/np.pi
            for i, w in enumerate(freq):
                f.write(f'{w:.7f} {spec[i]:.7f}\n')

        return output_file, folder
