from __future__ import annotations # For TYPE_CHECKING

import numpy as np
from copy import copy
from abc import ABC, abstractmethod
import netCDF4
import re

# Import abstract classes and needed instances of them
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .spc_mod import Spectra

class SpectralWriter(ABC):
    """Writes omnidirectional spectra spectra to a certain file format.

    This object is provided to the .export_spectra() method.
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
    def _convention(self) -> str:
        """Defines in which format the incoming spectra should be.

        The conventions to choose from are predetermined:

        'Ocean':    Oceanic convention (direction to)

        'Met':      Meteorological convention (direction from)
        """

    @abstractmethod
    def __call__(self, spectra: Spectra, filename: str) -> tuple[str, str]:
        """Write the data from the Spectra object and returns the file and
        folder where data were written."""

class DumpToNc(SpectralWriter):
    def __init__(self, convention: str='Met') -> None:
        self._convention_in = convention
        return

    def _extension(self) -> str:
        return 'nc'

    def _convention(self) -> str:
        """Convention of spectra"""
        return self._convention_in

    def __call__(self, spectra: Spectra, filename: str) -> tuple[str, str]:

        spectra.data.to_netcdf(filename)

        return filename

class REEF3D(SpectralWriter):
    def __init__(self, convention: str='Met') -> None:
        self._convention_in = convention
        return

    def _convention(self) -> str:
        """Convention of spectra"""
        return self._convention_in

    def _extension(self) -> str:
        return 'dat'

    def __call__(self, spectra: Spectra, filename: str) -> tuple[str, str]:

        # Take first spectra and first time step for now
        x = 0
        t = 0
        with open(filename, 'w') as f:
            spec = spectra.data.spec.values[x,t,:]
            freq = spectra.data.freq.values
            freq = freq*2*np.pi
            spec = spec/2/np.pi
            for i, w in enumerate(freq):
                f.write(f'{w:.7f} {spec[i]:.7f}\n')

        return filename
