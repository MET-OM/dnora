from __future__ import annotations # For TYPE_CHECKING

import numpy as np
from copy import copy
from abc import ABC, abstractmethod
import netCDF4
import re

# Import abstract classes and needed instances of them
from typing import TYPE_CHECKING, Union
if TYPE_CHECKING:
    from ...mdl.mdl_mod import ModelRun
    from ...file_module import FileNames

from ...bnd.conventions import SpectralConvention

from ...aux_funcs import write_monthly_nc_files

class SpectralWriter(ABC):
    """Writes omnidirectional spectra spectra to a certain file format.

    This object is provided to the .export_spectra() method.
    """

    _convention = None

    def _im_silent(self) -> bool:
        """Return False if you want to be responsible for printing out the
        file names."""
        return True

    def _clean_filename(self) -> bool:
        """If this is set to False, then the ModelRun object does not clean
        the filename, and possible placeholders (e.g. #T0) can still be
        present.
        """
        return True

    def convention(self) -> str:
        """Defines in which format the incoming spectra should be.

        The conventions to choose from are predetermined:

        OCEAN:    Oceanic convention
                    Direction to. North = 0, East = 90.

        MET:      Meteorological convention
                    Direction from. North = 0, East = 90.

        MATH:     Mathematical convention
                    Direction to. North = 90, East = 0.
        """
        if isinstance(self._convention, str):
            self._convention = SpectralConvention[self._convention.upper()]
        return self._convention

    @abstractmethod
    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> tuple[str, str]:
        """Write the data from the Spectra object and returns the file and
        folder where data were written."""

class Null(SpectralWriter):
    def convention(self):
        return SpectralConvention.OCEAN

    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs):
        return ''

class DnoraNc(SpectralWriter):
    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> tuple[str, str]:
        output_files = write_monthly_nc_files(model.spectra(), file_object)
        return output_files

class DumpToNc(SpectralWriter):
    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> tuple[str, str]:
        filename = file_object.get_filepath(extension='nc')
        model.spectra().ds().to_netcdf(filename)
        return filename

class REEF3D(SpectralWriter):
    def __call__(self, model: ModelRun, file_object: FileNames, convention: Union[SpectralConvention, str]=SpectralConvention.MET, **kwargs) -> tuple[str, str]:

        # Take first spectra and first time step for now
        x = 0
        t = 0
        filename = file_object.get_filepath()
        spectra = model.spectra()
        self._convention = convention
        spectra._set_convention(self.convention())
        with open(filename, 'w') as f:
            spec = spectra.spec(angular=True)[x,t,:]
            freq = spectra.freq(angular=True)
            for i, w in enumerate(freq):
                f.write(f'{w:.7f} {spec[i]:.7f}\n')

        return filename
