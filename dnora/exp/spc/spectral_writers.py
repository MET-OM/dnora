from __future__ import annotations  # For TYPE_CHECKING

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


class SpectralWriter(ABC):
    """Writes omnidirectional spectra spectra to a certain file format.

    This object is provided to the .export_spectra() method.
    """

    _convention = None

    @abstractmethod
    def __call__(
        self, model: ModelRun, file_object: FileNames, **kwargs
    ) -> Union[str, list[str]]:
        """Write the data from the Spectra object and returns the file and
        folder where data were written."""
        pass

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


class REEF3D(SpectralWriter):
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        convention: Union[SpectralConvention, str] = SpectralConvention.MET,
        **kwargs,
    ) -> tuple[str, str]:
        # Take first spectra and first time step for now
        x = 0
        t = 0
        filename = file_object.get_filepath()
        spectra = model.spectra()
        self._convention = convention
        spectra._set_convention(self.convention())
        with open(filename, "w") as f:
            spec = spectra.spec(angular=True)[x, t, :]
            freq = spectra.freq(angular=True)
            for i, w in enumerate(freq):
                f.write(f"{w:.7f} {spec[i]:.7f}\n")

        return filename
