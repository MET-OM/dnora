import numpy as np
from abc import ABC, abstractmethod

# Import aux_funcsiliry functions
from .. import msg

from ..spectral_conventions import SpectralConvention
from ..aux_funcs import check_that_spectra_are_consistent


class SpectralProcessor(ABC):
    def __init__(self):
        pass

    def _convention_out(self) -> str:
        """The convention of the spectra returned to the object.

        The conventions to choose from are predetermined and listed below.
        Return None if the convention is not changed (e.g. interpolation or
        multiplication etc.)

        OCEAN:    Oceanic convention
                    Direction to. North = 0, East = 90.

        MET:      Meteorological convention
                    Direction from. North = 0, East = 90.

        MATH:     Mathematical convention
                    Direction to. North = 90, East = 0.
        """
        return None

    def _convention_in(self) -> str:
        """The convention of the spectra passed to the processor.

        The conventions to choose from are predetermined and listed below.
        Return None if the convention does not matter (e.g. interpolation or
        multiplication etc.)

        OCEAN:    Oceanic convention
                    Direction to. North = 0, East = 90.

        MET:      Meteorological convention
                    Direction from. North = 0, East = 90.

        MATH:     Mathematical convention
                    Direction to. North = 90, East = 0.
        """
        return None

    @abstractmethod
    def __call__(self, spec, dirs, freq, spr, inds) -> tuple:
        """Processes individual spectra and returns them to object.

        In addition to the spectra, also the direction and frequency
        vectors can be modified. The Spectral object in dnora is
        3 dimensional: [station, time, freq].

        NB! It is recommended that the processor should be able to also deal
        with 1, and 2 dimensional objects so that it can be called by the user
        to modify a single spectra etc.
        """
        return spec, dirs, freq, spr, inds

    @abstractmethod
    def __str__(self):
        """Describes how the spectral values as processed"""
        pass


class CutFrequency(SpectralProcessor):
    """Cuts the spectrum down to a certain frequency range"""

    def __init__(self, freq: tuple):
        self._freq = freq

    def __call__(self, spec, dirs, freq, spr, inds) -> tuple:
        check_that_spectra_are_consistent(spec, dirs, freq, expected_dim=1)
        mask = np.logical_and(freq >= self._freq[0], freq <= self._freq[-1])
        new_freq = freq[mask]
        new_dirs = dirs[:, :, mask]
        new_spr = spr[:, :, mask]
        new_spec = spec[:, :, mask]
        check_that_spectra_are_consistent(new_spec, new_dirs, new_freq, expected_dim=1)
        check_that_spectra_are_consistent(new_spec, new_spr, new_freq, expected_dim=1)
        return new_spec, new_dirs, new_freq, new_spr, inds

    def __str__(self):
        return f"Cutting frequency range to {self._freq[0]-self._freq[-1]}..."


class OceanToMet(SpectralProcessor):
    """Changes all directions from Oceanic convention to Meteorological convention.

    OCEAN:    Oceanic convention
                Direction to. North = 0, East = 90.

    MET:      Meteorological convention
                Direction from. North = 0, East = 90.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return SpectralConvention.OCEAN

    def _convention_out(self) -> str:
        return SpectralConvention.MET

    def __call__(self, spec, dirs, freq, spr, inds) -> tuple:
        check_that_spectra_are_consistent(spec, dirs, freq, expected_dim=1)
        new_dirs = np.mod(dirs + 180, 360)
        check_that_spectra_are_consistent(spec, new_dirs, freq, expected_dim=1)
        return spec, new_dirs, freq, spr, inds

    def __str__(self):
        return "Shifting directions 180 degrees."


class MetToOcean(SpectralProcessor):
    """Changes all directions from Meteorological convention to Oceanic convention.

    OCEAN:    Oceanic convention
                Direction to. North = 0, East = 90.

    MET:      Meteorological convention
                Direction from. North = 0, East = 90.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return SpectralConvention.MET

    def _convention_out(self) -> str:
        return SpectralConvention.OCEAN

    def __call__(self, spec, dirs, freq, spr, inds) -> tuple:
        check_that_spectra_are_consistent(spec, dirs, freq, expected_dim=1)
        new_dirs = np.mod(dirs + 180, 360)
        check_that_spectra_are_consistent(spec, new_dirs, freq, expected_dim=1)
        return spec, new_dirs, freq, spr, inds

    def __str__(self):
        return "Shifting directions 180 degrees."


class OceanToMath(SpectralProcessor):
    """Changes all directions from Oceanic convention to Mathematical convention.

    OCEAN:    Oceanic convention
                Direction to. North = 0, East = 90.

    MATH:      Meteorological convention
                Direction from. North = 90, East = 0.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return SpectralConvention.OCEAN

    def _convention_out(self) -> str:
        return SpectralConvention.MATH

    def __call__(self, spec, dirs, freq, spr, inds) -> tuple:
        check_that_spectra_are_consistent(spec, dirs, freq, expected_dim=1)
        new_dirs = np.mod(-dirs + 90, 360)
        check_that_spectra_are_consistent(spec, new_dirs, freq, expected_dim=1)
        return spec, new_dirs, freq, spr, inds

    def __str__(self):
        return "Shifting 90 degrees and changing to anti-clockwise."


class MathToOcean(SpectralProcessor):
    """Changes all directions from Oceanic convention to Mathematical convention.

    OCEAN:    Oceanic convention
                Direction to. North = 0, East = 90.

    MATH:      Meteorological convention
                Direction from. North = 90, East = 0.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return SpectralConvention.MATH

    def _convention_out(self) -> str:
        return SpectralConvention.OCEAN

    def __call__(self, spec, dirs, freq, spr, inds) -> tuple:
        check_that_spectra_are_consistent(spec, dirs, freq, expected_dim=1)
        new_dirs = np.mod(-dirs + 90, 360)
        check_that_spectra_are_consistent(spec, new_dirs, freq, expected_dim=1)
        return spec, new_dirs, freq, spr, inds

    def __str__(self):
        return "Shifting 90 degrees and changing to clockwise."


def spectral_processor_for_convention_change(
    current_convention: str, wanted_convention: str
) -> SpectralProcessor:
    """Provides a BoundaryProcessor object for the .change_convention()
    method of the Boundary objects.

    The conventions to choose from are predetermined:

    OCEAN:    Oceanic convention
                Direction to. North = 0, East = 90.

    MET:      Meteorological convention
                Direction from. North = 0, East = 90.

    MATH:     Mathematical convention
                Direction to. North = 90, East = 0.
    """

    dict_of_processors = {
        SpectralConvention.OCEAN: {
            SpectralConvention.MET: OceanToMet(),
            SpectralConvention.MATH: OceanToMath(),
        },
        SpectralConvention.MET: {
            SpectralConvention.OCEAN: MetToOcean(),
            SpectralConvention.MATH: [MetToOcean(), OceanToMath()],
        },
        SpectralConvention.MATH: {
            SpectralConvention.OCEAN: MathToOcean(),
            SpectralConvention.MET: [MathToOcean(), OceanToMet()],
        },
    }

    if not wanted_convention or (current_convention == wanted_convention):
        return None
    if not current_convention in list(dict_of_processors.keys()):
        raise ValueError(
            f"Current convention {current_convention} not recognized! (should be {list(dict_of_processors.keys())})"
        )
    elif not wanted_convention in list(dict_of_processors[current_convention].keys()):
        raise ValueError(
            f"Wanted convention {wanted_convention} not recognized! (should be {list(dict_of_processors.keys())})"
        )
    elif dict_of_processors[current_convention][wanted_convention] is None:
        raise NotImplementedError(
            f"Can't process conversion {current_convention} >> {wanted_convention} yet!"
        )
    else:
        return dict_of_processors[current_convention][wanted_convention]
