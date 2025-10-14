import numpy as np
from copy import copy
from abc import ABC, abstractmethod
from typing import Optional
from dnora.utils.spec import (
    interp_spec,
    shift_spec,
    flip_spec,
    check_that_spectra_are_consistent,
)

from dnora.type_manager.spectral_conventions import SpectralConvention


class SpectralProcessor(ABC):
    def __init__(self):
        pass

    def _convention_out(self) -> str:
        """The convention of the spectra returned to the object.

        The conventions to choose from are predetermined and listed below.
        Return None if the convention is not changed (e.g. interpolation or
        multiplication etc.)

        OCEAN:    Oceanic convention
                    Directional vector monotonically increasing.
                    Direction to. North = 0, East = 90.

        MET:      Meteorological convention
                    Directional vector monotonically increasing.
                    Direction from. North = 0, East = 90.

        MATH:     Mathematical convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        MATHVEC:  Mathematical convention in vector
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        WW3:      WAVEWATCH III output convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 0, East = 90.
        """
        return None

    def _convention_in(self) -> str:
        """The convention of the spectra passed to the processor.

        The conventions to choose from are predetermined and listed below.
        Return None if the convention does not matter (e.g. interpolation or
        multiplication etc.)

        OCEAN:    Oceanic convention
                    Directional vector monotonically increasing.
                    Direction to. North = 0, East = 90.

        MET:      Meteorological convention
                    Directional vector monotonically increasing.
                    Direction from. North = 0, East = 90.

        MATH:     Mathematical convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        MATHVEC:  Mathematical convention in vector
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        WW3:      WAVEWATCH III output convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 0, East = 90.
        """
        return None

    @abstractmethod
    def __call__(self, spec, dirs, freq, inds, times, sprm) -> tuple:
        """Processes individual spectra and returns them to object.

        In addition to the spectra, also the direction and frequency
        vectors can be modified. The Boundary object in dnora is
        4 dimensional: [station, time, freq, dirs].

        NB! It is recommended that the processor should be able to also deal
        with 2, 3 and 4 dimensional objects so that it can be called by the user
        to modify a single spectra etc.
        """
        return spec, dirs, freq, inds

    @abstractmethod
    def __str__(self):
        """Describes how the spectral values as processed"""
        pass


class Multiply(SpectralProcessor):
    """Multiplies all spectra with a constant defined at initialization."""

    def __init__(self, calib_spec=1) -> None:
        self.calib_spec = calib_spec
        return

    def __call__(self, spec, dirs, freq, inds, times, spr=None) -> tuple:
        spec_dim = 1 if spr is not None else 2
        check_that_spectra_are_consistent(spec, dirs, freq, expected_dim=spec_dim)
        new_spec = spec * self.calib_spec
        check_that_spectra_are_consistent(new_spec, dirs, freq, expected_dim=spec_dim)

        if spec_dim == 1:
            return new_spec, dirs, freq, inds, times, spr
        else:
            return new_spec, dirs, freq, inds, times

    def __str__(self):
        return f"Multiplying spectral values with {self.calib_spec}"


class CutFrequency(SpectralProcessor):
    """Cuts the spectrum down to a certain frequency range"""

    def __init__(self, freq: tuple):
        self._freq = freq

    def __call__(self, spec, dirs, freq, inds, time, spr=None) -> tuple:
        spec_dim = 1 if spr is not None else 2
        check_that_spectra_are_consistent(spec, dirs, freq, expected_dim=spec_dim)
        mask = np.logical_and(freq >= self._freq[0], freq <= self._freq[-1])
        new_freq = freq[mask]
        if spr is not None:
            new_spr = spr[:, :, mask]
            new_dirs = dirs[:, :, mask]
        else:
            new_dirs = dirs
        new_spec = spec[:, :, mask]
        check_that_spectra_are_consistent(
            new_spec, new_dirs, new_freq, expected_dim=spec_dim
        )

        if spec_dim == 1:
            return new_spec, new_dirs, new_freq, inds, time, new_spr
        else:
            return new_spec, new_dirs, new_freq, inds, time

    def __str__(self):
        return f"Cutting frequency range to {self._freq[0]}-{self._freq[-1]}..."


class RemoveEmpty(SpectralProcessor):
    """Remove all empty spectra."""

    def __init__(self, threshold: float = 0.01) -> None:
        self.threshold = threshold

    def __call__(self, spec, dirs, freq, inds, times, spr=None) -> tuple:
        spec_dim = 1 if spr is not None else 2
        check_that_spectra_are_consistent(spec, dirs, freq, expected_dim=spec_dim)

        mask = np.full(len(inds), True)
        for n in inds:
            if (
                np.max(spec[:, n, ...]) < self.threshold
                or np.isnan(spec[:, n, ...]).any()
            ):
                mask[n] = False

        new_inds = inds[mask]
        new_spec = spec[:, mask, ...]
        if spec_dim == 1:
            new_dirs = dirs[:, mask, ...]
            new_spr = spr[:, mask, ...]

        check_that_spectra_are_consistent(spec, dirs, freq, expected_dim=spec_dim)

        if spec_dim == 1:
            return new_spec, new_dirs, freq, new_inds, times, new_spr
        else:
            return new_spec, dirs, freq, new_inds, times

    def __str__(self):
        return f"Removing spectra with all values less than {self.threshold} or with NaN's..."


class RemoveNanTimes(SpectralProcessor):
    """Remove times when some index has nan-values"""

    def __call__(self, spec, dirs, freq, inds, times, spr=None) -> tuple:
        spec_dim = 1 if spr is not None else 2
        check_that_spectra_are_consistent(spec, dirs, freq, expected_dim=spec_dim)

        mask = np.full(len(times), True)
        for n in range(len(times)):
            if spec_dim == 2:
                if np.isnan(spec[n, :, :, :]).any():
                    mask[n] = False
            elif spec_dim == 1:
                if np.isnan(spec[n, :, :]).any():
                    mask[n] = False

        new_times = times[mask]
        if spec_dim == 2:
            new_spec = spec[mask, :, :, :]
            new_dirs = dirs
        elif spec_dim == 1:
            new_spec = spec[mask, :, :]
            new_dirs = dirs[mask, :, :]
            new_spr = spr[mask, :, :]

        check_that_spectra_are_consistent(
            new_spec, new_dirs, freq, expected_dim=spec_dim
        )

        if spec_dim == 2:
            return new_spec, new_dirs, freq, inds, new_times
        else:
            return new_spec, new_dirs, freq, inds, new_times, new_spr

    def __str__(self):
        return f"Removing times where at least one point has NaN value..."


class ReGridDirs(SpectralProcessor):
    """Interpolates the spectra to have the same resoltuon but to start from
    a certain values, e.g. 0."""

    def __init__(
        self,
        res: Optional[int] = None,
        first_dir: Optional[int] = None,
        nbins: Optional[int] = None,
    ) -> None:
        self.first_dir = copy(first_dir)
        if res is not None and nbins is not None:
            raise ValueError(f"Provide EITHER 'res' or 'nbins'!")
        self.res = copy(res)
        self.nbins = copy(nbins)
        return

    def __call__(self, spec, dirs, freq, inds, times) -> tuple:
        check_that_spectra_are_consistent(spec, dirs, freq, expected_dim=2)

        nbins = self.nbins or len(dirs)
        dD = self.res or 360 / nbins
        dD = int(dD)

        first_dir = self.first_dir or dirs[0]

        new_dirs = np.array(range(0, 360, dD), dtype="float32") + first_dir
        new_spec = np.zeros((spec.shape[0], len(inds), len(freq), len(new_dirs)))
        N1 = spec.shape[0]
        N2 = spec.shape[1]
        for n1 in range(N1):
            for n2 in range(N2):
                m0_orig = np.trapz(np.sum(spec[n1, n2, :], axis=1) * 360 / nbins, freq)
                temp_spec = interp_spec(freq, dirs, spec[n1, n2, :], freq, new_dirs)
                m0_new = np.trapz(np.sum(temp_spec, axis=1) * dD, freq)
                new_spec[n1, n2, :, :] = temp_spec * m0_new / m0_orig
        check_that_spectra_are_consistent(new_spec, new_dirs, freq, expected_dim=2)
        return new_spec, new_dirs, freq, inds, times

    def __str__(self):
        return_string = ""
        if self.first_dir is not None:
            return_string += f"Interpolating spectra to start from {self.first_dir} "
        if self.res is not None:
            return_string += f"Setting resolution to {self.res}"

        return return_string


# class ClearNaN(SpectralProcessor):
#     def __init__(self):
#         pass
#
#     def __call__(self, spec, freq, dirs, time, x, lon, lat, mask):
#         new_spec = copy(spec)
#         new_mask = copy(mask)
#         new_freq = copy(freq)
#         new_dirs = copy(dirs)
#
#         for n in range(len(x)):
#             if np.isnan(spec[:,n,:,:]).any():
#                 msg.info(f"Point {n} ({lon[n]:10.7f}, {lat[n]:10.7f}) contains NaN's. Masking as False.")
#                 new_mask[n] = False
#
#         return new_spec, new_mask, new_freq, new_dirs


class OceanToWW3(SpectralProcessor):
    """Changes all spectra from Oceanic convention to WW3 convention.
    OCEAN:    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.

    WW3:      WAVEWATCH III output convention
                Directional vector of type: [90 80 ... 10 0 350 ... 100]
                Direction to. North = 0, East = 90.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return SpectralConvention.OCEAN

    def _convention_out(self) -> str:
        return SpectralConvention.WW3

    def __call__(self, spec, dirs, freq, inds, times) -> tuple:
        check_that_spectra_are_consistent(spec, dirs, freq, expected_dim=2)
        # Flip direction of the both spectra and directional vector
        spec_flip = flip_spec(spec, dirs)
        D_flip = flip_spec(dirs, dirs)

        # Shift 0 to be at 90
        new_spec = shift_spec(spec_flip, D_flip, -90)
        # Also shift direction
        new_dirs = shift_spec(D_flip, D_flip, -90)

        check_that_spectra_are_consistent(new_spec, new_dirs, freq, expected_dim=2)
        return new_spec, new_dirs, freq, inds, times

    def __str__(self):
        return "Flipping both spectra and direction to mathematical notation."


class WW3ToOcean(SpectralProcessor):
    """Changes all spectra from WW3 convention to Oceanic convention.

    WW3:      WAVEWATCH III output convention
                Directional vector of type: [90 80 ... 10 0 350 ... 100]
                Direction to. North = 0, East = 90.

    OCEAN:    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return SpectralConvention.WW3

    def _convention_out(self) -> str:
        return SpectralConvention.OCEAN

    def __call__(self, spec, dirs, freq, inds, times) -> tuple:
        check_that_spectra_are_consistent(spec, dirs, freq, expected_dim=2)
        # Flip direction of the both spectra and directional vector
        spec_flip = flip_spec(spec, dirs)
        D_flip = flip_spec(dirs, dirs)

        # Shift 0 to be at 90
        new_spec = shift_spec(spec_flip, D_flip, 90)
        # Also shift direction
        new_dirs = shift_spec(D_flip, D_flip, 90)

        check_that_spectra_are_consistent(new_spec, new_dirs, freq, expected_dim=2)

        return new_spec, new_dirs, freq, inds, times

    def __str__(self):
        return "Flipping both spectra and direction to oceanic notation."


class MathToOcean(SpectralProcessor):
    """Changes all spectra from Mathematical convention to Oceanic convention.

    MATH:     Mathematical convention
                Directional vector monotonically increasing.
                Direction to. North = 90, East = 0.

    OCEAN:    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return SpectralConvention.MATH

    def _convention_out(self) -> str:
        return SpectralConvention.OCEAN

    def __call__(self, spec, dirs, freq, inds, times) -> tuple:
        check_that_spectra_are_consistent(spec, dirs, freq, expected_dim=2)
        # Flip direction of the both spectra and directional vector
        spec_flip = flip_spec(spec, dirs)
        D_flip = flip_spec(dirs, dirs)

        # Shift 0 to be at 90
        new_spec = shift_spec(spec_flip, D_flip, -90)
        # Also shift direction
        # new_dirs = shift_spec(D_flip, D_flip, -270)

        check_that_spectra_are_consistent(new_spec, dirs, freq, expected_dim=2)
        return new_spec, dirs, freq, inds, times

    def __str__(self):
        return "Flipping spectra to oceanic notation."


class MathVecToOcean(SpectralProcessor):
    """Changes all spectra from MathVec convention to Oceanic convention.

    MATHVEC:  Mathematical convention in vector
                Directional vector of type: [90 80 ... 10 0 350 ... 100]
                Direction to. North = 90, East = 0.

    OCEAN:    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return SpectralConvention.MATHVEC

    def _convention_out(self) -> str:
        return SpectralConvention.OCEAN

    def __call__(self, spec, dirs, freq, inds, times) -> tuple:
        # Flip direction of the both spectra and directional vector
        spec_flip = flip_spec(spec, dirs)
        D_flip = flip_spec(dirs, dirs)

        # Shift 0 to be at 90
        # new_spec = shift_spec(spec_flip, D_flip, -270)
        # Also shift direction
        new_dirs = shift_spec(D_flip, D_flip, -270)

        # if freq is not None:
        return spec, new_dirs, freq, inds, times
        # else:
        #     return spec, new_dirs

    def __str__(self):
        return "Flipping direction to oceanic notation."


class OceanToMath(SpectralProcessor):
    """Changes all spectra from Oceanic convention to Mathematical convention.

    OCEAN:    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.

    MATH:     Mathematical convention
                Directional vector monotonically increasing.
                Direction to. North = 90, East = 0.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return SpectralConvention.OCEAN

    def _convention_out(self) -> str:
        return SpectralConvention.MATH

    def __call__(self, spec, dirs, freq, inds, times) -> tuple:
        # Flip direction of the both spectra and directional vector
        spec_flip = flip_spec(spec, dirs)
        D_flip = flip_spec(dirs, dirs)

        # Shift 0 to be at 90
        new_spec = shift_spec(spec_flip, D_flip, -90)

        return new_spec, dirs, freq, inds, times

    def __str__(self):
        return "Flipping only spectra to mathematical notation."


class OceanToMathVec(SpectralProcessor):
    """Changes all spectra from Oceanic convention to Mathematical convention.

    OCEAN:    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.

    MATHVEC:     Mathematical convention
                Directional vector of type: [90 80 ... 10 0 350 ... 100]
                Direction to. North = 90, East = 0.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return SpectralConvention.OCEAN

    def _convention_out(self) -> str:
        return SpectralConvention.MATHVEC

    def __call__(self, spec, dirs, freq, inds, times) -> tuple:
        # Flip direction of the both spectra and directional vector
        spec_flip = flip_spec(spec, dirs)
        D_flip = flip_spec(dirs, dirs)

        # Shift 0 to be at 90
        new_dirs = shift_spec(D_flip, D_flip, -90)

        return spec, new_dirs, freq, inds, times

    def __str__(self):
        return "Flipping only directions to mathematical notation."


class MetToOcean(SpectralProcessor):
    """Changes all spectra from Meteorological convention to Ocanic convention.

    OCEAN:    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.

    MET:      Meteorological convention
                Directional vector monotonically increasing.
                Direction from. North = 0, East = 90.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return SpectralConvention.MET

    def _convention_out(self) -> str:
        return SpectralConvention.OCEAN

    def __call__(self, spec, dirs, freq, inds, times) -> tuple:
        new_spec = shift_spec(spec, dirs, 180)

        return new_spec, dirs, freq, inds, times

    def __str__(self):
        return "Shifting spectrum 180 degrees."


class OceanToMet(SpectralProcessor):
    """Changes all spectra from Oceanic convention to Meteorological convention.

    OCEAN:    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.

    MET:      Meteorological convention
                Directional vector monotonically increasing.
                Direction from. North = 0, East = 90.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return SpectralConvention.OCEAN

    def _convention_out(self) -> str:
        return SpectralConvention.MET

    def __call__(self, spec, dirs, freq, inds, times) -> tuple:
        new_spec = shift_spec(spec, dirs, 180)

        return new_spec, dirs, freq, inds, times

    def __str__(self):
        return "Shifting spectrum 180 degrees."


def spectral_processor_for_convention_change(
    current_convention: str, wanted_convention: str
) -> SpectralProcessor:
    """Provides a SpectralProcessor object for the .change_convention()
    method of the Boundary objects.

    The conventions to choose from are predetermined:

    OCEAN:    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.

    MET:      Meteorological convention
                Directional vector monotonically increasing.
                Direction from. North = 0, East = 90.

    MATH:     Mathematical convention
                Directional vector of type: [90 80 ... 10 0 350 ... 100]
                Direction to. North = 90, East = 0.

    MATHVEC:  Mathematical convention in vector
                Directional vector of type: [90 80 ... 10 0 350 ... 100]
                Direction to. North = 90, East = 0.

    WW3:      WAVEWATCH III output convention
                Directional vector of type: [90 80 ... 10 0 350 ... 100]
                Direction to. North = 0, East = 90.
    """

    dict_of_processors = {
        SpectralConvention.OCEAN: {
            SpectralConvention.MET: OceanToMet(),
            SpectralConvention.WW3: OceanToWW3(),
            SpectralConvention.MATH: OceanToMath(),
            SpectralConvention.MATHVEC: OceanToMathVec(),
        },
        SpectralConvention.MET: {
            SpectralConvention.OCEAN: MetToOcean(),
            SpectralConvention.WW3: [MetToOcean(), OceanToWW3()],
            SpectralConvention.MATH: [MetToOcean(), OceanToMath()],
            SpectralConvention.MATHVEC: [MetToOcean(), OceanToMathVec()],
        },
        SpectralConvention.WW3: {
            SpectralConvention.OCEAN: WW3ToOcean(),
            SpectralConvention.MET: [WW3ToOcean(), OceanToMet()],
            SpectralConvention.MATH: [WW3ToOcean(), OceanToMath()],
            SpectralConvention.MATHVEC: [WW3ToOcean(), OceanToMathVec()],
        },
        SpectralConvention.MATH: {
            SpectralConvention.OCEAN: MathToOcean(),
            SpectralConvention.WW3: [MathToOcean(), OceanToWW3()],
            SpectralConvention.MET: [MathToOcean(), OceanToMet()],
            SpectralConvention.MATHVEC: [MathToOcean(), OceanToMathVec()],
        },
        SpectralConvention.MATHVEC: {
            SpectralConvention.OCEAN: MathVecToOcean(),
            SpectralConvention.MET: [MathVecToOcean(), OceanToMet()],
            SpectralConvention.WW3: [MathVecToOcean(), OceanToWW3()],
            SpectralConvention.MATH: [MathVecToOcean(), OceanToMath()],
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
