import numpy as np
from copy import copy
from abc import ABC, abstractmethod
from typing import Tuple

# Import aux_funcsiliry functions
from .. import msg
from ..aux_funcs import interp_spec, shift_spec, flip_spec

class BoundaryProcessor(ABC):
    def __init__(self):
        pass

    def _convention_out(self) -> str:
        """The convention of the spectra returned to the object.

        The conventions to choose from are predetermined and listed below.
        Return None if the convention is not changed (e.g. interpolation or
        multiplication etc.)

        'Ocean':    Oceanic convention
                    Directional vector monotonically increasing.
                    Direction to. North = 0, East = 90.

        'Met':      Meteorological convention
                    Directional vector monotonically increasing.
                    Direction from. North = 0, East = 90.

        'Math':     Mathematical convention
                    Directional vector monotonically increasing.
                    Direction to. North = 90, East = 0.

        'MathVec':  Mathematical convention in vector
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        'WW3':      WAVEWATCH III output convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 0, East = 90.
        """
        return None

    def _convention_in(self) -> str:
        """The convention of the spectra passed to the processor.

        The conventions to choose from are predetermined and listed below.
        Return None if the convention does not matter (e.g. interpolation or
        multiplication etc.)

        'Ocean':    Oceanic convention
                    Directional vector monotonically increasing.
                    Direction to. North = 0, East = 90.

        'Met':      Meteorological convention
                    Directional vector monotonically increasing.
                    Direction from. North = 0, East = 90.

        'Math':     Mathematical convention
                    Directional vector monotonically increasing.
                    Direction to. North = 90, East = 0.

        'MathVec':  Mathematical convention in vector
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        'WW3':      WAVEWATCH III output convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 0, East = 90.
        """
        return None

    @abstractmethod
    def __call__(self, spec, dirs, freq) -> Tuple:
        """Processes individual spectra and returns them to object.

        In addition to the spectra, also the direction and frequency
        vectors can be modified. The Boundary object in dnora is
        4 dimensional: [time, station, freq, dirs].

        NB! It is recommended that the processor should be able to also deal
        with 2, 3 and 4 dimensional objects so that it can be called by the user
        to modify a single spectra etc.
        """
        return spec, dirs, freq

    @abstractmethod
    def __str__(self):
        """Describes how the spectral values as processed"""
        pass

class Multiply(BoundaryProcessor):
    """Multiplies all spectra with a constant defined at initialization."""

    def __init__(self, calib_spec = 1) -> None:
        self.calib_spec = calib_spec
        return

    def __call__(self, spec, dirs, freq) -> Tuple:
        new_spec = copy(spec)*self.calib_spec
        return new_spec, dirs, freq

    def __str__(self):
        return(f"Multiplying spectral values with {self.calib_spec}")

class ReGridDirs(BoundaryProcessor):
    """Interpolates the spectra to have the same resoltuon but to start from
    a certain values, e.g. 0."""
    def __init__(self, first_dir = 0) -> None:
        self.first_dir = copy(first_dir)

        return

    def __call__(self, spec, dirs, freq) -> Tuple:

        if dirs[0] > 0:
            nbins = len(dirs)
            dD=int(360/nbins)

            new_dirs = np.array(range(0,360,dD), dtype='float32') + self.first_dir


            if len(spec.shape) == 2:
                new_spec = interp_spec(freq, dirs, spec, freq, new_dirs)
            elif len(spec.shape) == 3:
                new_spec = copy(spec)
                N = spec.shape[0]
                for n in range(N):
                    temp_spec = interp_spec(freq, dirs, spec[n,:,:], freq, new_dirs)
                    new_spec[n,:,:] = temp_spec
            elif len(spec.shape) == 4:
                new_spec = copy(spec)
                N1 = spec.shape[0]
                N2 = spec.shape[1]
                for n1 in range(N1):
                    for n2 in range(N2):
                        temp_spec = interp_spec(freq, dirs, spec[n1,n2,:], freq, new_dirs)
                        new_spec[n1,n2,:,:] = temp_spec

        return new_spec, new_dirs, freq

    def __str__(self):
        return(f"Interpolating spectra to start from {self.first_dir}")

# class ClearNaN(BoundaryProcessor):
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

class OceanToWW3(BoundaryProcessor):
    """Changes all spectra from Oceanic convention to WW3 convention.
    'Ocean':    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.

    'WW3':      WAVEWATCH III output convention
                Directional vector of type: [90 80 ... 10 0 350 ... 100]
                Direction to. North = 0, East = 90.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return 'Ocean'

    def _convention_out(self) -> str:
        return 'WW3'

    def __call__(self, spec, dirs, freq = None) -> Tuple:
        # Flip direction of the both spectra and directional vector
        spec_flip = flip_spec(spec, dirs)
        D_flip = flip_spec(dirs, dirs)

        # Shift 0 to be at 90
        new_spec = shift_spec(spec_flip, D_flip, -90)
        # Also shift direction
        new_dirs = shift_spec(D_flip, D_flip, -90)

        if freq is not None:
            return new_spec, new_dirs, freq
        else:
            return new_spec, new_dirs

    def __str__(self):
        return("Flipping both spectra and direction to mathematical notation.")

class WW3ToOcean(BoundaryProcessor):
    """Changes all spectra from WW3 convention to Oceanic convention.

    'WW3':      WAVEWATCH III output convention
                Directional vector of type: [90 80 ... 10 0 350 ... 100]
                Direction to. North = 0, East = 90.

    'Ocean':    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return 'WW3'

    def _convention_out(self) -> str:
        return 'Ocean'

    def __call__(self, spec, dirs, freq = None) -> Tuple:
        # Flip direction of the both spectra and directional vector
        spec_flip = flip_spec(spec, dirs)
        D_flip = flip_spec(dirs, dirs)

        # Shift 0 to be at 90
        new_spec = shift_spec(spec_flip, D_flip, 90)
        # Also shift direction
        new_dirs = shift_spec(D_flip, D_flip, 90)

        if freq is not None:
            return new_spec, new_dirs, freq
        else:
            return new_spec, new_dirs

    def __str__(self):
        return("Flipping both spectra and direction to oceanic notation.")

class MathToOcean(BoundaryProcessor):
    """Changes all spectra from Mathematical convention to Oceanic convention.

    'Math':     Mathematical convention
                Directional vector monotonically increasing.
                Direction to. North = 90, East = 0.

    'Ocean':    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.
    """
    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return 'Math'

    def _convention_out(self) -> str:
        return 'Ocean'

    def __call__(self, spec, dirs, freq = None) -> Tuple:
        # Flip direction of the both spectra and directional vector
        spec_flip = flip_spec(spec, dirs)
        D_flip = flip_spec(dirs, dirs)

        # Shift 0 to be at 90
        new_spec = shift_spec(spec_flip, D_flip, -90)
        # Also shift direction
        #new_dirs = shift_spec(D_flip, D_flip, -270)

        if freq is not None:
            return new_spec, dirs, freq
        else:
            return new_spec, dirs

    def __str__(self):
        return("Flipping spectra to oceanic notation.")

class MathVecToOcean(BoundaryProcessor):
    """Changes all spectra from MathVec convention to Oceanic convention.

    'MathVec':  Mathematical convention in vector
                Directional vector of type: [90 80 ... 10 0 350 ... 100]
                Direction to. North = 90, East = 0.

    'Ocean':    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return 'MathVec'

    def _convention_out(self) -> str:
        return 'Ocean'

    def __call__(self, spec, dirs, freq = None) -> Tuple:
        # Flip direction of the both spectra and directional vector
        spec_flip = flip_spec(spec, dirs)
        D_flip = flip_spec(dirs, dirs)

        # Shift 0 to be at 90
        #new_spec = shift_spec(spec_flip, D_flip, -270)
        # Also shift direction
        new_dirs = shift_spec(D_flip, D_flip, -270)

        if freq is not None:
            return spec, new_dirs, freq
        else:
            return spec, new_dirs

    def __str__(self):
        return("Flipping direction to oceanic notation.")


class OceanToMath(BoundaryProcessor):
    """Changes all spectra from Oceanic convention to Mathematical convention.

    'Ocean':    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.

    'Math':     Mathematical convention
                Directional vector monotonically increasing.
                Direction to. North = 90, East = 0.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return 'Ocean'

    def _convention_out(self) -> str:
        return 'Math'

    def __call__(self, spec, dirs, freq = None) -> Tuple:
        # Flip direction of the both spectra and directional vector
        spec_flip = flip_spec(spec, dirs)
        D_flip = flip_spec(dirs, dirs)

        # Shift 0 to be at 90
        new_spec = shift_spec(spec_flip, D_flip, -90)

        if freq is not None:
            return new_spec, dirs, freq
        else:
            return new_spec, dirs

    def __str__(self):
        return("Flipping only spectra to mathematical notation.")

class OceanToMathVec(BoundaryProcessor):
    """Changes all spectra from Oceanic convention to Mathematical convention.

    'Ocean':    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.

    'MathVec':     Mathematical convention
                Directional vector of type: [90 80 ... 10 0 350 ... 100]
                Direction to. North = 90, East = 0.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return 'Ocean'

    def _convention_out(self) -> str:
        return 'MathVec'

    def __call__(self, spec, dirs, freq = None) -> Tuple:
        # Flip direction of the both spectra and directional vector
        spec_flip = flip_spec(spec, dirs)
        D_flip = flip_spec(dirs, dirs)

        # Shift 0 to be at 90
        new_dirs = shift_spec(D_flip, D_flip, -90)

        if freq is not None:
            return spec, new_dirs, freq
        else:
            return spec, new_dirs

    def __str__(self):
        return("Flipping only directions to mathematical notation.")

class MetToOcean(BoundaryProcessor):
    """Changes all spectra from Meteorological convention to Ocanic convention.

    'Ocean':    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.

    'Met':      Meteorological convention
                Directional vector monotonically increasing.
                Direction from. North = 0, East = 90.
    """
    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return 'Met'

    def _convention_out(self) -> str:
        return 'Ocean'

    def __call__(self, spec, dirs, freq = None) -> Tuple:
        new_spec = shift_spec(spec, dirs, 180)

        if freq is not None:
            return new_spec, dirs, freq
        else:
            return new_spec, dirs

    def __str__(self):
        return("Shifting spectrum 180 degrees.")

class OceanToMet(BoundaryProcessor):
    """Changes all spectra from Oceanic convention to Meteorological convention.

    'Ocean':    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.

    'Met':      Meteorological convention
                Directional vector monotonically increasing.
                Direction from. North = 0, East = 90.
    """

    def __init__(self):
        pass

    def _convention_in(self) -> str:
        return 'Ocean'

    def _convention_out(self) -> str:
        return 'Met'

    def __call__(self, spec, dirs, freq = None) -> Tuple:
        new_spec = shift_spec(spec, dirs, 180)

        if freq is not None:
            return new_spec, dirs, freq
        else:
            return new_spec, dirs

    def __str__(self):
        return("Shifting spectrum 180 degrees.")


def processor_for_convention_change(current_convention: str, wanted_convention: str) -> BoundaryProcessor:
        """Provides a BoundaryProcessor object for the .change_convention()
        method of the Boundary objects.

        The conventions to choose from are predetermined:

        'Ocean':    Oceanic convention
                    Directional vector monotonically increasing.
                    Direction to. North = 0, East = 90.

        'Met':      Meteorological convention
                    Directional vector monotonically increasing.
                    Direction from. North = 0, East = 90.

        'Math':     Mathematical convention
                    Directional vector monotonically increasing.
                    Direction to. North = 90, East = 0.

        'MathVec':  Mathematical convention in vector
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        'WW3':      WAVEWATCH III output convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 0, East = 90.
        """

        dict_of_processors = {
        'Ocean':    {'Met': OceanToMet(),
                    'WW3': OceanToWW3(),
                    'Math': OceanToMath(),
                    'MathVec': OceanToMathVec()},
        'Met':      {'Ocean': MetToOcean(),
                    'WW3': [MetToOcean(), OceanToWW3()],
                    'Math': [MetToOcean(), OceanToMath()],
                    'MathVec': [MetToOcean(), OceanToMathVec()]},
        'WW3':      {'Ocean': WW3ToOcean(),
                    'Met': [WW3ToOcean(), OceanToMet()],
                    'Math': [WW3ToOcean(), OceanToMath()],
                    'MathVec': [WW3ToOcean(), OceanToMathVec()]},
        'Math':     {'Ocean': MathToOcean(),
                    'WW3': [MathToOcean(), OceanToWW3()],
                    'Met': [MathToOcean(), OceanToMet()],
                    'MathVec': [MathToOcean(), OceanToMathVec()]},
        'MathVec':  {'Ocean': MathVecToOcean(),
                    'Met': [MathVecToOcean(), OceanToMet()],
                    'WW3': [MathVecToOcean(), OceanToWW3()],
                    'Math': [MathVecToOcean(), OceanToMath()]}
        }

        if not wanted_convention or (current_convention == wanted_convention):
            return None
        if not current_convention in list(dict_of_processors.keys()):
            raise ValueError (f"Current convention {current_convention} not recognized! (should be {list(dict_of_processors.keys())})")
        elif not wanted_convention in list(dict_of_processors[current_convention].keys()):
            raise ValueError (f"Wanted convention {wanted_convention} not recognized! (should be {list(dict_of_processors.keys())})")
        elif dict_of_processors[current_convention][wanted_convention] is None:
            raise NotImplementedError(f"Can't process conversion {current_convention} >> {wanted_convention} yet!")
        else:
            return dict_of_processors[current_convention][wanted_convention]
