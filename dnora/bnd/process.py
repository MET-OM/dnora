import numpy as np
from copy import copy
from abc import ABC, abstractmethod
from typing import Tuple

from .. import msg
from ..aux import interp_spec, shift_spec, flip_spec

#from . import bnd_mod
class BoundaryProcessor(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def convention(self) -> str:
        """Return the convention of the spectra returned to the object.

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
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        'WW3':      WAVEWATCH III output convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 0, East = 90.
        """
        return

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

    def convention(self) -> str:
        return None

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

    def convention(self) -> str:
        return None

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

    def convention(self) -> str:
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
        return("Flipping both spectra and direction to mathematical notation. Convention is unchenged.")

class WW3ToOcean(BoundaryProcessor):
    """Changes all spectra from WW3 convention to Oceanic convention.

    'Ocean':    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.

    'WW3':      WAVEWATCH III output convention
                Directional vector of type: [90 80 ... 10 0 350 ... 100]
                Direction to. North = 0, East = 90.
    """
    def __init__(self):
        pass

    def convention(self) -> str:
        return 'Ocean'

    def __call__(self, spec, dirs, freq = None) -> Tuple:
        # Flip direction of the both spectra and directional vector
        spec_flip = flip_spec(spec, dirs)
        D_flip = flip_spec(dirs, dirs)

        # Shift 0 to be at 90
        new_spec = shift_spec(spec_flip, D_flip, -270)
        # Also shift direction
        new_dirs = shift_spec(D_flip, D_flip, -270)

        if freq is not None:
            return new_spec, new_dirs, freq
        else:
            return new_spec, new_dirs

    def __str__(self):
        return("Flipping both spectra and direction to oceanic notation. Convention is unchanged.")


class OceanToMath(BoundaryProcessor):
    """Changes all spectra from Oceanic convention to Mathematical convention.

    'Ocean':    Oceanic convention
                Directional vector monotonically increasing.
                Direction to. North = 0, East = 90.

    'Math':     Mathematical convention
                Directional vector of type: [90 80 ... 10 0 350 ... 100]
                Direction to. North = 90, East = 0.
    """

    def __init__(self):
        pass

    def convention(self) -> str:
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
        return("Flipping only spectra to mathematical notation. Convention is now Mathematical.")

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

    def convention(self) -> str:
        return 'Ocean'

    def __call__(self, spec, dirs, freq = None) -> Tuple:
        new_spec = shift_spec(spec, dirs, 180)

        if freq is not None:
            return new_spec, dirs, freq
        else:
            return new_spec, dirs

    def __str__(self):
        return("Shifting spectrum 180 degrees. Convention is now Oceanic (direction to).")

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

    def convention(self) -> str:
        return 'Met'

    def __call__(self, spec, dirs, freq = None) -> Tuple:
        if freq is not None:
            new_spec,  new_dirs, new_freq = MetToOcean()(spec, dirs, freq)
            return new_spec, new_dirs, freq
        else:
            new_spec,  new_dirs = MetToOcean()(spec, dirs)
            return new_spec, new_dirs

    def __str__(self):
        return("Shifting spectrum 180 degrees. Convention is not Meteorological (direction from).")


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
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        'WW3':      WAVEWATCH III output convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 0, East = 90.
        """

        if not wanted_convention or (current_convention == wanted_convention):
            msg.info(f'Convention ({current_convention}) equals wanted convention({wanted_convention}). Doing nothing.')
            return None
        else:
            if wanted_convention == 'Ocean':
                if current_convention == 'WW3':
                    return WW3ToOcean()
                elif current_convention == 'Met':
                    return MetToOcean()
                else:
                    raise Exception (msg.info(f"Can't process conversion {current_convention} >> {wanted_convention} yet!"))

            elif wanted_convention == 'Met':
                if current_convention == 'Ocean':
                    return OceanToMet()
                else:
                    raise Exception (msg.info(f"Can't process conversion {current_convention} >> {wanted_convention} yet!"))

            elif wanted_convention == 'WW3':
                if current_convention == 'Ocean':
                    return OceanToWW3()
                else:
                    raise Exception (msg.info(f"Can't process conversion {current_convention} >> {wanted_convention} yet!"))

            elif wanted_convention == 'Math':
                if current_convention == 'Ocean':
                    return OceanToMath()
                else:
                    raise Exception (msg.info(f"Can't process conversion {current_convention} >> {wanted_convention} yet!"))
            else:
                raise Exception (msg.info(f"Wanted convention {wanted_convention} not recognized! (should be Ocean/Met/Math/WW3)"))
