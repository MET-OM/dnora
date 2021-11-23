import numpy as np
from copy import copy
from abc import ABC, abstractmethod
from .. import msg
#from .bnd_mod import SpectralProcessor
from ..aux import interp_spec, shift_spec, flip_spec

#from . import bnd_mod
class SpectralProcessor(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def get_convention(self) -> str:
        pass

    @abstractmethod
    def __call__(self, spec, freqs, dirs):
        pass

    @abstractmethod
    def __str__(self):
        """Describes how the spectral values as processed"""
        pass

class Multiply(SpectralProcessor):
    def __init__(self, calib_spec = 1):
        self.calib_spec = calib_spec
        return

    def get_convention(self) -> str:
        return None

    def __call__(self, spec, dirs, freq):
        new_spec = copy(spec)*self.calib_spec
        return new_spec, dirs, freq

    def __str__(self):
        return(f"Multiplying spectral values with {self.calib_spec}")

class ReGridDirs(SpectralProcessor):
    def __init__(self, first_dir = 0):
        self.first_dir = copy(first_dir)

        return

    def get_convention(self) -> str:
        return None

    def __call__(self, spec, dirs, freq):

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
            #    Nx
            #Nx = len(self.x())
            # Nt = len(self.time())
            # for x in range(Nx):
            #     for t in range(Nt):


        return new_spec, new_dirs, freq

    def __str__(self):
        return(f"Interpolating spectra to start from {self.first_dir}")

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
    def __init__(self):
        pass

    def get_convention(self) -> str:
        return 'WW3'

    def __call__(self, spec, dirs, freq = None):
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

class WW3ToOcean(SpectralProcessor):
    def __init__(self):
        pass

    def get_convention(self) -> str:
        return 'Ocean'

    def __call__(self, spec, dirs, freq = None):
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


class OceanToMath(SpectralProcessor):
    def __init__(self):
        pass

    def get_convention(self) -> str:
        return 'Math'

    def __call__(self, spec, dirs, freq = None):
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

class MetToOcean(SpectralProcessor):
    def __init__(self):
        pass

    def get_convention(self) -> str:
        return 'Ocean'

    def __call__(self, spec, dirs, freq = None):
        new_spec = shift_spec(spec, dirs, 180)

        if freq is not None:
            return new_spec, dirs, freq
        else:
            return new_spec, dirs

    def __str__(self):
        return("Shifting spectrum 180 degrees. Convention is now Oceanic (direction to).")

class OceanToMet(SpectralProcessor):
    def __init__(self):
        pass

    def get_convention(self) -> str:
        return 'Met'

    def __call__(self, spec, dirs, freq = None):
        if freq is not None:
            new_spec,  new_dirs, new_freq = MetToOcean()(spec, dirs, freq)
            return new_spec, new_dirs, freq
        else:
            new_spec,  new_dirs = MetToOcean()(spec, dirs)
            return new_spec, new_dirs

    def __str__(self):
        return("Shifting spectrum 180 degrees. Convention is not Meteorological (direction from).")


def processor_for_convention_change(current_convention: str, wanted_convention: str) -> SpectralProcessor:
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
