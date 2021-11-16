import numpy as np
from copy import copy

from .. import msg

#from ..spec import ocean_to_math, naut_to_ocean, shift_spec
from ..aux import interp_spec, shift_spec, flip_spec

from .bnd_mod import SpectralProcessor

# SpectralProcessors used as defaults in the Boundary object methods
from .bnd_mod import Multiply



class ReGridDirs(SpectralProcessor):
    def __init__(self, first_dir = 0):
        self.first_dir = copy(first_dir)

        return

    def __call__(self, spec, dirs, freq):

        if dirs[0] > 0:
            nbins = len(dirs)
            dD=int(360/nbins)

            new_dirs = np.array(range(0,360,dD), dtype='float32') + self.first_dir

            #msg.info(f"Interpolating spectra to directional grid {new_dirs[0]:.0f}:{dD}:{new_dirs[-1]:.0f}")

            new_spec = interp_spec(freq, dirs, spec, freq, new_dirs)

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

    def __call__(self, spec, dirs, freq = None):
        if freq is not None:
            new_spec,  new_dirs, new_freq = MetToOcean()(spec, dirs, freq)
            return new_spec, new_dirs, freq
        else:
            new_spec,  new_dirs = MetToOcean()(spec, dirs)
            return new_spec, new_dirs

    def __str__(self):
        return("Shifting spectrum 180 degrees. Convention is not Meteorological (direction from).")
