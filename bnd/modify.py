from abc import ABC, abstractmethod
from copy import copy
import numpy as np
from .. import msg
from ..spec import interp_spec, ocean_to_math, naut_to_ocean

class SpectralProcessor(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, bnd_in, bnd_mask):
        pass


class Trivial(SpectralProcessor):
    def __init__(self, calib_spec = 1):
        self.calib_spec = calib_spec
        return

    def __call__(self, spec, freq, dirs, time, x, lon, lat, mask):
        new_spec = copy(spec)*self.calib_spec
        new_mask = copy(mask)
        new_freq = copy(freq)
        new_dirs = copy(dirs)
        return new_spec, new_mask, new_freq, new_dirs

class Interp(SpectralProcessor):
    def __init__(self, first_dir = 0):
        self.first_dir = copy(first_dir)

        return

    def __call__(self, spec, freq, dirs, time, x, lon, lat, mask):
        new_spec = copy(spec)
        new_mask = copy(mask)
        new_freq = copy(freq)
        new_dirs = copy(dirs)

        if dirs[0] > 0:
            nbins = len(dirs)
            dD=int(360/nbins)


            new_dirs = np.array(range(0,360,dD), dtype='float32') + self.first_dir

            msg.info(f"Interpolating spectra to directional grid {new_dirs[0]:.0f}:{dD}:{new_dirs[-1]:.0f}")

            for n in range(len(x)):
                for k in range(len(time)):
                    new_spec[k,n,:,:] = interp_spec(freq, dirs, spec[k,n,:,:], new_freq, new_dirs)

        return new_spec, new_mask, new_freq, new_dirs

class NaNCleaner(SpectralProcessor):
    def __init__(self):
        pass

    def __call__(self, spec, freq, dirs, time, x, lon, lat, mask):
        new_spec = copy(spec)
        new_mask = copy(mask)
        new_freq = copy(freq)
        new_dirs = copy(dirs)

        for n in range(len(x)):
            if np.isnan(spec[:,n,:,:]).any():
                msg.info(f"Point {n} ({lon[n]:10.7f}, {lat[n]:10.7f}) contains NaN's. Masking as False.")
                new_mask[n] = False

        return new_spec, new_mask, new_freq, new_dirs

class OceanToWW3(SpectralProcessor):
    def __init__(self):
        pass

    def __call__(self, spec, freq, dirs, time, x, lon, lat, mask):
        new_spec = copy(spec)
        new_mask = copy(mask)
        new_dirs = copy(dirs)
        new_freq = copy(freq)

        for n in range(len(x)):
            for k in range(len(time)):
                new_spec[k,n,:,:] = ocean_to_math(spec[k,n,:,:], dirs)

        new_dirs = ocean_to_math(dirs, dirs)
        return new_spec, new_mask, new_freq, new_dirs

class NautToOcean(SpectralProcessor):
    def __init__(self):
        pass

    def __call__(self, spec, freq, dirs, time, x, lon, lat, mask):
        new_spec = copy(spec)
        new_mask = copy(mask)
        new_dirs = copy(dirs)
        new_freq = copy(freq)

        for n in range(len(x)):
            for k in range(len(time)):
                new_spec[k,n,:,:] = naut_to_ocean(spec[k,n,:,:], dirs)

        return new_spec, new_mask, new_freq, new_dirs

class OceanToNaut(SpectralProcessor):
    def __init__(self):
        pass
    def __call__(self, spec, freq, dirs, time, x, lon, lat, mask):
        new_spec, new_mask, new_freq, new_dirs = NautToOcean(spec, freq, dirs, time, x, lon, lat, mask)
        return new_spec, new_mask, new_freq, new_dirs
