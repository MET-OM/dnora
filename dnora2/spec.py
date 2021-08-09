#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 15:55:06 2021

@author: janvb
"""
import numpy as np
from scipy import interpolate
from statistics import mode
from copy import copy
from dnora2 import msg
from abc import ABC, abstractmethod

# =============================================================================
# STAND ALONE FUNCTIONS
# =============================================================================
def flip_spec(spec,D):
    # This check enables us to flip directions with flip_spec(D,D)
    
    if len(spec.shape) == 1:
        flipping_dir = True
        spec = np.array([spec])
    else:
        flipping_dir = False
    spec_flip = np.zeros(spec.shape)

    ind = np.arange(0,len(D), dtype='int')
    dD = np.diff(D).mean()
    steps = D/dD # How many delta-D from 0
    
    ind_flip = ((ind - 2*steps).astype(int) + len(D)) % len(D)
    
    spec_flip=spec[:, list(ind_flip)]
    
    if flipping_dir:
        spec_flip = spec_flip[0]
    return spec_flip


def shift_spec(spec, D, shift = 0):
    # This check enables us to flip directions with flip_spec(D,D)
    if len(spec.shape) == 1:
        shifting_dir = True
        spec = np.array([spec])
    else:
        shifting_dir = False
    spec_shift = np.zeros(spec.shape)

    D = np.round(D)
    ind = np.arange(0,len(D), dtype='int')
    dD = mode(abs(np.diff(D)))
    
    if not (shift/dD).is_integer():
        print('aa')
        raise Exception ('Shift needs to be multiple of frequency resolution! Otherwise interpolation would be needed.')
      
    ind_flip = ((ind + int(shift/dD)).astype(int) + len(D)) % len(D)
    
    spec_shift=spec[:, list(ind_flip)]
    if shifting_dir:
        spec_shift = spec_shift[0]
    return spec_shift
    

def ocean_to_naut(oceanspec, D):
    """Convert spectrum in nautical convention (0 north, 90 east, direction from) to oceanic convention (0 north, 90 east, direction to)"""
    nautspec = shift_spec(oceanspec,D, 180)

    return nautspec


def naut_to_ocean(nautspec, D): # Just defined separately to not make for confusing code
    """Convert spectrum in oceanic convention (0 north, 90 east, direction to) to nautical convention (0 north, 90 east, direction from)"""    
    return ocean_to_naut(nautspec, D)


def ocean_to_math(oceanspec, D):
    """Convert spectrum in oceanic convention (0 north, 90 east, direction to) to mathematical convention (90 north, 0 east, direction to)"""
    
    # Flip direction
    spec_flip = flip_spec(oceanspec, D)
    D_flip = flip_spec(D,D)            

    # Shift 0 to be at 90    
    mathspec = shift_spec(spec_flip, D_flip, -90)

    return mathspec

def interp_spec(f, D, S, fi, Di):
    Sleft = S
    Sright = S
    Dleft = -D[::-1] 
    Dright = D + 360
    
    bigS = np.concatenate((Sleft, S, Sright),axis=1)
    bigD = np.concatenate((Dleft, D, Dright))
        
    Finterpolator = interpolate.RectBivariateSpline(f, bigD, bigS, kx=1, ky=1, s=0)
    
    Si = Finterpolator(fi,Di)
    
    return Si
# =============================================================================


# =============================================================================
# SPECTRAL PROCESSOR CLASSES FOR PROCESSING SPECTA OF BOUNDARY OBJECT
# =============================================================================

class SpectralProcessor(ABC):
    def __init__(self):
        pass
    
    @abstractmethod
    def __call__(self, bnd_in, bnd_mask):
        pass


class TrivialSpectralProcessor(SpectralProcessor):
    def __init__(self, calib_spec = 1):
        self.calib_spec = calib_spec
        return
    
    def __call__(self, spec, freq, dirs, time, x, lon, lat, mask):
        new_spec = copy(spec)*self.calib_spec
        new_mask = copy(mask)
        new_freq = copy(freq)
        new_dirs = copy(dirs)
        return new_spec, new_mask, new_freq, new_dirs

class InterpSpectralProcessor(SpectralProcessor):
    def __init__(self, start = 0):
        self.start = copy(start)
             
        return
    
    def __call__(self, spec, freq, dirs, time, x, lon, lat, mask):
        new_spec = copy(spec)
        new_mask = copy(mask)
        new_freq = copy(freq)
        
        nbins = len(dirs)
        dD=int(360/nbins)
        #if self.start > dD:
        #    msg.info("First bin {start} is defined as larger than the directional resolution {dD}. This might spell trouble!")
        
        new_dirs = np.array(range(0,360,dD), dtype='float32') + self.start
        
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
    
class NautToWW3(SpectralProcessor):
    def __init__(self, calib_spec = 1):
        self.calib_spec = calib_spec
        return
    
    def __call__(self, spec, freq, dirs, time, x, lon, lat, mask):
        new_spec = copy(spec)
        new_mask = copy(mask)
        new_dirs = copy(dirs)
        new_freq = copy(freq)
        
        for n in range(len(x)):
            for k in range(len(time)):
                S = naut_to_ocean(new_spec[k,n,:,:], dirs)
                new_spec[k,n,:,:] = ocean_to_math(S, dirs)
                
        new_dirs = ocean_to_math(dirs, dirs)
        return new_spec, new_mask, new_freq, new_dirs
# =============================================================================