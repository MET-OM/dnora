#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 12:05:30 2021

@author: janvb
"""
import sys
sys.path.append("../")
import numpy as np
import dnora2.bnd as bnd

# =============================================================================
# Tutorial on the function that can be used to transform spectral arrays
# The transformations can be used to move between spectral conventions
# These function are a part of the DNORA software package
# -----------------------------------------------------------------------------
# 6.8.2021 
# Jan-Victor Björkqvist & Konstantinos Christakos
# Norwegian Meteorological Institute
# =============================================================================


# =============================================================================
# Functions only assumes that the directional vector is equally spaced
# It doesn't have to contain 0, and it doesn't have to start at 0 or 90 etc.
# This is cruicial for piping the elementary functions presented below
# =============================================================================
shift = 5 


# =============================================================================
# The directional resolution can be anything
# However, note that we can only shift the spectra an even amount of bins
# Anything else would require interpolating, not just manipulating the array
# The functions throw an expection if the condition is violated
# 
# Specifically, we cannot use the function to shift from a 
# [5, 15, 25, ...] grid to a
# [0, 10, 30, ...] grid.
# =============================================================================
dD = 10

rr = range(0,360,dD)
# Directional vector
D = (np.array(rr) + shift) % 360

# Spectra with same values as direction (divided by 10 to avoid confusion)
# Similar values help us track what is mapped where
S = (np.array([rr, rr, rr]) + shift)/10 % 36


# =============================================================================
# Lets shift both the spectra and directional vector 90 degrees and 180 degrees
# If the pair (S,D) was in OCEANIC convention, then:
# The pairs (D_shift180, D) and (D, S_shift180) are now in NAUTICAL convention
# The pair (D_shift180, S_shift180) is still in OCEANIC convention:
#   Every spectral component is still mapped to itself!
# =============================================================================
S_shift90 = bnd.shift_spec(S,D, 90)
D_shift90 = bnd.shift_spec(D,D, 90)

S_shift180 = bnd.shift_spec(S,D, 180)
D_shift180 = bnd.shift_spec(D,D, 180)

# =============================================================================
# Flip the spectra and directional vector (from clockwise to anticlockwise)
# =============================================================================
S_flip = bnd.flip_spec(S,D)
D_flip = bnd.flip_spec(D,D)            

# =============================================================================
# If we shift the flipped spectra 90 degrees we do the "ocean2math" conversion
# Again, if the pair (S,D) was in OCEANIC convention, then:
# The pairs (D_math, S) and (D, S_math) are now in MATHEMATICAL convention
# The pair (D_math, S_math) is still in OCEANIC convention:
#   Every spectral component is still mapped to itself!
# =============================================================================
S_math = bnd.shift_spec(S_flip, D_flip, 90)
D_math = bnd.shift_spec(D_flip, D_flip, 90)


# =============================================================================
# We also have convenient wrappers for the most common transformations:
# =============================================================================
S0 = bnd.ocean_to_math(S, D)
D0 = bnd.ocean_to_math(D, D)

S1 = bnd.ocean_to_naut(S, D)
D1 = bnd.ocean_to_naut(D, D)

S2 = bnd.naut_to_ocean(S, D)
D2 = bnd.naut_to_ocean(D, D)

