#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 11:56:40 2021

@author: janvb
"""
from dnora2 import bnd
import numpy as np
#import pandas as pd
#start_date = '2019-01-09T00:00' ; end_date = '2019-01-10T00:00'
start_date = '2021-01-09T00:00' ; end_date = '2021-01-12T00:00'
#start_date = '2020-01-01T00:00' ; end_date = '2020-12-31T23:00'
project_name = 'Sulafjorden'; dgm = 250 ; dbm = 10000 
bnd_points = np.loadtxt(project_name+str(dgm)+'_Boundaries.txt')
#bnd_points = np.loadtxt(project_name+str(dgm)+'_1point.txt')

# Cut data to the desired points
point_picker = bnd.PPNearestGridPoint(bnd_points)
#point_picker = bnd.PPLegacyPicker(bnd_points) # This is my best implementation of what was in the original DNORA.py

# Read input
#read_boundary_spectra3 = bnd.InputNORA3(point_picker)
read_boundary_spectra4 = bnd.InputWAM4(point_picker)
#bnd_spec3 = read_boundary_spectra3(start_date, end_date)
bnd_spec4, bnd_mask = read_boundary_spectra4(start_date, end_date)

# Write output
#write_boundary_spectra = bnd.OutputSWANascii(factor = 1E-4, calib_spec = 1)
write_boundary_spectra = bnd.OutputWW3nc(calib_spec = 1)

#write_boundary_spectra(bnd_spec3, bnd_points)
write_boundary_spectra(bnd_spec4, bnd_mask)

