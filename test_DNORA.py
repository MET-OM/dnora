#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 11:56:40 2021

@author: janvb
"""
from dnora2 import bnd
import numpy as np
import pandas as pd
start_date = '2019-01-09T00:00' ; end_date = '2019-01-10T00:00'
project_name = 'Sulafjorden'; dgm = 250 ; dbm = 10000 
bnd_points = np.loadtxt(project_name+str(dgm)+'_Boundaries.txt')

# Cut data to the desired points
point_picker = bnd.PPNearestGridPoint(bnd_points)
#point_picker = bnd.PPLegacyPicker(bnd_points) # This is my best implementation of what was in the original DNORA.py

# Read input
#read_boundary_spectra = bnd.InputNORA3(point_picker)
read_boundary_spectra = bnd.InputWAM4(point_picker)
bnd_spec = read_boundary_spectra(start_date, end_date)

# Write output
#write_boundary_spectra = bnd.OutputSWANascii()
#write_boundary_spectra(bnd_spec)

