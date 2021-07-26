#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 11:56:40 2021

@author: janvb
"""
from dnora2 import bnd
import numpy as np
start_date = '2019-01-09T00:00' ; end_date = '2019-01-10T00:00'
project_name = 'Sulafjorden'; dgm = 250 ; dbm = 10000 
bnd_points = np.loadtxt(project_name+str(dgm)+'_Boundaries.txt')

# Read input
input_model = bnd.InputNORA3(start_date, end_date)
bnd_in = input_model.read_bnd()

# Cut data to the desired points
s2g = bnd.NearestSpectraToGrid(bnd_in, bnd_points)
bnd_out = s2g.fit_spectra_to_grid()

# Write output
output_model = bnd.OutputWW3nc(bnd_out)
output_model.output_spec()
