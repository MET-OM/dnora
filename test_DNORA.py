#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 11:56:40 2021

@author: janvb
"""
from dnora2 import bnd
from dnora2 import grd
import numpy as np
#import pandas as pd
start_date = '2019-01-09T00:00' ; end_date = '2019-01-10T00:00'
#start_date = '2021-01-09T00:00' ; end_date = '2021-01-12T00:00'
#start_date = '2020-01-01T00:00' ; end_date = '2020-12-31T23:00'
#project_name = 'Sulafjorden'; dgm = 250 ; dbm = 10000 
#bnd_points = np.loadtxt(project_name+str(dgm)+'_Boundaries.txt')
#bnd_points = np.loadtxt(project_name+str(dgm)+'_1point.txt')

gridname = 'Sulafjorden250'
grid = grd.regenerate_ww3(gridname)

grid.set_boundary(bounN = 25, edges = ['N', 'W']) 


#p_picker = bnd.PPNearestGridPoint()
p_picker = bnd.PPAreaPicker()

#read_boundary_spectra = bnd.InputNORA3()
read_boundary_spectra = bnd.InputWAM4()
bnd_spec, bnd_mask = read_boundary_spectra(start_date, end_date, grid, point_picker = p_picker)

spectral_processor = bnd.TrivialSpectralProcessor(calib_spec = 1)
bnd_spec, bnd_mask = spectral_processor(bnd_spec, bnd_mask)

spectral_processor = bnd.NaNCleanerSpectralProcessor()
bnd_spec, bnd_mask = spectral_processor(bnd_spec, bnd_mask)

# Write output
#write_boundary_spectra = bnd.OutputSWANascii(grid, factor = 1E-4)
write_boundary_spectra = bnd.OutputWW3nc()

write_boundary_spectra(bnd_spec, bnd_mask)


