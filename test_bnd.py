#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 11:56:40 2021

@author: janvb
"""
from dnora2 import bnd
from dnora2 import grd
from dnora2 import spec
#import numpy as np
#import pandas as pd
#start_date = '2019-01-09T00:00' ; end_date = '2019-01-10T00:00'
start_date = '2021-01-09T00:00' ; end_date = '2021-01-12T00:00'
#start_date = '2020-01-01T00:00' ; end_date = '2020-12-31T23:00'
#project_name = 'Sulafjorden'; dgm = 250 ; dbm = 10000 
#bnd_points = np.loadtxt(project_name+str(dgm)+'_Boundaries.txt')
#bnd_points = np.loadtxt(project_name+str(dgm)+'_1point.txt')

gridname = 'Sulafjorden250'
grid = grd.regenerate_ww3(gridname)

grid.set_boundary(bounN = 40, edges = ['N', 'W']) 




# =============================================================================
# lon_min=5.39; lat_min=62.05; lon_max=6.8; lat_max=62.61
# grid = grd.Grid(lon_min, lon_max, lat_min, lat_max, name = 'Sulafjorden250')
# 
# grid.set_spacing(dm = 250)
# =============================================================================


boundary = bnd.Boundary(grid)
boundary_fetcher = bnd.BoundaryWAM4()
#boundary_fetcher = bnd.BoundaryNORA3()
#point_picker = bnd.NearestGridPointPicker()
point_picker = bnd.AreaPicker()
boundary.import_boundary(start_date, end_date, boundary_fetcher, point_picker)
boundary.process_spectra([spec.NaNCleaner(), spec.InterpSpectralProcessor(start = 0)])
write_boundary_spectra = bnd.OutputWW3nc()
write_boundary_spectra(boundary)


# =============================================================================
# 
# #p_picker = bnd.PPNearestGridPoint()
# p_picker = bnd.PPAreaPicker()
# 
# #read_boundary_spectra = bnd.InputNORA3()
# #bnd_spec3, bnd_mask3 = read_boundary_spectra(start_date, end_date, grid, point_picker = p_picker)
# read_boundary_spectra = bnd.InputWAM4()
# bnd_spec4, bnd_mask4 = read_boundary_spectra(start_date, end_date, grid, point_picker = p_picker)
# 
# 
# #spectral_processor = bnd.InterpSpectralProcessor(nbins = 24, start = 0)
# spectral_processor = bnd.InterpSpectralProcessor(data_set = bnd_spec4[0]) # Automatic reading of number of directional bins
# bnd_speci, bnd_maski = spectral_processor(bnd_spec4, bnd_mask4)
# # =============================================================================
# # spectral_processor = bnd.TrivialSpectralProcessor(calib_spec = 1)
# # bnd_spec, bnd_mask = spectral_processor(bnd_spec, bnd_mask)
# # 
# spectral_processor = bnd.NaNCleanerSpectralProcessor()
# bnd_speci, bnd_maski = spectral_processor(bnd_speci, bnd_maski)
# # =============================================================================
# 
# # Write output
# # =============================================================================
# # write_boundary_spectra = bnd.OutputSWANascii(grid, factor = 1E-4)
# # write_boundary_spectra(bnd_spec3, bnd_mask3)
# # write_boundary_spectra(bnd_spec4, bnd_mask4)
# # 
# # 
# 
# #write_boundary_spectra(bnd_spec3, bnd_mask3)
# write_boundary_spectra(bnd_speci, bnd_maski)
# #write_boundary_spectra(bnd_spec4, bnd_mask4)
# # 
# # =============================================================================
# 
# 
# 
# =============================================================================
