#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 15:33:41 2021

@author: janvb
"""
import sys
sys.path.append("../")
from dnora2 import grd, bnd, spec

# =============================================================================
# GRID CREATION
# =============================================================================
# Create a grid for Sulafjorden
lon_min=5.39; lat_min=62.05; lon_max=6.8; lat_max=62.61
grid = grd.Grid(lon_min, lon_max, lat_min, lat_max, name = 'Sulafjorden250')
grid.set_spacing(dm = 250) # Setting grid spacing to approximately 250 metres

# This is still empty, so lets import a topography from EMODNET
grid.import_topo(grd.TopoEMODNET2018()) # Import topography
grid.mesh_grid() # Default mesher bilinear, but can be changed

# Let's define where we want boundary points in the grid
grid.set_boundary(bounN = 1, edges = ['N', 'W']) # Every wet point on north and west edges

# Change grid points under 2 m depth to land
filter_function = grd.SetMinDepth(-2, to_land = 0)
grid.filter_grid(filter_function)

# Impose a minimum depth of 5 m on the rest of the sea points
filter_function = grd.SetMinDepth(-5)
grid.filter_grid(filter_function)
# =============================================================================


# =============================================================================
# TESTING BOUNDARY SPECTRA
# =============================================================================
# Let's check what type of boundary spectra are available for this grid
bnd_spec = bnd.Boundary(grid) # Create a boundary object for this grid

# Define boundary fetchers 
boundary_fetcherWAM4 = bnd.BoundaryWAM4(ignore_nan = True)
boundary_fetcherNORA3 = bnd.BoundaryNORA3()

# How do we want to pick the spectra around the grid / boundary points?
point_picker = bnd.AreaPicker()

# Import boundary only for one instant (faster loading)
start_date = '2019-01-09T00:00' ; end_date = '2019-01-09T00:00'
bnd_spec.import_boundary(start_date, end_date, boundary_fetcherWAM4, point_picker)

# Plot the grid and the boundary points
grid.plot(boundary = bnd_spec)

# We didn't get any points on the western side. Need to read in a bigger area
point_picker = bnd.AreaPicker(expansion_factor = 2.5) # default 1.5
bnd_spec.import_boundary(start_date, end_date, boundary_fetcherWAM4, point_picker)
grid.plot(boundary = bnd_spec)

# Much better!

# What about NORA3?
bnd_spec.import_boundary(start_date, end_date, boundary_fetcherNORA3, point_picker)
grid.plot(boundary = bnd_spec)
# Original area would have probably been ok for NORA3
# =============================================================================


# =============================================================================
# GETTING ACTUAL BOUNDARY SPECTRA
# =============================================================================
# Choose WAM4 because NORA3 is not available for 2021
start_date = '2021-01-09T00:00' ; end_date = '2021-01-10T00:00'
bnd_spec.import_boundary(start_date, end_date, boundary_fetcherWAM4, point_picker)

# We can get basic info about the spectra
bnd_spec.lon()
bnd_spec.lat()
bnd_spec.dirs()
bnd_spec.freq()
bnd_spec.spec(start_time = '2021-01-09T10:00', end_time = '2021-01-09T10:00', x = [0, 1]) # Get two first spectra in list for one isntance

# Lets write this for WAVEWATCH III
# WW3 needs the spectra to start from 0 in directions (not e.g. [5, 15, 25, ...])
bnd_spec.process_spectra([spec.InterpSpectralProcessor(first_dir = 0)])

# Choose WW3 output routine which rotates spectra for WW3 convention
write_boundary_spectra = bnd.OutputWW3nc()
write_boundary_spectra(bnd_spec)

# Don't forget to write the grid also
write_grid = grd.OutputModelWW3()
write_grid(grid)

# =============================================================================
# SAME FOR SWAN
# =============================================================================
# SWAN output function needs information about the grid points we marked as being boundary points
write_boundary_spectra = bnd.OutputSWANascii(grid)
write_boundary_spectra(bnd_spec)

# SWAN grid output not yet implemented, but shouldn't be hard

# =============================================================================
# RELOAD A GRID
# =============================================================================
# You can reload a saved WW3 grid from the files saved above
grid2 = grd.regenerate_ww3('Sulafjorden250')