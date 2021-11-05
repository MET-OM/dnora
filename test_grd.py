#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 14:26:52 2021

@author: janvb
"""
import dnora2.grd as grd


lon_min=5.39; lat_min=62.05; lon_max=6.8; lat_max=62.61
grid = grd.Grid(lon_min, lon_max, lat_min, lat_max, name = 'Sulafjorden250')

grid.set_spacing(dm = 250)

grid.set_boundary(bounN = 1, edges = ['N', 'W']) 

topo_fetcher = grd.TopoEMODNET2018()
grid.import_topo(topo_fetcher)
grid.mesh_grid()

# Change grid points under 2 m depth to land
grid.filter_grid(grd.SetMinDepth(-2, to_land = 0))

# We can check the grid status with print(grid)
print(grid)

# Impose a minimum depth of 5 m on the rest of the sea points
grid.filter_grid(grd.SetMinDepth(-5))

print(grid)

grid.plot(save_fig=True)
grid.plot_mask()

ww3_output = grd.OutputModelWW3()

ww3_output(grid)

grid.write_status() ## This writes the status to a file named after the grid name
grid.write_status(filename = 'another_file.temp') ## We can override the default name like this

print("################# STARTING NEW GRID #########################")
# We can regenerate the grid from above from its output files that wew created by write_output
gridname = 'Sulafjorden250'
grid_regen = grd.regenerate_ww3(gridname)


print("################# STARTING NEW GRID #########################")

grid2 = grd.Grid(lon_min, lon_max, lat_min, lat_max, name = 'Sulafjorden250v2')

# This keeps the grid edge definitions fixed and takes dlon/dlat as close as possible
grid2.set_spacing(dlon=1/240, dlat= 1/480) # dlat = 1/480 is one eight of a nautical mile

grid2.set_boundary(bounN = 1, edges = ['N', 'W']) # SEt every third point to boundary point

print("################# STARTING NEW GRID #########################")

grid3 = grd.Grid(lon_min, lon_max, lat_min, lat_max, name = 'Sulafjorden250v3')

# This uses exactly the given dlon and dlat, and changes the edges of the grid slightly
grid3.set_spacing(dlon=1/240, dlat= 1/480, floating_edge = True) # dlat = 1/480 is one eight of a nautical mile

# We only have the automatically created trivial grid, but lets write it out anyway
ww3_output(grid3, matrix = True) # Write in more human readable format



print("################# STARTING NEW GRID #########################")

grid4 = grd.Grid(lon_min, lon_max, lat_min, lat_max, name = 'Sulafjorden250v4')

# We can also specify the number of grid points
grid4.set_spacing(nx = 291, ny = 249)

# We only have the automatically created trivial grid, but lets write it out anyway
ww3_output(grid4, matrix = True) # Write in more human readable format
