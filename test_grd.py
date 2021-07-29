#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 14:26:52 2021

@author: janvb
"""
import dnora2.grd as grd

grid = grd.WW3Grid(6.,9.,63.,64.,0.2,0.1,'dummy')
topo_fetcher = grd.TopoEMODNET2018()
grid.import_topo(topo_fetcher)

grid.set_boundary(3) # SEt every third point to boundary point

print(grid)

grid.write_topo()
grid.write_topo(matrix = True) # Write in more human readable format