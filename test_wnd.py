#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 11:15:48 2021

@author: janvb
"""
from dnora2 import grd, wnd

gridname = 'Sulafjorden250'
grid = grd.regenerate_ww3(gridname)

#start_time = '2018-12-29T00:00' ; end_time = '2019-01-16T00:00'
start_time = '2021-08-17T00:00' ; end_time = '2021-08-20T00:00'
#start_time = '2018-11-01T00:00:00' ; end_time = '2018-11-01T23:00:00'

forcing = wnd.Forcing(grid, name = 'MEPS')
forcing_fetcher = wnd.ForcingMEPS(prefix = 'det', hours_per_file=73, lead_time=0)

forcing.import_forcing(start_time,end_time, forcing_fetcher)

write_output = wnd.DumpToNc()
write_output(forcing)