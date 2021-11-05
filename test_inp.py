#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 16:07:42 2021

@author: janvb
"""
from dnora2 import grd, inp

gridname = 'Sulafjorden250'
grid = grd.regenerate_ww3(gridname)

start_time = '2021-08-17T00:00' ; end_time = '2021-08-20T00:00'
input_writer = inp.SWANInputFile(grid)
a = input_writer(start_time, end_time)
