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
input_model = bnd.InputNORA3(start_date, end_date)


bnd_in = input_model.read_bnd()

bnd_points = np.loadtxt(project_name+str(dgm)+'_Boundaries.txt')

#output_model = bnd.OutputSWANascii(bnd_in, bnd_points)
output_model = bnd.OutputWW3nc(bnd_in, bnd_points)
#days = input_model.days()
output_model.output_spec()
