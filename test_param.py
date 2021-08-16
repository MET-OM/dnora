#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 10:41:40 2021

@author: janvb
"""
from dnora2 import param
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
time = np.array([0,1,2])
a= np.array([[0.1, 0.2 , 0.3]])
b= np.array([[0.1, 0.2 , 0.3]])

# /thredds/dodsC/obs/buoy-svv-e39/2020/05/202005_E39_F_Vartdalsfjorden_wave.nc


start_time = '2018-11-01T00:00:00' ; end_time = '2019-12-31T23:59:59'
wave = param.Parameter()

#param_fetcher = param.ParameterForceFeed(time = time,  wavedata = {"hs": a, "tp": b}, lon = np.array([0]), lat = np.array([0]))
#param_fetcher = param.ParameterFromNc('202101_E39_D_Breisundet_wave.nc', {'Hm0': 'hs', 'tm01': 'tm'}, dump_rest = True)
param_fetcher = param.ParameterFromThredds('D_Breisundet', {'Hm0': 'hs', 'tm01': 'tm'}, dump_rest = True)
wave.import_parameter(start_time, end_time, param_fetcher)

wave.plot(['hs', 'mdir'])


#plt.xlim(np.datetime64('2020-01-01T00:00:00'),np.datetime64('2020-01-10T00:00:00'))
#ax2=plt.plot(wave.data.time.values,wave.data.tp.values)