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



#start_time = '2018-11-01T00:00:00' ; end_time = '2019-12-31T23:59:59'
#start_time = '2018-12-31T00:00:00' ; end_time = '2019-01-16T00:00:00'
start_time = '2020-01-06T00:00:00' ; end_time = '2020-01-09T23:00:00'
waveD = param.Parameter()
waveF = param.Parameter()

#param_fetcher = param.ParameterForceFeed(time = time,  wavedata = {"hs": a, "tp": b}, lon = np.array([0]), lat = np.array([0]))
#param_fetcher = param.ParameterFromNc('202101_E39_D_Breisundet_wave.nc', {'Hm0': 'hs', 'tm01': 'tm'}, dump_rest = True)
param_fetcherD = param.ParameterFromThredds('D_Breisundet', {'Hm0': 'Hm0', 'tm01': 'tm01'}, dump_rest = True)
param_fetcherF = param.ParameterFromThredds('F_Vartdalsfjorden', {'Hm0': 'Hm0', 'tm01': 'tm01'}, dump_rest = True)
waveD.import_parameter(start_time, end_time, param_fetcherD)
waveF.import_parameter(start_time, end_time, param_fetcherF)

#waveD.data.to_netcdf('D_Breisundet_201901.nc')
#waveF.data.to_netcdf('F_Vartdalsfjorden_201901.nc')
#wave.plot(['hs', 'mdir'])
#waveD.statistics([param.WaveStatMean(),])
#waveD.statistics([np.nanmean, np.nanstd])

#plt.xlim(np.datetime64('2020-01-01T00:00:00'),np.datetime64('2020-01-10T00:00:00'))
#ax2=plt.plot(wave.data.time.values,wave.data.tp.values)
