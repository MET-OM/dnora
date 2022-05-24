#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 11:35:29 2022

@author: emiliebyermoen
"""

#%% Example E39BuoyFetcher Sulafjorden

tf = E39BuoyFetcher('C')

available = tf.is_available(('2021/08/18T03:00','2021/08/18T05:00'), ('Hm0','tp','latitude','longitude','WindSpeed'))


if available['time'] is not None: 
    data1 = tf.fetch_data(available['time'], wind_elements=available['wind_elements'])
    
#%% Example E39BuoyFetcher Halsafjorden 18.08.2021 

tf = E39BuoyFetcher('G')

available = tf.is_available(('2021/08/18','2021/08/18'), elements=('Hm0','tp','latitude','longitude'))

if available['time'] is not None:
    data2 = tf.fetch_data(available['time'],wind_elements=available['wind_elements'], wave_elements=available['wave_elements'])


#%% Example E39BuoyFetcher Breisundet 1. feb - 05. mar

tf = E39BuoyFetcher('D')

available = tf.is_available(('2021/02/01','2021/03/05'),elements=('Hm0','hm0a','ff'),lon=(0,40),lat=(50,70))

if available['time'] is not None:
    data3 = tf.fetch_data(available['time'],wind_elements=available['wind_elements'], wave_elements=available['wave_elements'])

#%% Example KystverketBuoyFetcher 
        
kystFetch = KystverketBuoyFetcher('F')

available = kystFetch.is_available(('2019','2019/03'), ['latitude','longitude','Average_Wind_Speed','Average_Wind_Direction','Significant_Wave_Height_Hm0'])

data4 = kystFetch.fetch_data(available['time'], wave_elements=available['wave_elements'], wind_elements=available['wind_elements'])
   
#%% Example KystverketBuoyFetcher

kf = KystverketBuoyFetcher('F')

available = kf.is_available(('2018/12','2018/12'), ('latitude','longitude','Average_Wind_Speed'))

data5 = kystFetch.fetch_data(available['time'], wave_elements=available['wave_elements'], wind_elements=available['wind_elements'])

    
    