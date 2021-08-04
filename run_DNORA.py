#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from DNORA import  generate_grid, generate_input_swn_file_SWAN, generate_input_WAM4spec_to_SWAN , generate_input_NORA3spec_to_SWAN, generate_wind_forcing_SWAN,run_SWAN, generate_boundary_input                                 
    
#project_name = 'Sulafjorden'; file_emodnet = '/home/konstantinosc/PhD/github/DNORA/bathy/C5_2018.dtm'; dgm = 250 ; dbm = 10000 
project_name = 'Sulafjorden'; file_emodnet = '/lustre/storeB/project/fou/om/WW3/bathy/emodnet_115m_x_115m/C5_2018.dtm'; dgm = 250 ; dbm = 250
lon_min=5.39; lat_min=62.05; lon_max=6.8; lat_max=62.61
#start_date = '2019-01-09T00:00' ; end_date = '2019-01-10T00:00'
#cds=0.5000E-04

#swan_directory='/home/konstantinosc/Programs/swan4120'

NX, NY = generate_grid(project_name=project_name, file_emodnet=file_emodnet, lon_min=lon_min, lat_min=lat_min, lon_max=lon_max, lat_max=lat_max, dgm=dgm, dbm=dbm, model='WW3')
#input_file = generate_input_swn_file_SWAN(project_name = project_name,dgm=dgm,swan_directory =swan_directory, calib_wind=1, calib_wcap=cds, start_date = start_date , end_date = end_date, lon_min=lon_min, lat_min=lat_min, lon_max=lon_max, lat_max=lat_max, NX = NX, NY = NY )
#generate_boundary_input(project_name, dgm, start_date, end_date, input_model = 'WAM3', output_model = 'SWAN')
#generate_input_NORA3spec_to_SWAN(project_name = project_name,dgm=dgm,calib_spec=1, start_date = start_date , end_date = end_date , nr_spec_interpolate = 0 )
#generate_input_WAM4spec_to_SWAN(project_name = project_name,dgm=dgm,calib_spec=1, start_date = start_date , end_date = end_date , nr_spec_interpolate = 0 )
#generate_wind_forcing_SWAN(project_name = project_name, dgm=dgm, start_date=start_date, end_date = end_date, lon_min=lon_min, lat_min=lat_min, lon_max=lon_max, lat_max=lat_max, NX = NX , NY = NY)
#run_SWAN(input_file=input_file, swan_directory = swan_directory)
