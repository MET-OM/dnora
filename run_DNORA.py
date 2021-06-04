#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from DNORA import  generate_SWAN_grid, generate_input_swn_file_SWAN , generate_input_NORA3spec_to_SWAN, generate_wind_forcing_SWAN,run_SWAN                                      
    
project_name = 'Vestland'; file_emodnet = '/home/konstantinosc/PhD/github/DNORA/bathy/D5_2018.dtm'; dgm = 500 ; dbm = 10000 
lon_min=4.317; lat_min=58.955; lon_max=5.619; lat_max=59.609
start_date = '2001-01-01T03:00' ; end_date = '2001-01-03T06:00'
swan_directory='/home/konstantinosc/Programs/swan4120'

NX, NY = generate_SWAN_grid(project_name=project_name, file_emodnet=file_emodnet, lon_min=lon_min, lat_min=lat_min, lon_max=lon_max, lat_max=lat_max, dgm=dgm, dbm=dbm)
input_file = generate_input_swn_file_SWAN(project_name = project_name,dgm=dgm,swan_directory =swan_directory, start_date = start_date , end_date = end_date, lon_min=lon_min, lat_min=lat_min, lon_max=lon_max, lat_max=lat_max, NX = NX, NY = NY )
generate_input_NORA3spec_to_SWAN(project_name = project_name,dgm=dgm, start_date = start_date , end_date = end_date , nr_spec_interpolate = 0 )
generate_wind_forcing_SWAN(project_name = project_name, dgm=dgm, start_date=start_date, end_date = end_date, lon_min=lon_min, lat_min=lat_min, lon_max=lon_max, lat_max=lat_max, NX = NX , NY = NY)
run_SWAN(input_file=input_file, swan_directory = swan_directory)

