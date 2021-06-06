#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#created by: Konstantinos Christakos 
##########################################################################
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from scipy.interpolate import griddata
from subprocess import call, Popen
import pandas as pd
import os



def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def distance_2points(lat1,lon1,lat2,lon2):
    R = 6371.0
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    distance = R * c # in km
    return distance


def generate_SWAN_grid(project_name, file_emodnet, lon_min, lat_min, lon_max, lat_max, dgm, dbm):
    f_min = 0.04
    lonlim = [lon_min, lon_max] 
    latlim = [lat_min, lat_max]
    print('Generate  grid:',project_name+str(dgm))
    ###########################################################################
    distance_x = distance_2points(latlim[0],lonlim[0],latlim[0],lonlim[1])
    distance_y = distance_2points(latlim[0],lonlim[0],latlim[1],lonlim[0])
    nr_points_x = distance_x*1000/dgm
    nr_points_y = distance_y*1000/dgm
    #DX:dlon, DY:dlat,
    DX = np.diff(lonlim)[0]/nr_points_x
    DY = np.diff(latlim)[0]/nr_points_y 
    #Boundaries
    XUR=max(lonlim); YUR=max(latlim);
    XLL=min(lonlim); YLL=min(latlim); 
    YUL=YUR; XLR=XUR
    #
    #Grid size - number of grid positions
    NX=np.round((XLR-XLL)/DX); 
    NY=np.round((YUL-YLL)/DY);
    #
    #Grid lat/lon 
    #newgrid_lon=np.arange(XLL,XLL+((NX-1)*DX),DX)
    #newgrid_lat=np.arange(YLL,YLL+((NY-1)*DY),DY)  
    newgrid_lon=np.linspace(XLL,XUR,int(NX))
    newgrid_lat=np.linspace(YLL,YUR,int(NY))
    newlon,newlat=np.meshgrid(newgrid_lon,newgrid_lat)
    del newgrid_lon,newgrid_lat
    ###########################################################################
    ##    Bathymetry
    ###########################################################################
    print('Loading Emodnet data....')
    ds = xr.open_dataset(file_emodnet)
    lon0,lat0=np.meshgrid(ds.COLUMNS.values,ds.LINES.values)
    M = np.column_stack((ds.DEPTH.values.ravel(), lon0.ravel(),lat0.ravel()))
    # remove data outside the domain of interest
    M = M[M[:,1] > XLL]
    M = M[M[:,1] < XUR]
    M = M[M[:,2] > YLL]
    M = M[M[:,2] < YUR]
    M[M[:,0] > 0] = 0 # remove depths greater than 0 
    print('Generating Grid Bathymetry....')
    topografi=griddata(M[:,1:3],M[:,0],(newlon,newlat), method='nearest')
    topografi[np.isnan(topografi)] = 32767  # replace nan with 32767, land points
    del M, ds, lon0, lat0
    mask_map = np.copy(topografi)
    mask_map[mask_map >= 0] = 0 # land points
    mask_map[mask_map < 0] = 1 # sea points
    
    ############################################################################
    print('Estimate Boundaries....')
    bounN =int(dbm/dgm) #subset of bondary spectra  
    # --------- North boundary ----------
    maskN = mask_map[-1,::bounN]
    maskN[maskN > 0]=2 
    # --------- East boundary ----------
    #maskE = mask_map[::bounN,-1][:-1]
    #maskE[maskE > 0 ] = 2
    # --------- West boundary ----------
    maskW = mask_map[::bounN,0][:-1]
    maskW[maskW > 0 ] = 2
    # --------- South boundary ----------
    maskS = mask_map[0,::bounN]
    maskS[maskS > 0 ] = 2
    # Estimate Active Boundaries (=2 in ww3)
    mask_flat = mask_map.ravel()
    lonlat_flat = np.column_stack((newlon.ravel(),newlat.ravel()))
    BOUND = np.column_stack((lonlat_flat,mask_flat))
    BOUND = BOUND[BOUND[:,2] == 2]
    ##############################################################################

    print('Plotting topografi...')
    levels = np.linspace(np.min(topografi), 0, 100, endpoint=True)
    plt.contourf(newlon,newlat,topografi,levels)
    plt.colorbar()
    plt.savefig(project_name+str(dgm)+'_depth.pdf',dpi=300)
    plt.clf()
    
    print('Plotting mask_map...')
    plt.contourf(newlon,newlat,mask_map)
    plt.plot(BOUND[:,0],BOUND[:,1],'r*')
    plt.colorbar()
    plt.savefig(project_name+str(dgm)+'_mask.pdf')
    plt.clf()
    
    print('Create files for regular grid')
    #boundaries,bathy.txt,mapsta.txt, lat.txt and lon.txt
    np.savetxt(project_name+str(dgm)+'_Boundaries.txt', BOUND[:,::-1][:,1:3], delimiter='\t',fmt='%1.6f')
    np.savetxt(project_name+str(dgm)+'_SWAN.bot', topografi, delimiter='\t',fmt='%1.0f')
    #######################################################################################################
    dt_xy = np.round(123766*DX*np.cos(np.radians(YUR))*f_min)
    dt_global = 2*dt_xy # 2 or 3 times dt_xy
    dt_k = 0.5*dt_global
    dt_s = 15 # 15s is the  minimum allowed by the model 
    
    print('Dt_global:',dt_global,'Dt_xy:',dt_xy,'Dt_k:',dt_k,'Dt_s:',dt_s)
    ######################################################################################
    NX = np.shape(topografi)[1]-1
    NY = np.shape(topografi)[0]-1
    return NX, NY 

def generate_input_swn_file_SWAN(project_name,dgm,swan_directory,calib_wind, calib_wcap, start_date , end_date, lon_min, lat_min, lon_max, lat_max, NX, NY):
    path_forcing = os.getcwd() +'/' # path for directory where forcing and boundaries are saved, here it is used the current directory
    DATE_START = start_date.replace('-','').replace('T','.').replace(':','')+'00'
    DATE_END   = end_date.replace('-','').replace('T','.').replace(':','')+'00'
    delta_X = np.round(np.abs(lon_max - lon_min),5)
    delta_Y = np.round(np.abs(lat_max - lat_min),5)
    factor_wind = calib_wind*0.001
    input_file = swan_directory +'/input_'+DATE_START.split('.')[0]+'_'+project_name+str(dgm)+'.swn'
    with open(input_file, 'w') as file_out:
        file_out.write('$************************HEADING************************\n')
        file_out.write('$ \n')
        file_out.write(' PROJ \'' +project_name+str(dgm)+ '\' \'T24\' \n')
        file_out.write('$ \n')
        file_out.write('$ Topography - Emodnet 2018 \n')
        file_out.write('$ Time of setup: 2021 may \n')
        file_out.write('$ \n')
        file_out.write('$*******************MODEL INPUT*************************\n')
        file_out.write('$ \n')
        file_out.write('SET NAUT \n')
        file_out.write('$ \n')
        file_out.write('MODE NONSTATIONARY TWOD \n')
        file_out.write('COORD SPHE CCM \n')
        file_out.write('CGRID '+str(lon_min)+' '+str(lat_min)+' 0. '+str(delta_X)+' '+str(delta_Y)+' '+str(NX)+' '+str(NY)+' CIRCLE 36 0.04 1.0 31 \n')
        file_out.write('$ \n')
        file_out.write('INPGRID BOTTOM '  +str(lon_min)+' '+str(lat_min)+' 0. '+str(NX)+' '+str(NY)+' '+ str((delta_X/NX).round(4)) +' '+ str((delta_Y/NY).round(4)) +' EXC 32767\n')
        file_out.write('READINP BOTTOM -1 \''+path_forcing+project_name+str(dgm)+'_SWAN.bot\' 3 0 FREE \n')
        file_out.write('$ \n')
        file_out.write('BOU NEST \''+path_forcing+project_name+str(dgm)+'_spec_'+DATE_START.split('.')[0]+'_'+DATE_END.split('.')[0]+'.asc\' OPEN \n')
        file_out.write('$ \n')
        file_out.write('INPGRID WIND '+str(lon_min)+' '+str(lat_min)+' 0. '+str(NX)+' '+str(NY)+' '+str((delta_X/NX).round(4)) +' '+str((delta_Y/NY).round(4)) +' NONSTATIONARY '+ DATE_START +' 1 HR ' + DATE_END +'\n')
        file_out.write('READINP WIND '+str(factor_wind)+'  \''+path_forcing+project_name+str(dgm)+'_wind_'+DATE_START.split('.')[0]+'_'+DATE_END.split('.')[0]+'.asc\' 3 0 0 1 FREE \n')
        file_out.write('$ \n')
        file_out.write('GEN3 WESTH cds2='+str(calib_wcap) +'\n')
        file_out.write('FRICTION JON 0.067 \n')
        #file_out.write('OFF QUAD \n')
        file_out.write('PROP BSBT \n')
        file_out.write('NUM ACCUR NONST 1 \n')
        file_out.write('$ \n')
        file_out.write('$*******************************************************\n')
        file_out.write('$ Generate block-output \n')
        file_out.write('BLOCK \'COMPGRID\' HEAD \''+path_forcing+project_name+str(dgm)+'_'+DATE_START.split('.')[0]+'.nc\' & \n')
        file_out.write('LAY 1 HSIGN RTP TPS PDIR TM01 DIR DSPR WIND DEP OUTPUT '+ DATE_START +' 1 HR \n')
        file_out.write('$ \n')
        file_out.write('COMPUTE '+DATE_START+' 10 MIN '+ DATE_END + '\n')
        file_out.write('STOP \n')
    return input_file


def generate_input_NORA3spec_to_SWAN(project_name, dgm, calib_spec, start_date, end_date, nr_spec_interpolate):    
    factor = 1E-4
    points = np.loadtxt(project_name+str(dgm)+'_Boundaries.txt')
    days = pd.date_range(start=start_date.split('T')[0], end=end_date.split('T')[0], freq='D')
    url = 'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/'+days[0].strftime('%Y') +'/'+days[0].strftime('%m')+'/SPC'+days[0].strftime('%Y%m%d')+'00.nc'
    if len(days)>1:
        data = xr.open_dataset(url).sel(time=slice(start_date, start_date.split('T')[0]+'T23:00'))
    else:
        data = xr.open_dataset(url).sel(time=slice(start_date, end_date))
    delth = 360/data.direction.shape[0]
    SPEC_naut_convection = np.zeros((data.freq.shape[0],data.direction.shape[0]))
    with open(project_name+str(dgm)+'_spec_'+days[0].strftime('%Y%m%d')+'_'+days[-1].strftime('%Y%m%d')+'.asc', 'w') as file_out:
                    file_out.write('SWAN   1\n')
                    file_out.write('$ Data produced by WAM3\n')
                    file_out.write('TIME\n')
                    file_out.write('          1\n')
                    file_out.write('LONLAT\n')    
                    file_out.write('          '+format(points.shape[0])+'\n')     
                    for k in range((points.shape[0])):
                        file_out.write('   '+format(points[k][1],'.4f')+'  '+format(points[k][0],'.4f')+'\n')
                    file_out.write('AFREQ\n')
                    file_out.write('          '+str(data.freq.shape[0])+'\n')
                    for l in range((data.freq.shape[0])):
                        file_out.write('   '+format(data.freq[l].values,'.4f')+'\n')
                    file_out.write('NDIR\n')
                    file_out.write('          '+format(data.direction.shape[0])+'\n')
                    for m in range((data.direction.shape[0])):
                        file_out.write('   '+format(data.direction[m].values,'.1f')+'\n') 
                    file_out.write('QUANT\n')
                    file_out.write('          1\n')
                    file_out.write('VaDens\n')
                    file_out.write('m2/Hz/degr \n')
                    file_out.write('-32767\n')
                    #first day
                    print('Generating 2d spectra at boundaries:')
                    print(days[0].strftime('%Y-%m-%d'))
                    for time_step in range(data.time.shape[0]):
                        file_out.write(str(data.time.time[time_step].values).split('-')[0]+str(data.time.time[time_step].values).split('-')[1]+\
                                       str(data.time.time[time_step].values).split('-')[2][:2]+'.'+str(data.time.time[time_step].values).split('-')[2][3:5]+'0000\n')
                        for p in range(points.shape[0]):
                            file_out.write('FACTOR\n')    
                            file_out.write(format(factor,'1.0E')+'\n')
                            # find the nearest grid point
                            diff_lat = np.abs(points[p][0]-data.latitude[0,:])
                            diff_lon = np.abs(points[p][1]-data.longitude[0,:])
                            add_lonlat = diff_lon + diff_lat
                            index_min_dinstance = np.where((add_lonlat >= np.sort(add_lonlat)[0]) & (add_lonlat <= np.sort(add_lonlat)[nr_spec_interpolate]))[0]
                            #print('Time,Point:'+str(time_step)+','+str(p))
                            #print(points[p][0],points[p][1])
                            #print(data.latitude[0,index_min_dinstance].values,data.longitude[0,index_min_dinstance].values)
                            #print(distance[index_min_dinstance],'km')
                            #print('Diff in km:',distance_2points(points[p][0],points[p][1],data.latitude[0,index_min_dinstance].values,data.longitude[0,index_min_dinstance].values))
                            #print('--------------------------')
                            #SPEC_ocean_convection = data.SPEC[time_step,0,index_min_dinstance,:,:].values
                            SPEC_ocean_convection = data.SPEC[time_step,0,index_min_dinstance,:,:].mean('x').values
                            SPEC_naut_convection[:,0:data.direction.shape[0]//2] = SPEC_ocean_convection[:,data.direction.shape[0]//2:] # Step 1a: 180..355 to start of array
                            SPEC_naut_convection[:,data.direction.shape[0]//2:]  = SPEC_ocean_convection[:,0:data.direction.shape[0]//2] # Step 1b: 0..175 to end of array
                            np.savetxt(file_out,SPEC_naut_convection/(delth*factor), fmt='%-10.0f') #         
                    #days excluding first and last days:
                    for i in range(1,len(days)-1):
                        print(days[i].strftime('%Y-%m-%d'))
                        url = 'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/'+days[i].strftime('%Y') +'/'+days[i].strftime('%m')+'/SPC'+days[i].strftime('%Y%m%d')+'00.nc'
                        data = xr.open_dataset(url)
                        for time_step in range(data.time.shape[0]):
                            file_out.write(str(data.time.time[time_step].values).split('-')[0]+str(data.time.time[time_step].values).split('-')[1]+\
                                           str(data.time.time[time_step].values).split('-')[2][:2]+'.'+str(data.time.time[time_step].values).split('-')[2][3:5]+'0000\n')
                            for p in range(points.shape[0]):
                                file_out.write('FACTOR\n')    
                                file_out.write(format(factor,'1.0E')+'\n')
                                # find the nearest grid point
                                diff_lat = np.abs(points[p][0]-data.latitude[0,:])
                                diff_lon = np.abs(points[p][1]-data.longitude[0,:])
                                add_lonlat = diff_lon + diff_lat
                                index_min_dinstance = np.where((add_lonlat >= np.sort(add_lonlat)[0]) & (add_lonlat <= np.sort(add_lonlat)[nr_spec_interpolate]))[0]
                                #print('Time,Point:'+str(time_step)+','+str(p))
                                #print(points[p][0],points[p][1])
                                #print(data.latitude[0,index_min_dinstance].values,data.longitude[0,index_min_dinstance].values)
                                #print(distance[index_min_dinstance],'km')
                                #print('Diff in km:',distance_2points(points[p][0],points[p][1],data.latitude[0,index_min_dinstance].values,data.longitude[0,index_min_dinstance].values))
                                #print('--------------------------')
                                #SPEC_ocean_convection = data.SPEC[time_step,0,index_min_dinstance,:,:].values
                                SPEC_ocean_convection = data.SPEC[time_step,0,index_min_dinstance,:,:].mean('x').values
                                SPEC_naut_convection[:,0:data.direction.shape[0]//2] = SPEC_ocean_convection[:,data.direction.shape[0]//2:] # Step 1a: 180..355 to start of array
                                SPEC_naut_convection[:,data.direction.shape[0]//2:]  = SPEC_ocean_convection[:,0:data.direction.shape[0]//2] # Step 1b: 0..175 to end of array
                                np.savetxt(file_out,calib_spec*SPEC_naut_convection/(delth*factor), fmt='%-10.0f') #     
                    #last day
                    if len(days)>1:
                        print(days[-1].strftime('%Y-%m-%d'))
                        url = 'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/'+days[-1].strftime('%Y') +'/'+days[-1].strftime('%m')+'/SPC'+days[-1].strftime('%Y%m%d')+'00.nc'
                        data = xr.open_dataset(url).sel(time=slice(days[-1].strftime('%Y-%m-%d')+ 'T00:00', end_date))
                        for time_step in range(data.time.shape[0]):
                            file_out.write(str(data.time.time[time_step].values).split('-')[0]+str(data.time.time[time_step].values).split('-')[1]+\
                                           str(data.time.time[time_step].values).split('-')[2][:2]+'.'+str(data.time.time[time_step].values).split('-')[2][3:5]+'0000\n')
                            for p in range(points.shape[0]):
                                file_out.write('FACTOR\n')    
                                file_out.write(format(factor,'1.0E')+'\n')
                                # find the nearest grid point
                                diff_lat = np.abs(points[p][0]-data.latitude[0,:])
                                diff_lon = np.abs(points[p][1]-data.longitude[0,:])
                                add_lonlat = diff_lon + diff_lat
                                index_min_dinstance = np.where((add_lonlat >= np.sort(add_lonlat)[0]) & (add_lonlat <= np.sort(add_lonlat)[nr_spec_interpolate]))[0]
                                #print('Time,Point:'+str(time_step)+','+str(p))
                                #print(points[p][0],points[p][1])
                                #print(data.latitude[0,index_min_dinstance].values,data.longitude[0,index_min_dinstance].values)
                                #print(distance[index_min_dinstance],'km')
                                #print('Diff in km:',distance_2points(points[p][0],points[p][1],data.latitude[0,index_min_dinstance].values,data.longitude[0,index_min_dinstance].values))
                                #print('--------------------------')
                                #SPEC_ocean_convection = data.SPEC[time_step,0,index_min_dinstance,:,:].values
                                SPEC_ocean_convection = data.SPEC[time_step,0,index_min_dinstance,:,:].mean('x').values
                                SPEC_naut_convection[:,0:data.direction.shape[0]//2] = SPEC_ocean_convection[:,data.direction.shape[0]//2:] # Step 1a: 180..355 to start of array
                                SPEC_naut_convection[:,data.direction.shape[0]//2:]  = SPEC_ocean_convection[:,0:data.direction.shape[0]//2] # Step 1b: 0..175 to end of array
                                np.savetxt(file_out,SPEC_naut_convection/(delth*factor), fmt='%-10.0f') # 



def u_v_from_dir(ws,wdir):
# see http://tornado.sfsu.edu/geosciences/classes/m430/Wind/WindDirection.html
    u = -ws * (np.sin(np.deg2rad(wdir)))
    v = -ws * (np.cos(np.deg2rad(wdir)))
    return u,v


def generate_wind_forcing_SWAN(project_name,dgm, start_date , end_date, lon_min, lat_min, lon_max, lat_max, NX , NY):
    days = pd.date_range(start=start_date.split('T')[0], end=end_date.split('T')[0], freq='D')
    d_lon = ((lon_max-lon_min)/NX)
    d_lat = ((lat_max-lat_min)/NY)
    nc_fimex ='temp.nc'
    print('Generation of wind forcing:')
    with open(project_name + str(dgm)+'_wind_'+days[0].strftime('%Y%m%d')+'_'+days[-1].strftime('%Y%m%d')+'.asc', 'w') as file_out:
        for i in range(len(days)):
            print(days[i].strftime('%Y-%m-%d'))
            url = 'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/'+days[i].strftime('%Y')+'/'+days[i].strftime('%m')+'/'+days[i].strftime('%Y%m%d')+'_MyWam3km_hindcast.nc' 
            if i==0:
                start_date_fimex = start_date
                if len(days)>1:
                    end_date_fimex = days[i].strftime('%Y-%m-%dT23:00')
                else:
                    end_date_fimex = end_date
            elif i==len(days)-1:
                if len(days)>1:
                    start_date_fimex = days[i].strftime('%Y-%m-%dT00:00')
                    end_date_fimex   = end_date
            else:    
                start_date_fimex = days[i].strftime('%Y-%m-%dT00:00')
                end_date_fimex   = days[i].strftime('%Y-%m-%dT23:00')
    
            call(['fimex-1.6', '--input.file='+url,
                         '--interpolate.method=bilinear',
                         '--interpolate.projString=+proj=latlong +ellps=sphere +a=6371000 +e=0',
                         '--interpolate.xAxisValues='+str(lon_min)+','+str(lon_min+d_lon)+',...,'+str(lon_max)+'',
                         '--interpolate.yAxisValues='+str(lat_min)+','+str(lat_min+d_lat)+',...,'+str(lat_max)+'',
                         '--interpolate.xAxisUnit=degree', '--interpolate.yAxisUnit=degree',
                         '--process.rotateVector.all',
                         '--extract.reduceTime.start='+start_date_fimex,'--extract.reduceTime.end='+end_date_fimex,
                         '--extract.selectVariables=ff','--extract.selectVariables=dd',
                         '--extract.reduceDimension.start=1',
                         '--extract.reduceDimension.end=1',
                         '--process.rotateVector.direction=latlon',
                         '--output.file='+nc_fimex])
    
            data = xr.open_dataset(nc_fimex)
            u,v = u_v_from_dir(data.ff, data.dd) # factor 1000
            u = u.fillna(0)
            v = v.fillna(0)
    
            for time_step in range(data.time.shape[0]):
                file_out.write(str(data.time.time[time_step].values).split('-')[0]+str(data.time.time[time_step].values).split('-')[1]+\
                               str(data.time.time[time_step].values).split('-')[2][:2]+'.'+str(data.time.time[time_step].values).split('-')[2][3:5]+'0000\n')
                np.savetxt(file_out,u[time_step]*1000, fmt='%i') #
                file_out.write(str(data.time.time[time_step].values).split('-')[0]+str(data.time.time[time_step].values).split('-')[1]+\
                               str(data.time.time[time_step].values).split('-')[2][:2]+'.'+str(data.time.time[time_step].values).split('-')[2][3:5]+'0000\n')
                np.savetxt(file_out,v[time_step]*1000, fmt='%i') #

def run_SWAN(input_file,swan_directory):
    print('Running SWAN----------------------->>>>>>>>>>>>>>>>>>>>>>>>>>')
    p = Popen(['./swanrun','-input', input_file], cwd=swan_directory)
    p.wait()                                               
    

