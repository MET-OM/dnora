#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 10:27:59 2022

@author: emiliebyermoen
"""
import xarray as xr
from matplotlib import pyplot as plt
import numpy as np
import cmocean.cm
from netCDF4 import Dataset
from netCDF4 import num2date
from cartopy import crs as ccrs
from cartopy import feature as cfeature
import pandas as pd

#%%

def plotArrows(ax, lon, lat, xdata, ydata, var, scale=100, reduce_arrows: int=None):
    """

    Parameters
    ----------
    ax : ax
        DESCRIPTION.
    lon : np.array or alike. len=m
        DESCRIPTION.
    lat : np.array or alike. len=n
        DESCRIPTION.
    xdata : gridded data [mxn]
        DESCRIPTION.
    ydata : gridded data [mxn]
        DESCRIPTION.
    reduce_arrows : int, optional
        The default is 2. Higher numer=reduce number of arrows more
    scale : int, optional
        Scale length of the arrow, the width and the length of arrow head. The default is 100.

    Returns
    -------
    ax : ax, with arrows

    """
    if reduce_arrows is not None:
        step_lat = reduce_arrows
        step_lon = reduce_arrows
    else:
        step_lat = round(len(lat)/10)
        step_lon = round(len(lon)/10)

    for m in range(0,len(lon), step_lon):
        for n in range(0,len(lat), step_lat):
            ax.arrow(lon[m], lat[n], xdata[n][m]/scale, ydata[n][m]/scale, color='white',
                     linewidth=0.15, head_width=2/scale, head_length=2/scale, overhang=1) #linewidth=.02, head_width=.01, head_length=.01
    return ax

﻿def plot_magnitude(fig, ax, lon, lat, data, var, vmin=None, vmax=None, cbar=False):
    """
    Takes lon and lat positions, and data points. Plots countourplot on given ax.
    """
    if vmin is None:
        vmin = np.min(data)
    if vmax is None:
        vmax = np.max(data)
    unit = default[var]['unit']
    cmap = default[var]['cmap']

    xx, yy = np.meshgrid(lon, lat)
    if not cbar:
        levels = np.linspace(vmin,vmax,8)
    else:
        levels = cbar.boundaries

    cont = ax.contourf(xx, yy, data, cmap=cmap, levels=levels)

    if not cbar:
        if var != 'mask':
            # creating axes for the colorbar
            if ax.get_position().height >= ax.get_position().width:
                orientation = 'vertical'
                cax = fig.add_axes([ax.get_position().x1+0.04,ax.get_position().y0,0.02,ax.get_position().height])
            elif ax.get_position().height < ax.get_position().width:
                orientation= 'horizontal'
                cax = fig.add_axes([ax.get_position().x0, ax.get_position().y0-0.1,ax.get_position().width,0.03])
            cbar = fig.colorbar(cont, orientation=orientation, cax=cax, label=f"[{unit}]")

    return fig, ax, cbar

def mask_land_from_topo(ax, lon, lat, topo):
    topo_lon = topo[:,0]
    loto_lat = topo[0,:]
    lon_min_inx = np.where(lon>min(lon))[0][0]
    lon_max_inx = np.where(lon<max(lon))[0][-1]
    lat_min_inx = np.where(lat>min(lat))[0][0]
    lat_max_inx = np.where(lat<max(lat))[0][-1]
    lon = lon[lon_min_inx:lon_max_inx]
    lat = lat[lat_min_inx:lat_max_inx]
    xdata = xdata[lat_min_inx:lat_max_inx,lon_min_inx:lon_max_inx]

    topo=pd.DataFrame(topo)
    topo[topo>=0] = np.nan #ocean
    topo[topo<0] = 1 #land
    xx, yy = np.meshgrid(lon, lat)
    ax.contourf(xx,yy,topo, cmap='gray')
    return ax


def plot_field(lon, lat, xdata, ydata=None, var: str='ff', titleStr=' ', boundary=None, scale=100,
                vmin=None, vmax=None, cbar=None,
               obs_lon: float=None, obs_lat: float=None, obs_value: float=None):
    """
    Parameters
    ----------
    lon : np.array or alike
        logitude
    lat : np.array or alike
        latitude
    xdata : gridded or listed data
        can be total wind, wind in x-direction, significant wave height, topography.
    ydata : optional additional data
        can be wind in y-direction
    var : str, optional
        Modifies the plot according to the type of variable.
        wind: 'ff'
        significant wave height: 'hs'
        topography/bathymetry: 'topo'
    title : str. to add to the title of the plot
    boundary : list of [lon_min, lon_max, lat_min, lat_max]
    scale : used to scale arrow length if plotting wind. default is 100

    Returns
    -------
    None. Plots a field with the given data.

    """

    #INSERT list of data --> gridded data (only for data that comes in list)

    # Slicing lat, lon and data with the boundary conditions if given
    if boundary is not None:
        lon_min = boundary[0]
        lon_max = boundary[1]
        lat_min = boundary[2]
        lat_max = boundary[3]
        lon_min_inx = np.where(lon>lon_min)[0][0]
        lon_max_inx = np.where(lon<lon_max)[0][-1]
        lat_min_inx = np.where(lat>lat_min)[0][0]
        lat_max_inx = np.where(lat<lat_max)[0][-1]
        lon = lon[lon_min_inx:lon_max_inx]
        lat = lat[lat_min_inx:lat_max_inx]
        xdata = xdata[lat_min_inx:lat_max_inx,lon_min_inx:lon_max_inx]
        if ydata is not None:
            ydata = ydata[lat_min_inx:lat_max_inx,lon_min_inx:lon_max_inx]

    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}) #AlbersEqualArea didnt work
    ax.set(title='{} {}'.format(default[var]['name'], titleStr))
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False

    if ydata is not None and var == 'ff': #interpretes the data as wind in x- and y-direction
        # plots magnitude of wind
        windMagnitude = (xdata**2 + ydata**2)**0.5
        ax = plotMagnitude(ax, fig, lon, lat, windMagnitude, var, vmin=vmin, vmax=vmax, cbar=cbar)
        # Plot arrows
        ax = plotArrows(ax, lon, lat, xdata, ydata, var, scale)
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m', facecolor='none', edgecolor='black'))
        #ax.plot()
    elif ydata is not None and var == 'vel': #interpretes the data as wind in x- and y-direction
        # plots magnitude of vel
        velMagnitude = (xdata**2 + ydata**2)**0.5
        velMagnitude = np.where(velMagnitude!=0, velMagnitude, np.nan) # only wet points to plots
        ax = plotMagnitude(ax, fig, lon, lat, velMagnitude, var, vmin=vmin, vmax=vmax, cbar=cbar)
        # Plot arrows
        ax = plotArrows(ax, lon, lat, xdata, ydata, var, scale)
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m', facecolor='none', edgecolor='black'))
        #ax.plot()
    else:
        # Plot only magnutide of field
        ax = plotMagnitude(ax, fig, lon, lat, xdata, var)
        if var =='hs':
            #xdata = pd.DataFrame(xdata)
            #xdata[xdata>=0]=np.nan
            #xdata[xdata<0]=1
            #plotMagnitude(ax, fig,lon,lat,xdata,var='mask')
            ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m', facecolor=[.8,.8,.8], edgecolor='black'))
        #ax.plot()
    if obs_value is not None:
        ax=observation(ax, obs_lon, obs_lat, obs_value, var, levels) #need to get the levels from plotMagnutide function
    return {'fig':fig,'ax':ax, 'cbar':cbar}


#%% Add observation plot

def observation(ax, lon, lat, xdata, var, levels, ydata=None):
    if ydata is not None:
        magnitude = (xdata**2+ydata**2)**0.5
    else:
        magnitude = xdata
    ax.scatter(lon, lat, s=50, c=magnitude, cmap=default[var]['cmap'],levels=levels, edgecolor='red') # cmap=default[var]['cmap'])
    #ax.arrow(lon, lat, xdata/scale, ydata/scale, color='black',
             #linewidth=0.15, head_width=2/scale, head_length=2/scale, overhang=1)
    return ax


#%%

default = {'hs': {'name': 'Significant wave height', 'unit':'m', 'cmap': cmocean.cm.amp},
           'ff': {'name': 'Wind', 'unit':'m/s', 'cmap': cmocean.cm.tempo},
           'topo':{'name': 'Topography', 'unit':'m', 'cmap': cmocean.cm.topo_r},
           'mask':{'name': ' ', 'unit':' ', 'cmap': 'gray'}
           }


#%% Example plot, dataset with wind
ds = xr.open_dataset('/Users/emiliebyermoen/Dropbox/Energi/ENERGI240/wind_0003_MetNo_MEPS.nc')

plot_field(ds['x'], ds['y'],ds['x_wind_10m'][0][0],ds['y_wind_10m'][0][0],var='ff', scale=300, boundary=[5,5.5,61,61.2])




#%% Example plot, dataset with wave height
data = Dataset('https://thredds.met.no/thredds/dodsC/e39_models/SWAN250/Sula/swanSula202001.nc','r')
lat = data.variables['latitude'][:] # degrees north, 222 points
lon = data.variables['longitude'][:] # degrees east, 280 points
time = data.variables['time'][:] # 744
actual_time = num2date(time, 'seconds since 1970-01-01', 'gregorian')
hs = data.variables['hs'][0][:][:] # sea surface wave significant height [m] (time, latitude, longitude)


plot_field(lon, lat ,hs, var='hs', titleStr=str(actual_time[0]))

#%% Exaple plot, dataset wind 6. april
data2 = Dataset('https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc','r')

time = data2.variables['time'][:] #62
actual_time = num2date(time, 'seconds since 1970-01-01', 'gregorian')
lat = data2.variables['latitude'][:,0] #shape (1069, 949)
lon = data2.variables['longitude'][0,:] #shape (1069, 949)
wind_x10m = data2.variables['x_wind_10m']
wind_y10m = data2.variables['y_wind_10m']

time_step = 30
plot_field(lon, lat, wind_x10m[time_step,0,0,:,:], wind_y10m[time_step,0,0,:,:], var='ff',scale=50, titleStr=str(actual_time[time_step]),boundary=[3,7,58,62])

#%% Example plot, observation from Sulafjorden

obs = Dataset('https://thredds.met.no/thredds/dodsC/obs/buoy-svv-e39/2022/04/202204_E39_C_Sulafjorden_wind.nc','r')
obs_lat = obs['latitude'][0]
obs_lon = obs['longitude'][0]
obs_wind_speed = obs['WindSpeed'][0]

ax = plot_field(lon, lat, wind_x10m[time_step,0,0,:,:], wind_y10m[time_step,0,0,:,:], var='ff',scale=300,
                titleStr=str(actual_time[time_step]),boundary=[5.8,6.3,62.2,62.5], obs_lon=obs_lon, obs_lat=obs_lat, obs_value=obs_wind_speed)
#observation(ax, obs_lon, obs_lat, 'ff',obs_wind_speed)


#%% Plotting topograpgy example

topo = Dataset('/Users/emiliebyermoen/Dropbox/Energi/ENERGI240/dnora-lokal/data/GEBCO_09_Mar_2022_d1ac3aae417d/gebco_2021_n61.0_s59.0_w4.0_e6.0.nc')

lon = topo['lon'][:]
lat = topo['lat'][:]
elev = topo['elevation'][:]*(-1)

plot_field(lon, lat, elev, var='topo')

#%% Plotting coastline from topography dataset example

elev=pd.DataFrame(elev)

elev[elev>=0] = np.nan
elev[elev<0] = 1


plot_field(lon,lat,elev,var='mask')

#%%¨
fig, ax = plt.subplots()

ax=mask_land_from_topo(ax,lon,lat, elev)

ax.plot()

#%%

#reduce_arrows=10
time_step=0
"""
if reduce_arrows is not None:
    step_lat = reduce_arrows
    step_lon = reduce_arrows
else:"""
step_lat = round(len(lat)/10)
step_lon = round(len(lon)/10)

mask_lon = range(0,len(lon), step_lon)
mask_lat = range(0,len(lat), step_lat)

xx, yy = np.meshgrid(lon[mask_lon], lat[mask_lat])
u = wind_x10m[time_step,0,0,mask_lat, mask_lon]
v = wind_y10m[time_step,0,0,mask_lat, mask_lon]
plt.quiver(xx,yy,u,v)


"""
mask_lon = []
for x in xx[0,:]:
    mask_lon.append(np.where(lon==x)[0][0])

mask_lat = []
for y in yy[0,:]:
    mask_lat.append(np.where(lat==y)[0][0])

u = wind_x10m[time_step,0,0,mask_lat, mask_lon]
v = wind_y10m[time_step,0,0,mask_lat, mask_lon]
plt.quiver(xx,yy,u,v)
"""
