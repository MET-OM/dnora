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
from .defaults import default
#%%

def plot_barbs(fig, ax, lon, lat, xdata, ydata, var, reduce_arrows: int=None):
    # from m/s to knots
    xdata = xdata*1.9438
    ydata = ydata*1.9438

    """
    if reduce_arrows is not None:
        step_lat = reduce_arrows
        step_lon = reduce_arrows
    else:
        step_lat = np.floor(len(lat)/10)
        step_lon = np.floor(len(lon)/10)

    lon_red = []
    for n in range(0,len(lon), step_lon):
        lon_red.append(lon[n])

    lat_red = []
    for m in range(0,len(lat), step_lat):
        lat_red.append(lat[m])

    x = []
    for n in range(0,len(lon), step_lon):
        list = []
        for m in range(0, len(lat), step_lat):
            list.append(xdata[n][m])
        x.append(list)

    y = []
    for n in range(0,len(lat), step_lat):
        list = []
        for m in range(0, len(lon), step_lon):
            list.append(ydata[n][m])
        y.append(list)
        """
    ax.barbs(lon, lat, xdata, ydata,
             sizes=dict(emptybarb=0.2, spacing=0.2, height=0.5),
             linewidth=0.5, transform= ccrs.PlateCarree(),
             barbcolor='white', flagcolor='white', length=4)
    #ax.barbs(lon_red,lat_red,x,y, barbcolor='white', flagcolor='white', length=3)


    return fig, ax



def plot_arrows(fig, ax, lon, lat, xdata, ydata, var, scale=100, reduce_arrows: int=None):
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
    return fig, ax

def plot_magnitude(fig, ax, lon, lat, data, var, vmin=None, vmax=None, cbar=True):
    """
    Takes lon and lat positions, and data points. Plots countourplot on given ax.
    """
    unit = default[var]['unit']
    cmap = default[var]['cmap']
    if vmin is None:
        vmin = np.min(data)
    if vmax is None:
        vmax = np.max(data)
    vmin = np.floor(vmin).astype(int)
    vmax = np.ceil(vmax).astype(int)
    if vmax == vmin:
        vmax += 1
    levels = np.linspace(vmin, vmax, vmax-vmin+1)
    xx, yy = np.meshgrid(lon, lat)
    cont = ax.contourf(xx, yy, data, cmap=cmap, levels=levels)
    if cbar:
        if var != 'mask':
            # creating axes for the colorbar
            if ax.get_position().height >= ax.get_position().width:
                orientation = 'vertical'
                cax = fig.add_axes([ax.get_position().x1+0.04,ax.get_position().y0,0.02,ax.get_position().height])
            elif ax.get_position().height < ax.get_position().width:
                orientation= 'horizontal'
                cax = fig.add_axes([ax.get_position().x0, ax.get_position().y0-0.1,ax.get_position().width,0.03])
            cbar = fig.colorbar(cont, orientation=orientation, cax=cax, label=f"[{unit}]")
    else:
        cbar = None
    return fig, ax, cbar

def mask_land_from_topo(fig, ax, lon, lat, topo):
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
    return fig, ax


def plot_field(lon, lat, xdata, ydata=None, fig=None, ax=None, position=111,
               var: str='ff', title_str=' ', boundary=None, scale=100,
               vmin=None, vmax=None, cbar=True, barbs=False,
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
    # Slicing lat, lon and data with the boundary conditions if given
    if boundary is not None:
        lon_min, lon_max, lat_min, lat_max = boundary
        lon_min_inx = np.where(lon>lon_min)[0][0]
        lon_max_inx = np.where(lon<lon_max)[0][-1]
        lat_min_inx = np.where(lat>lat_min)[0][0]
        lat_max_inx = np.where(lat<lat_max)[0][-1]
        lon = lon[lon_min_inx:lon_max_inx]
        lat = lat[lat_min_inx:lat_max_inx]
        xdata = xdata[lat_min_inx:lat_max_inx,lon_min_inx:lon_max_inx]
        if ydata is not None:
            ydata = ydata[lat_min_inx:lat_max_inx,lon_min_inx:lon_max_inx]

    if fig is None:
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    else:
        fig.delaxes(ax)
        ax = fig.add_subplot(position, projection=ccrs.PlateCarree())
    ax.set(title=f"{default[var]['name']} {title_str}")
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    if ydata is not None and var == 'ff': #interpretes the data as wind in x- and y-direction
        windMagnitude = (xdata**2 + ydata**2)**0.5
        fig, ax, cbar = plot_magnitude(fig, ax, lon, lat, windMagnitude, var, vmin=vmin, vmax=vmax, cbar=cbar)
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m', facecolor='none', edgecolor='black'))
        if barbs:
            fig, ax = plot_barbs(fig, ax, lon, lat, xdata, ydata, var, scale)
        else:
            fig, ax = plot_arrows(fig, ax, lon, lat, xdata, ydata, var, scale)

    else:
        fig, ax, cbar = plot_magnitude(fig, ax, lon, lat, xdata, var, vmin=vmin, vmax=vmax, cbar=cbar)
        if var =='hs':
            #xdata = pd.DataFrame(xdata)
            #xdata[xdata>=0]=np.nan
            #xdata[xdata<0]=1
            #plotMagnitude(ax, fig,lon,lat,xdata,var='mask')
            ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m', facecolor=[.8,.8,.8], edgecolor='black'))
    if obs_value is not None:
       fig, ax =observation(fig, ax, obs_lon, obs_lat, obs_value, var, cbar) #need to get the levels from plotMagnutide function



    return {'fig':fig,'ax':ax, 'cbar':cbar}


#%% Add observation plot

def observation(fig, ax, lon, lat, xdata, var, levels, ydata=None):
    if ydata is not None:
        magnitude = (xdata**2+ydata**2)**0.5
    else:
        magnitude = xdata
    ax.scatter(lon, lat, s=50, c=magnitude, cmap=default[var]['cmap'], edgecolor='red') # cmap=default[var]['cmap'])
    #ax.arrow(lon, lat, xdata/scale, ydata/scale, color='black',
             #linewidth=0.15, head_width=2/scale, head_length=2/scale, overhang=1)
    return fig, ax
