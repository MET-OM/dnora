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

def plot_spectra(freq, spec, mdir, spr, title_str, fig_dict, ymax):
    if fig_dict is None:
        fig_dict = {}

    fig = fig_dict.get('fig')
    ax = fig_dict.get('ax')
    ax2 = fig_dict.get('ax2')
    if fig is None:
        fig, ax = plt.subplots()
        fig_dict['fig']=fig
        fig_dict['ax']=ax
    else:
        fig.delaxes(ax)
        if ax2 is not None:
            fig.delaxes(ax2)
        ax = fig.add_subplot()
        fig_dict['ax']=ax

    if ymax is None:
        ymax = np.max(spec)

    ax.plot(freq, spec, color='black')
    ax.set_ylim([0,ymax])
    ax.set_title(title_str)
    ax.set_xlabel('f (Hz)')
    ax.set_ylabel('E(f) (m**2/Hz)')

    ax2 = ax.twinx()
    fig_dict['ax2']=ax2
    ax2.plot(freq, np.mod(mdir+180, 360), color='gray')
    ax2.plot(freq, np.mod(mdir+180+spr, 360), color='gray', linestyle='--')
    ax2.plot(freq, np.mod(mdir+180-spr, 360), color='gray', linestyle='--')
    ax2.set_ylim([0, 360])
    ax2.set_ylabel('Mean direction (deg)')

    return fig_dict

def plot_polar_spectra(freq, dirs, spec, title_str, fig_dict, vmax, vmin, cbar=True):

    if fig_dict is None:
        fig_dict = {}

    fig = fig_dict.get('fig')
    ax = fig_dict.get('ax')
    if fig is None:
        fig, ax = plt.subplots(polar=True)
        fig_dict['fig']=fig
        fig_dict['ax']=ax
    else:
        fig.delaxes(ax)
        ax = fig.add_subplot(polar=True)
        fig_dict['ax']=ax

    if vmax is None:
        vmax = np.max(spec)

    if vmin is None:
        vmin = np.min(spec)
    #ds.isel(site=0,time=i).efth.spec.split(fmin=0.04).spec.plot
    last_row = np.transpose([spec[:,0]])
    spec_plot = np.hstack([spec, last_row])
    dir_plot = np.hstack([dirs, dirs[0]+360])

    # if vmax-vmin<20:
    #     levels = np.linspace(vmin, vmax, np.floor(vmax-vmin+1).astype(int))
    # # else:
    levels = np.linspace(vmin, vmax, 20)
    print(vmin)
    cont = ax.contourf(np.deg2rad(dir_plot), freq, spec_plot, cmap="ocean_r",vmin=vmin, vmax=vmax, levels=levels)
    fig_dict['cont'] = cont
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_title(title_str)

    orientation = 'vertical'

    if cbar:
        cax = fig.add_axes([ax.get_position().x1+0.04,ax.get_position().y0,0.02,ax.get_position().height])
        cbar = fig.colorbar(cont, orientation=orientation, cax=cax, label=f"E(f, theta) (m**2/Hz/rad)")
        fig_dict['cbar'] = cbar
        fig_dict['cax'] = cax
    return fig_dict


def plot_barbs(fig_dict, lon, lat, xdata, ydata, var, reduce_arrows: int=None):
    # from m/s to knots
    xdata = xdata*1.9438
    ydata = ydata*1.9438

    ax = fig_dict.get('ax')
    ax.barbs(lon, lat, xdata, ydata,
             sizes=dict(emptybarb=0.2, spacing=0.2, height=0.5),
             linewidth=0.5, transform= ccrs.PlateCarree(),
             barbcolor='white', flagcolor='white', length=4)
    #ax.barbs(lon_red,lat_red,x,y, barbcolor='white', flagcolor='white', length=3)

    return fig_dict



def plot_arrows(fig_dict, lon, lat, xdata, ydata, var, scale=100, reduce_arrows: int=None):
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
        step_lat = max(round(len(lat)/10),1)
        step_lon = max(round(len(lon)/10),1)

    ax = fig_dict.get('ax')
    for m in range(0,len(lon), step_lon):
        for n in range(0,len(lat), step_lat):
            ax.arrow(lon[m], lat[n], xdata[n][m]/scale, ydata[n][m]/scale, color='white',
                     linewidth=0.15, head_width=2/scale, head_length=2/scale, overhang=1) #linewidth=.02, head_width=.01, head_length=.01
    return fig_dict


def plot_magnitude(fig_dict, lon, lat, data, var, vmin=None, vmax=None, cbar=True):
    """
    Takes lon and lat positions, and data points. Plots countourplot on given ax.
    """
    fig = fig_dict.get('fig')
    ax = fig_dict.get('ax')
    unit = default[var]['unit']
    cmap = default[var]['cmap']

    if vmin is None:
        vmin = np.min(data)
    if vmax is None:
        vmax = np.max(data)
    vmin = np.floor(vmin).astype(int)
    vmax = np.ceil(vmax).astype(int)

    if vmax-vmin<20:
        levels = np.linspace(vmin, vmax, np.floor(vmax-vmin+1).astype(int))
    else:
        levels = np.linspace(vmin, vmax, 11)

    xx, yy = np.meshgrid(lon, lat)
    if len(levels) > 1:
        cont = ax.contourf(xx, yy, data, cmap=cmap, levels=levels)
    else:
        cont = ax.pcolor(lon, lat, data, cmap=cmap)
    fig_dict['levels']=levels

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
            fig_dict['cbar'] = cbar
    else:
        cbar = None
    return fig_dict

def mask_land_from_topo(fig_dict, lon, lat, topo):
    ax = fig_dict.get('ax')

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
    return fig_dict


def plot_field(lon, lat, xdata, ydata=None, fig_dict=None, position=111,
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
    fig_dict : dictionary with figure properties
        can be 'ax', 'fig', 'cbar', 'levels' etc.
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
    Dictionary with given parameters for plot, such as 'fig', 'ax' and other optional such as 'cbar' and 'levels'

    """
    if fig_dict is None:
        fig_dict = {}

    fig = fig_dict.get('fig')
    ax = fig_dict.get('ax')


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
        fig_dict['fig']=fig
        fig_dict['ax']=ax
    else:
        fig.delaxes(ax)
        ax = fig.add_subplot(position, projection=ccrs.PlateCarree())
        fig_dict['ax']=ax
    ax.set(title=f"{default[var]['name']} {title_str}")
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = None
    gl.right_labels = None

    if ydata is not None and var == 'ff': #interpretes the data as wind in x- and y-direction
        windMagnitude = (xdata**2 + ydata**2)**0.5

        fig_dict = plot_magnitude(fig_dict, lon, lat, windMagnitude, var, vmin=vmin, vmax=vmax, cbar=cbar)
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m', facecolor='none', edgecolor='black'))
        if barbs:
            fig_dict = plot_barbs(fig_dict, lon, lat, xdata, ydata, var, scale)
        else:
            fig_dict = plot_arrows(fig_dict, lon, lat, xdata, ydata, var, scale)
    elif ydata is not None and var == 'vel': #interpretes the data as wind in x- and y-direction
         velMagnitude = (xdata**2 + ydata**2)**0.5
         velMagnitude = np.where(velMagnitude!=0, velMagnitude, np.nan) # only wet points to plots

         fig_dict = plot_magnitude(fig_dict, lon, lat, velMagnitude, var, vmin=vmin, vmax=vmax, cbar=cbar)
         ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m', facecolor='none', edgecolor='black'))
         if barbs:
             fig_dict = plot_barbs(fig_dict, lon, lat, xdata, ydata, var, scale)
         else:
             fig_dict = plot_arrows(fig_dict, lon, lat, xdata, ydata, var, scale)

    else:
        fig_dict = plot_magnitude(fig_dict, lon, lat, xdata, var, vmin=vmin, vmax=vmax, cbar=cbar)
        if var =='hs':
            #xdata = pd.DataFrame(xdata)
            #xdata[xdata>=0]=np.nan
            #xdata[xdata<0]=1
            #plotMagnitude(ax, fig,lon,lat,xdata,var='mask')
            ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m', facecolor=[.8,.8,.8], edgecolor='black'))
    if obs_value is not None:
       fig_dict =observation(fig_dict, obs_lon, obs_lat, obs_value, var, cbar) #need to get the levels from plotMagnutide function

    ax.set_xlim([min(lon),max(lon)])
    ax.set_ylim([min(lat),max(lat)])

    return fig_dict


#%% Add observation plot

def observation(fig_dict, lon, lat, xdata, var, levels, ydata=None):
    ax = fig_dict.get('ax')
    if ydata is not None:
        magnitude = (xdata**2+ydata**2)**0.5
    else:
        magnitude = xdata
    ax.scatter(lon, lat, s=50, c=magnitude, cmap=default[var]['cmap'], edgecolor='red') # cmap=default[var]['cmap'])
    #ax.arrow(lon, lat, xdata/scale, ydata/scale, color='black',
             #linewidth=0.15, head_width=2/scale, head_length=2/scale, overhang=1)
    return fig_dict
