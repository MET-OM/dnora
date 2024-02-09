"""
Created on Mon Nov  8 08:58:23 2021

@author: KonstantinChri
"""

import numpy as np
import xarray as xr
import scipy.io


def convert_HOS_3Ddat_to_netcdf(input_file,output_file):
    """
    Function to convert the HOS-ocean file e.g., 3d.dat
    to a netcdf file
    Parameters:
    input_file = HOS-ocean output of eta, e.g., 3d.dat
    output_file = filename.nc
    ----------
    eta  : 3D free surface elevation [time,y, x]
   Returns
   -------
    ds : xr.Dataset
        eta: 3D surface elevation [time,y, x]
    """
    with open(input_file) as f:
        lines = f.readlines()

    # remove lines with # comments
    lines = [x for x in lines if not x.startswith('#')]
    lines = [s.replace("\n", "") for s in lines]

    I = int(lines[2].split(',')[1].split('=')[1])
    J = int(lines[2].split()[-1])

    title = lines[0].split('"')[1]

    #remove lines with titles
    lines = [x for x in lines if not x.startswith('T')]
    lines = [x for x in lines if not x.startswith('V')]
    lines = [x for x in lines if not x.startswith('Z')]
    lines = [s.split() for s in lines]

    timestep = int(len(lines)/(I*J))

    x = np.zeros(I*J)
    y = np.zeros(I*J)
    eta = np.zeros(len(lines))

    for i in range(len(lines)):
        if i < I*J:
            x[i] = float(lines[i][0])
            y[i] = float(lines[i][1])
            eta[i] = float(lines[i][2])
        else:
            eta[i] = float(lines[i][0])

    eta_3d = np.zeros((timestep,J,I))
    eta_3d[0,:,:] = eta[0:I*J].reshape((-J, I))

    # fill x, y
    x = x[0:I]
    y = y[0::I]

    # fill eta_3d
    for t in range(timestep):
        eta_3d[t,:,:] = eta[t*(I*J):(t+1)*(I*J)].reshape((-J, I))

    # create xarray
    ds = xr.Dataset({'eta': xr.DataArray(eta_3d,
                            coords={'time': np.arange(timestep),'y': y, 'x': x},
                            dims=["time", "y", "x"],
                            attrs  = {'units': 'm','long_name':title})})
    #save xarray ro netcdf
    ds.to_netcdf(output_file)

    return ds

