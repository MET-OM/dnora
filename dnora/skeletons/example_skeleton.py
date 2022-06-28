from spc_mod import Spectra
from grd_mod import Grid
import numpy as np
import xarray as xr
from dnora import grd
# psktn = PointSkeleton(lon=np.arange(10,20), lat=np.arange(60,70))
#
# ds=psktn.compile_to_xr(np.arange(80,90),'test_data')
#
#
#
# gsktn = GriddedSkeleton(lon=np.arange(10,20), lat=np.arange(60,63), time=[0,1])
# ds2=gsktn.compile_to_xr(np.zeros((3,10,2)),'test_data')

# spec = Spectra(lon=np.arange(10,20), lat=np.arange(60,70), time=[0,1])
# ds2=spec.compile_to_xr(np.zeros((10,2)),'test_data')
#
# ds2=spec.compile_to_xr(np.zeros((10,2,6)),'test_spec', additional_coords={'freq': np.array([0,1,2,3,4,5])})

grid = Grid(lon=(4.00, 5.73), lat=(60.53, 61.25), name='Skjerjehamn')
grid.set_spacing(dm=10000)
grid.import_topo(topo_reader=grd.read.EMODNET2020(tile='*'))
grid.mesh_grid()
