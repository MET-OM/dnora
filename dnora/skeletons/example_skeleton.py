from spc_mod import Spectra, Boundary
from grd_mod import Grid, UnstrGrid
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

grid = Grid(lon=(22.00, 22.73), lat=(60.00, 60.53), name='Skjerjehamn')
#ugrid = UnstrGrid(lon=(22.00, 22.73), lat=(60.00, 60.53), name='Skjerjehamn')
#ugrid.import_topo(topo_reader=grd.read.EMODNET2020(tile='*', folder='~/Documents/EMODNET2020'))
grid.set_spacing(nx=5, ny=5)
#grid.set_spacing(dlon=0.1, dlat=0.1)
#bnd_set = grd.boundary.EdgesAsBoundary(edges=['N', 'W', 'S'])
#grid.set_boundary(boundary_setter=bnd_set)
#grid.import_topo(topo_reader=grd.read.EMODNET2020(tile='*', folder='~/Documents/EMODNET2020'))
#grid.import_topo(grid)
#grid.mesh_grid()
#
#spec = Spectra(lon=(4.00, 5.73), lat=(60.53, 61.25), name='Skjerjehamn')
#bnd = Boundary(lon=(4.00), lat=(60.53), name='Skjerjehamn')

#spec = GriddedSpectra(lon=(4.00, 5.73), lat=(60.53, 61.25), name='Skjerjehamn')
#
# gspec = GriddedSpectra(lon=(4.00, 5.73), lat=(60.53, 61.25), name='Skjerjehamn')
