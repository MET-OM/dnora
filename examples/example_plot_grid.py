# =============================================================================
# IMPORT dnora
# =============================================================================
import sys
dnora_directory = '../'
sys.path.insert(0, dnora_directory)
from dnora import bnd, grd, mdl, wnd
from dnora.dnplot import grd_topo
from dnora import dnplot

lon_min = 4.0
lat_min = 60.53
lon_max = 5.73
lat_max = 61.25

grid = grd.Grid(lon_min, lon_max, lat_min, lat_max, name='Skjerjehamn')

grid.set_spacing(dm=1000)
bnd_set = grd.boundary.EdgesAsBoundary(edges = ['N', 'W', 'S'], step = 20)
grid.set_boundary(boundary_setter = bnd_set)
# Import topography and mesh it down to the grid definitions
grid.import_topo(topo_reader=grd.read.EMODNET2018(tile='D5'))
grid.mesh_grid()


model = mdl.SWAN_NORA3(grid, start_time = '2018-08-25T00:00', end_time = '2018-08-25T01:00')

model.import_boundary(name='NORA3')
model.import_forcing(name='NORA3')

# All three properties optional
# Default filestring defined in defaults.py
# The defaults plotter grd_topo is defined in the ModelRun object and it adds a suffix "_topo" to any filename given to it.
model.plot_grid(save_fig=True, show_fig=False, filestring='#Grid_#Forcing.pdf')

# If one wants to use a different plotter
model.plot_grid(grid_plotter=dnplot.grd_mask())

# If we clear the boundary points in the grid, it will reflect to the ModelRun-object!
grid.set_boundary(grd.boundary.ClearBoundary())
model.plot_grid(grid_plotter=dnplot.grd_mask()) # << No boundary points visible!
