# =============================================================================
# IMPORT dnora
# =============================================================================
import sys
dnora_directory = '../'
sys.path.insert(0, dnora_directory)
from dnora import bnd, grd, mdl, wnd

lon_min = 4.0
lat_min = 60.53
lon_max = 5.73
lat_max = 61.25

grid = grd.Grid(lon_min, lon_max, lat_min, lat_max, name='Skjerjehamn250')

grid.set_spacing(dm=1000)
bnd_set = grd.boundary.EdgesAsBoundary(edges = ['N', 'W', 'S'], step = 20)
grid.set_boundary(boundary_setter = bnd_set)


model = mdl.SWAN_NORA3()
model.set_grid(grid)

model.set_run_time(start_time = '2018-08-25T00:00', end_time = '2018-08-25T18:00')

model.import_boundary()
model.import_forcing()
model.export_boundary()
model.export_forcing()
