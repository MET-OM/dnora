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


model = mdl.SWAN_NORA3(grid, start_time = '2018-08-25T00:00', end_time = '2018-08-25T01:00')

model.export_grid()

model.import_boundary(name='NORA3')
model.import_forcing(name='NORA3')

model.export_boundary()
model.export_forcing()

# Uses the folder defined in defaults.py
model.write_input_file()
# Assums the model is located in the same folder as the input file was written to
# This can be overridden with normal folder='#Forcing_#Grid' etc.
model.run_model()
