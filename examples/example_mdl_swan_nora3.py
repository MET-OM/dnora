# =============================================================================
# IMPORT dnora
# =============================================================================
import sys
dnora_directory = '../'
sys.path.insert(0, dnora_directory)
from dnora import bnd, grd, mdl, wnd
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
lon_min = 4.0
lat_min = 60.53
lon_max = 5.73
lat_max = 61.25
grid = grd.Grid(lon_min, lon_max, lat_min, lat_max, name='Skjerjehamn250')
# Set spacing and boundary points
grid.set_spacing(dm=1000)
# Import topography and mesh it down to the grid definitions
grid.import_topo(topo_reader=grd.read.EMODNET2018(tile='D5',
                                                  folder='/home/konstantinosc/PhD/github/DNORA/bathy/'))
grid.mesh_grid()
# Set the boundaries
bnd_set = grd.boundary.EdgesAsBoundary(edges = ['N', 'W', 'S'], step = 10)
grid.set_boundary(boundary_setter = bnd_set)


model = mdl.SWAN_NORA3(grid, start_time = '2018-08-25T00:00', end_time = '2018-08-25T01:00')

# =============================================================================
# IMPORT BOUNDARIES AND FORCING
# =============================================================================
model.import_boundary(name='NORA3')
model.import_forcing(name='NORA3')
# =============================================================================
# WRITE OUTPUT FOR SWAN RUN
# =============================================================================
model.export_grid()
model.export_boundary()
model.export_forcing()
model.write_input_file()
# =============================================================================
# SWASH RUN
# =============================================================================
model.run_model()
