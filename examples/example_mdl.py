# =============================================================================
# IMPORT dnora
# =============================================================================
import sys
dnora_directory = '../'
sys.path.insert(0, dnora_directory)
from dnora import grd, mdl, bnd, wnd
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
grid = grd.Grid(lon_min=4.0, lon_max=60.53, lat_min=5.73, lat_max=61.25, name='Skjerjehamn')
grid.set_spacing(dm=1000)

# Import topography and mesh it down to the grid definitions
grid.import_topo(topo_reader=grd.read.EMODNET2018(tile='D5'))
grid.mesh_grid()

# Set the boundaries
bnd_set = grd.boundary.EdgesAsBoundary(edges=['N', 'W', 'S'])
grid.set_boundary(boundary_setter=bnd_set)

# Create a ModelRun-object
model = mdl.ModelRun(grid, start_time='2018-08-25T00:00',
                       end_time='2018-08-25T01:00')

model.export_grid(grd.write.WW3())

model.import_boundary(bnd.read_metno.NORA3(), point_picker=bnd.pick.Area())
model.import_forcing(wnd.read_metno.NORA3())

model.export_boundary(bnd.write.WW3())
model.export_forcing(wnd.write.WW3())
