# =============================================================================
# IMPORT dnora
# =============================================================================
from dnora import grd, mdl, bnd, wnd
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
grid = grd.Grid(lon=(4.00, 5.73), lat=(60.53, 61.25), name='Skjerjehamn')
grid.set_spacing(dm=10000)

# Import topography and mesh it down to the grid definitions
#grid.import_topo(topo_reader=grd.read.EMODNET2020(tile='*'))
grid.mesh_grid()

# Set the boundaries
bnd_set = grd.boundary.EdgesAsBoundary(edges=['N', 'W', 'S'])
grid.set_boundary(boundary_setter=bnd_set)

# Create a ModelRun-object
model = mdl.WW3_NORA3(grid, start_time='2018-08-25T00:00',
                       end_time='2018-08-25T01:00')

model.export_grid()

model.import_boundary()
model.import_forcing()

model.export_boundary()
model.export_forcing()
