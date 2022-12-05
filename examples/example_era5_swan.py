# =============================================================================
# IMPORT dnora
# =============================================================================
from dnora import grd, mdl, wnd, bnd
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
grid = grd.Grid(lon=(4.00, 5.73), lat=(60.53, 61.25), name='Skjerjehamn250')

# Set spacing and boundary points
grid.set_spacing(dm=250)

# Import topography and mesh it down to the grid definitions
topo_reader=grd.read.EMODNET2020(tile='*',folder='/home/janvb/Documents/EMODNET2020/')
grid.import_topo(topo_reader=topo_reader)

# This can be used to get an empty topography for testing
#grid.import_topo(topo_reader=grd.read.EmptyTopo(grid=grid))
#
grid.mesh_grid()
#
# Set the boundaries
bnd_set = grd.boundary.EdgesAsBoundary(edges=['N', 'W', 'S'])
grid.set_boundary(boundary_setter=bnd_set)
#
#
# =============================================================================
# DEFINE MODEL OBJECT
# =============================================================================
model = mdl.SWAN_ERA5(grid, start_time='2018-08-25T00:00',
                              end_time='2018-08-26T00:00')
# =============================================================================
# IMPORT BOUNDARIES AND FORCING
# =============================================================================
# Use write_cache the first time
model.import_boundary(write_cache=True, read_cache=True)
model.import_forcing(write_cache=True, read_cache=True)
# =============================================================================
# PLOT GRID, FORCING AND BOUNDARIES
# =============================================================================
model.plot_grid(save_fig=True, show_fig=False)
# =============================================================================
# WRITE OUTPUT FOR SWAN RUN
# =============================================================================
model.export_grid()
model.export_boundary()
model.export_forcing()
model.write_input_file()
# =============================================================================
# SWAN RUN
# =============================================================================
#model.run_model()
