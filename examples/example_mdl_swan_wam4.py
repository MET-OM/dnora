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
grid = grd.Grid(lon=(5.0, 7.0), lat=(62.3, 62.9), name='Alesund300')

# Set spacing and boundary points
grid.set_spacing(dm=300)

# Import topography and mesh it down to the grid definitions
topo_reader=grd.read.KartverketNo50m(tile='*',folder='/home/konstac/bathy/')
grid.import_topo(topo_reader=topo_reader)
#
grid.mesh_grid()
#
# Set the boundaries
bnd_set = grd.boundary.EdgesAsBoundary(edges=['N', 'W', 'S'], step=20)
grid.set_boundary(boundary_setter=bnd_set)
#
#
# =============================================================================
# DEFINE MODEL OBJECT
# =============================================================================
model = mdl.SWAN(grid, start_time='2022-04-18T00:00',
                         end_time='2022-04-18T05:00')
# =============================================================================
# IMPORT BOUNDARIES AND FORCING
# =============================================================================
model.import_boundary(bnd.read_metno.WAM4km())
model.import_forcing(wnd.read_metno.MEPS())
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
model.run_model()
