# =============================================================================
# IMPORT dnora
# =============================================================================
import sys
dnora_directory = '../'
sys.path.insert(0, dnora_directory)
from dnora import grd, mdl
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
grid = grd.Grid(lon=(4.00, 5.73), lat=(60.53, 61.25), name='Skjerjehamn300')

# Set spacing and boundary points
grid.set_spacing(dm=300)

# Import topography and mesh it down to the grid definitions
grid.import_topo(topo_reader=grd.read.EMODNET2018(tile='D5',
                                    folder='/home/konstantinosc/bathy/'))
grid.mesh_grid()

# Set the boundaries
bnd_set = grd.boundary.EdgesAsBoundary(edges=['N', 'W', 'S'], step=20)
grid.set_boundary(boundary_setter=bnd_set)


# =============================================================================
# DEFINE MODEL OBJECT
# =============================================================================
model = mdl.SWAN_NORA3(grid, start_time='2018-08-25T00:00',
                             end_time='2018-08-25T03:00')
# =============================================================================
# IMPORT BOUNDARIES AND FORCING
# =============================================================================
model.import_boundary()
model.import_forcing()
# =============================================================================
# PLOT GRID, FORCING AND BOUNDARIES
# =============================================================================
model.plot_grid()
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
