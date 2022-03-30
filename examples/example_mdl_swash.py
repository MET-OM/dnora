# =============================================================================
# IMPORT dnora
# =============================================================================
import sys
dnora_directory = '../'
sys.path.insert(0, dnora_directory)
from dnora import grd, mdl, inp
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
grid = grd.Grid(lon=(12.25, 12.39), lat=(44.75, 44.85), name='Gorino25')

# Set spacing and boundary points
grid.set_spacing(dm=25)

# Import topography and mesh it down to the grid definitions
grid.import_topo(topo_reader=grd.read.EMODNET2020(tile='E5',
                                    folder='/home/konstac/bathy/'))
grid.mesh_grid()

# =============================================================================
# DEFINE MODEL OBJECT
# =============================================================================
model = mdl.SWASH(grid, start_time='2018-08-25T00:00', end_time='2018-08-25T02:00')

# =============================================================================
# PLOT GRID, FORCING AND BOUNDARIES
# =============================================================================
model.plot_grid(save_fig=True,show_fig=False)

# =============================================================================
# WRITE OUTPUT FOR SWASH RUN
# =============================================================================
model.export_grid()
model.write_input_file(input_file_writer=inp.SWASH(
    bound_side_command='BOU SIDE S CCW CON REG 1.0 20 180 '))
# =============================================================================
# SWASH RUN
# =============================================================================
model.run_model()
