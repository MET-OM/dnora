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
lon_min = 5.028305
lat_min = 59.4087
lon_max = 5.1560
lat_max = 59.4635
grid = grd.Grid(lon_min, lon_max, lat_min, lat_max, name='Roevaer5')

# Set spacing and boundary points
grid.set_spacing(dm=5)

# Import topography and mesh it down to the grid definitions
grid.import_topo(topo_reader=grd.read.EMODNET2018(tile='D5',
                                    folder='/home/konstantinosc/bathy/'))
grid.mesh_grid()

# =============================================================================
# DEFINE MODEL OBJECT
# =============================================================================
model = mdl.SWASH(grid, start_time='2018-08-25T00:00', end_time='2018-08-25T02:00')

# =============================================================================
# PLOT GRID, FORCING AND BOUNDARIES
# =============================================================================
model.plot_grid()

# =============================================================================
# WRITE OUTPUT FOR SWASH RUN
# =============================================================================
model.export_grid() 
model.write_input_file(input_file_writer=inp.SWASH(
    bound_side_command='BOU SIDE N CCW CON REG 0.5 14 0 '))
# =============================================================================
# SWASH RUN
# =============================================================================
model.run_model()


