# =============================================================================
# IMPORT dnora
# =============================================================================
import sys
dnora_directory = '../'
sys.path.insert(0, dnora_directory)
from dnora import grd, inp, run, dnplot
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================

# Set grid definitions
lon_min = 5.933
lat_min = 62.383
lon_max = 6.110
lat_max = 62.463
grid = grd.Grid(lon_min, lon_max, lat_min, lat_max, name='Sula')

# Set spacing and boundary points
grid.set_spacing(dm=20)

# Import topography and mesh it down to the grid definitions
grid.import_topo(topo_reader=grd.read.EMODNET2018(tile='C5',
                                                  folder='/home/konstantinosc/PhD/github/DNORA/bathy/'))


grid.mesh_grid()

# =============================================================================
# PLOT TOPOGRAPHY
# =============================================================================
dnplot.grd_topo(grid)

# =============================================================================
# WRITE OUTPUT FOR SWAN RUN
# =============================================================================
output_folder = '/home/konstantinosc/test/'
swash_directory = '/home/konstantinosc/Programs/swash/'

# Grid
grid.export_grid(grid_writer=grd.write.SWASH(folder=output_folder))




# Write input file for SWASH model run
write_input_file = inp.SWASHInputFile(grid=grid)
input_file_name = write_input_file(start_time='2020-01-13T18:00',
                                   end_time='2020-01-13T18:05',
                                   bound_side_command='BOU SIDE W CCW CON REG 0.5 14 270 ' ,
                                   folder=swash_directory)


# =============================================================================
# SWASH RUN
# =============================================================================
#run.run_SWASH(input_file_name, swash_directory=swash_directory)
