import sys

# =============================================================================
# IMPORT dnora
# =============================================================================
dnora_directory = '../'
sys.path.insert(0, dnora_directory)
from dnora import grd, bnd, wnd, inp, run

# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================

# Set grid definitions
lon_min = 4.93
lat_min = 58.7
lon_max = 6.22
lat_max = 59.45
grid = grd.Grid(lon_min, lon_max, lat_min, lat_max, name='Bokn')

# Set spacing and boundary points
grid.set_spacing(dm=1000)

bnd_set = grd.boundary.EdgesAsBoundary(edges=['N', 'W', 'S'], step=1)
#grid.set_boundary(bounN=1, edges=['N', 'W', 'S'])
grid.set_boundary(boundary_setter=bnd_set)

# Import topography and mesh it down to the grid definitions
#grid.import_topo(topo_reader=grd.read.EMODNET2018(tile='D5',
#                                                  folder='/home/konstantinosc/PhD/github/DNORA/bathy/'))
grid.import_topo(topo_reader=grd.read.EMODNET2018(tile='D5')

grid.mesh_grid()


# =============================================================================
# DEFINE BOUNDARY OBJECT
# =============================================================================

# Initialize boundary object with the grid
boundary = bnd.Boundary(grid, name='NORA3')

# Fetch the boundary spectra
time0 = '2020-01-13T18:00'
time1 = '2020-01-13T19:00'
boundary.import_boundary(start_time=time0, end_time=time1, boundary_reader=bnd.read.MetNo_NORA3(
), point_picker=bnd.pick.NearestGridPoint())

# =============================================================================
# DEFINE WIND FORCING OBJECT
# =============================================================================

forcing = wnd.Forcing(grid, name='NORA3')

# # Fetch the wind forcing
forcing.import_forcing(start_time=time0, end_time=time1,
                       forcing_reader=wnd.read.MetNo_NORA3())


# =============================================================================
# WRITE OUTPUT FOR SWAN RUN
# =============================================================================
output_folder = '/home/konstantinosc/test/'
swan_directory = '/home/konstantinosc/Programs/swan4120/'
# Grid
write_grid = grd.write.SWAN(folder=output_folder)
write_grid(grid)

# Boundary
write_boundary = bnd.write.SWAN(
    folder=output_folder, filestring='Bokn_specT0-T1', datestring='%Y%m%d')
write_boundary(boundary)

# Wind forcing
write_forcing = wnd.write.SWAN(
    folder=output_folder, filestring='Bokn_windT0-T1', datestring='%Y%m%d')(boundary)
write_forcing(forcing)

# Write input file for SWAN model run
write_input_file = inp.SWANInputFile(grid=grid, forcing=forcing)
input_file_name = write_input_file(start_time=time0, end_time=time1,
                                   swan_directory=swan_directory, forcing_folder=output_folder, wind=True)


# =============================================================================
# SWAN RUN
# =============================================================================
#run.run_SWAN(input_file_name, swan_directory=swan_directory)
