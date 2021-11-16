from dnora import grd, bnd, wnd, inp, spec, run
import sys

# =============================================================================
# IMPORT dnora
# =============================================================================
dnora_directory = '/home/konstantinosc/github/dnora/'
sys.path.insert(0, dnora_directory)


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
grid.set_spacing(dm=250)
grid.set_boundary(bounN=1, edges=['N', 'W', 'S'])

# Import topography and mesh it down to the grid definitions
topo_fetcher = grd.TopoEMODNET2018()
grid.import_topo(topo_fetcher)
grid.mesh_grid()


# =============================================================================
# DEFINE BOUNDARY OBJECT
# =============================================================================

# Initialize boundary object with the grid
boundary = bnd.Boundary(grid)

# Fetch the boundary spectra
boundary_fetcher = bnd.BoundaryWAM4(ignore_nan=True)
point_picker = bnd.NearestGridPointPicker()

start_time = '2021-08-15T00:00'
end_time = '2021-08-18T00:00'
boundary.import_boundary(start_time, end_time, boundary_fetcher, point_picker)

# =============================================================================
# DEFINE WIND FORCING OBJECT
# =============================================================================

forcing = wnd.Forcing(grid, name='MEPS')
forcing_fetcher = wnd.ForcingMEPS(prefix='det', hours_per_file=66)

# # Fetch the wind forcing
forcing.import_forcing(start_time, end_time, forcing_fetcher)


# =============================================================================
# WRITE OUTPUT FOR SWAN RUN
# =============================================================================

# Grid
write_grid = grd.OutputModelSWAN()
write_grid(grid)

# Boundary
write_boundary = bnd.OutputSWANascii(grid)
write_boundary(boundary)

# Wind forcing
write_forcing = wnd.OutputSWANascii()
write_forcing(forcing)

# Write input file for SWAN model run
swan_directory = '/home/konstantinosc/Programs/swan4120'
write_input_file = inp.SWANInputFile(grid)
input_file_name = write_input_file(
    start_time, end_time, swan_directory=swan_directory, wind=False)


# =============================================================================
# SWAN RUN
# =============================================================================
run.run_SWAN(input_file_name, swan_directory=swan_directory)
