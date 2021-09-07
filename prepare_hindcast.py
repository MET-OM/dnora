from dnora2 import grd, bnd, wnd, inp, spec, run
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================

# Set grid definitions
lon_min = 4.93
lat_min = 58.7
lon_max = 6.22
lat_max = 59.45
grid = grd.Grid(lon_min, lon_max, lat_min, lat_max, name='Boknafjorden250')

# Set spacing and boundary points
grid.set_spacing(dm=250)
grid.set_boundary(bounN=30, edges=['N', 'W', 'S'])

# Import topography and mesh it down to the grid definitions
topo_fetcher = grd.TopoEMODNET2018(tile='D5',
                                   folder='/home/konstantinosc/PhD/github/DNORA/bathy/')
grid.import_topo(topo_fetcher)
grid.mesh_grid()


# Plot grid
grid.plot(save_fig=True)
grid.plot_mask()

# =============================================================================
# # Provide dates
start_time = '2020-01-13T18:00'
end_time = '2020-01-14T18:00'
#
# =============================================================================
# DEFINE BOUNDARY OBJECT
# =============================================================================

# Initialize boundary object with the grid
boundary = bnd.Boundary(grid, name='NORA3_boundary')

# Fetch the boundary spectra
boundary_fetcher = bnd.BoundaryNORA3()
point_picker = bnd.NearestGridPointPicker()

boundary.import_boundary(start_time, end_time, boundary_fetcher, point_picker)

# =============================================================================
# DEFINE WIND FORCING OBJECT
# =============================================================================
forcing = wnd.Forcing(grid, name='MEPS')
forcing_fetcher = wnd.ForcingMEPS()

# Fetch the wind forcing
forcing.import_forcing(start_time, end_time, forcing_fetcher)

# =============================================================================
# WRITE OUTPUT FOR SWAN RUN
# ===============================================================================

# Grid
write_grid = grd.OutputModelSWAN()
write_grid(grid)

# Boundary
boundary.process_spectra([spec.InterpSpectralProcessor(
    first_dir=0), spec.TrivialSpectralProcessor(calib_spec=1)])
write_boundary = bnd.OutputSWANascii()
write_boundary(boundary)

# Wind forcing
write_output = wnd.OutputSWANascii()
write_output(forcing)

# Write input file for SWAN model run
# =============================================================================
swan_directory = '/home/konstantinosc/Programs/swan4120'
write_input_file = inp.SWANInputFile(grid,forcing)
input_file_name = write_input_file(
    start_time, end_time, swan_directory=swan_directory, wind=True)
# =============================================================================

# =============================================================================
# SWAN RUN
# =============================================================================
run.run_SWAN(input_file_name, swan_directory=swan_directory)
