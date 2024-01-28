# =============================================================================
# IMPORT dnora
# =============================================================================
from dnora import grd, mdl, ice
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
grid = grd.Grid(lon=(30.67, 34.24), lat=(79.92, 80.54), name='Kvitoy')

# Set spacing and boundary points
grid.set_spacing(dm=500)

# Import topography and mesh it down to the grid definitions
# Options for grd.readers: EMODNET2020, KartverketNo50m, GEBCO2021
grid.import_topo(grd.read.EMODNET(), tile='B*', folder='/home/konstantinosc/bathy', source='internal', year=2022)

# This can be used to get an empty topography for testing
#grid.import_topo(topo_reader=grd.read.EmptyTopo(grid=grid))
#
grid.mesh_grid()
#
# Set the boundaries
bnd_set = grd.mask.Edges(edges=['N', 'W', 'S'])
grid.set_boundary_points(bnd_set)
#
#
# =============================================================================
# DEFINE MODEL OBJECT (mdl.)
# =============================================================================
# Options for mdl: SWAN_NORA3, SWAN_ERA5, SWAN_WW3_4km, SWAN_WAM4km
model = mdl.NORA3(grid, start_time='2021-11-21T12:00',
                               end_time='2021-11-21T13:00')
# =============================================================================
# IMPORT BOUNDARIES AND FORCING
# =============================================================================
# Boundary Spectra
#model.import_boundary(write_cache=True)
#model.import_boundary(bnd.read_metno.NORA3(source='lustre'),write_cache=True, read_cache=False) # for internal

# Wind Forcing
model.import_forcing(write_cache=True)
#model.import_forcing(wnd.read_metno.NORA3(source='lustre'),write_cache=True, read_cache=False)  # for internal

# Ocean Current
#model.import_oceancurrent(ocr.read_metno.NorKyst800())

# Water Level
#model.import_waterlevel(wlv.read_ec.GTSM_ERA5())

# Ice Forcing
model.import_ice()
#model.import_ice(ice.read_metno.Barents25())
# =============================================================================
# PLOT GRID, FORCING AND BOUNDARIES
# =============================================================================
#model.plot_grid(save_fig=True, show_fig=False)
# =============================================================================
# WRITE OUTPUT FOR SWAN RUN
# =============================================================================
model.export_grid()
#model.export_boundary()
#model.export_forcing()
#model.export_oceancurrent()
#model.export_waterlevel()
model.export_ice()
#model.write_input_file(input_file_writer=inp.SWAN(
#    spec_points=[(5.50, 59.16), (5.55, 59.15)]))
# =============================================================================
# SWAN RUN
# =============================================================================
model.run_model(model_executer = run.SWAN(nproc=15))
