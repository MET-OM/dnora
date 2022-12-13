# =============================================================================
# IMPORT dnora
# =============================================================================
from dnora import grd, mdl, wnd, bnd, inp, wlv, ocr
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
grid = grd.Grid(lon=(5.35, 5.6), lat=(59.00, 59.17), name='Boknafjorden')

# Set spacing and boundary points
grid.set_spacing(dm=150)

# Import topography and mesh it down to the grid definitions
# Options for grd.readers: EMODNET2018, EMODNET2020, KartverketNo50m, GEBCO2021
topo_reader = grd.read.EMODNET2020(
    tile='*', folder='/home/konstac/bathy')
grid.import_topo(topo_reader=topo_reader)

# This can be used to get an empty topography for testing
#grid.import_topo(topo_reader=grd.read.EmptyTopo(grid=grid))
#
grid.mesh_grid()
#
# Set the boundaries
bnd_set = grd.boundary.EdgesAsBoundary(edges=[ 'W', 'S'])
grid.set_boundary(boundary_setter=bnd_set)
#
#
# =============================================================================
# DEFINE MODEL OBJECT (mdl.)
# =============================================================================
# Options for mdl: SWAN_NORA3, SWAN_ERA5, SWAN
model = mdl.SWAN_NORA3(grid, start_time='2018-01-14T12:00',
                               end_time='2018-01-14T18:00')
# =============================================================================
# IMPORT BOUNDARIES AND FORCING ()
# =============================================================================
model.import_boundary(write_cache=True, read_cache=False)
model.import_forcing(write_cache=True, read_cache=False)
#model.import_oceancurrent(ocr.read_metno.NorKyst800())
#model.import_waterlevel(wlv.read_ec.GTSM_ERA5())
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
#model.export_waterlevel()
model.write_input_file(input_file_writer=inp.SWAN(
    spec_points=[(5.50, 59.16), (5.55, 59.15)]))
# =============================================================================
# SWAN RUN
# =============================================================================
model.run_model()
