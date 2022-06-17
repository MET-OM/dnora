# =============================================================================
# IMPORT dnora
# =============================================================================
from dnora import grd, mdl, inp, bnd
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
grid = grd.Grid(lon=(109.150, 109.6), lat=(12.430, 12.815), name='Vietnam')

# Set spacing and boundary points
grid.set_spacing(dm=400)
grid.set_boundary(grd.boundary.MidPointAsBoundary(edges='E'))
# Import topography and mesh it down to the grid definitions
grid.import_topo(topo_reader=grd.read.GEBCO2021(tile='Vietnam',
                                                folder='/home/konstantinosc/bathy'))
grid.mesh_grid()

# =============================================================================
# DEFINE MODEL OBJECT
# =============================================================================
model = mdl.SWASH(grid, start_time='2018-08-25T00:00',
                  end_time='2018-08-25T05:00')

# =============================================================================
# PLOT GRID, FORCING AND BOUNDARIES
# =============================================================================
model.plot_grid(save_fig=True, show_fig=False)
# =============================================================================
# WRITE OUTPUT FOR SWASH RUN
# =============================================================================
model.export_grid()
model.write_input_file(input_file_writer=inp.SWASH(
    bound_side_command='BOU SIDE S CCW CON REG 0.5 20 180 '))  # Vietnam
# =============================================================================
# SWASH RUN
# =============================================================================
model.run_model(dry_run=False, mat_to_nc=True)
