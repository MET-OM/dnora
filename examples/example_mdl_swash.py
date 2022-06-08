# =============================================================================
# IMPORT dnora
# =============================================================================
from dnora import grd, mdl, inp, bnd
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
grid = grd.Grid(lon=(4.84, 4.92), lat=(59.27, 59.33), name='Utsira')
#grid = grd.Grid(lon=(12.26, 12.34), lat=(44.78, 44.85), name='IsoladellAmore20')
#grid = grd.Grid(lon=(12.3, 12.425), lat=(45.20, 45.44), name='Venezia')


# Set spacing and boundary points
grid.set_spacing(dm=25)
grid.set_boundary(grd.boundary.MidPointAsBoundary(edges='N'))
# Import topography and mesh it down to the grid definitions
grid.import_topo(topo_reader=grd.read.EMODNET2020(tile='*',
                                    folder='/home/konstac/bathy/'))
grid.mesh_grid()

# =============================================================================
# DEFINE MODEL OBJECT
# =============================================================================
model = mdl.SWASH(grid, start_time='2018-08-25T00:00', end_time='2018-08-25T00:01')

# =============================================================================
# PLOT GRID, FORCING AND BOUNDARIES
# =============================================================================
model.plot_grid(save_fig=True,show_fig=False)
# =============================================================================
# WRITE OUTPUT FOR SWASH RUN
# =============================================================================
model.export_grid()
model.write_input_file(input_file_writer=inp.SWASH(
    bound_side_command='BOU SIDE N CCW CON REG 0.5 20 0 ')) # Utsira
#model.write_input_file(input_file_writer=inp.SWASH(
#    bound_side_command='BOU SIDE S CCW CON REG 0.5 6 180 ')) # IsoladellAmore
#model.write_input_file(input_file_writer=inp.SWASH(
#    bound_side_command='BOU SIDE S CCW CON REG 0.5 16 180 ')) # Venezia
# =============================================================================
# SWASH RUN
# =============================================================================
model.run_model(dry_run=False, mat_to_nc=True)
