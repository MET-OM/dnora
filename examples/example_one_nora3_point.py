# =============================================================================
# IMPORT dnora
# =============================================================================
from dnora import grd, mdl
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
grid = grd.Grid(lon=(4.00,5.00), lat=(60.53,61.53), name='Skjerjehamn')

# Create a ModelRun-object

# model = mdl.SWAN_NORA3(grid, start_time='2018-08-25T00:00', end_time='2018-08-25T06:00', dry_run=False)
#
# model.import_boundary(read_cache=False)

grid2 = grd.Grid(x=grid.edges('x'), y=grid.edges('y'))
grid2.set_utm(grid.utm()[0], grid.utm()[1])
model2 = mdl.SWAN_NORA3(grid2, start_time='2018-08-25T00:00', end_time='2018-08-25T06:00', dry_run=False)
model2.import_boundary(read_cache=False)
#model.import_spectra(read_cache=True, name='NORA3')
#model.import_forcing(read_cache=True, write_cache=False)
#model.boundary_to_waveseries()

#model.export_boundary()
#model.export_spectra()
#model.export_waveseries()
