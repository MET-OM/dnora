# =============================================================================
# IMPORT dnora
# =============================================================================
from dnora import grd, mdl
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
grid = grd.Grid(lon=4.00, lat=60.53, name='Skjerjehamn')

# Create a ModelRun-object
model = mdl.OnePoint_NORA3(grid, start_time='2018-08-25T00:00',
                       end_time='2018-08-25T12:00', dry_run=False)

model.import_boundary()
model.boundary_to_waveseries()

#model.export_boundary()
#model.export_spectra()
#model.export_waveseries()
