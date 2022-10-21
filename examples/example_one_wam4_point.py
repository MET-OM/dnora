# =============================================================================
# IMPORT dnora
# =============================================================================
from dnora import grd, mdl, bnd
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
#grid = grd.Grid(lon=(4.7799916, 4.7799916), lat=(59.5, 59.5), name='Haugesund_SE')
#grid = grd.Grid(lon=(4.8999844, 4.8999844), lat=(59.55, 59.55), name='Haugesund_NW')
grid = grd.Grid(lon=4.8999844, lat=59.55, name='Haugesund_NW')

breakpoint()

# 59.55
# 4.8999844
# Create a ModelRun-object
model = mdl.OnePoint(grid, start_time='2022-10-05T00:00',
                       end_time='2022-10-07T00:00', dry_run=False)

model.import_boundary(bnd.read_metno.WAM4km())
#model.boundary_to_waveseries()

model.export_boundary()
#model.export_spectra()
#model.export_waveseries()
