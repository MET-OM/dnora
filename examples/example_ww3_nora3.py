# =============================================================================
# IMPORT dnora
# =============================================================================
from dnora import grd, mdl
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
#grid = grd.Grid(lon=(4.00, 5.73), lat=(60.53, 61.25), name='Skjerjehamn')
grid = grd.Grid(lon=(4.00, 11.), lat=(62., 65.), name='sula_trondheim')


grid.set_spacing(dm=1000)
#grid.import_topo(grd.read.EMODNET2020(tile='C5'))
grid.import_topo(grd.read.KartverketNo50m(folder='/home/janvb/Documents/Kartverket50m'))
grid.mesh_grid()
breakpoint()
# Create a ModelRun-object
model = mdl.WW3_NORA3(grid, start_time='2018-08-25T00:00',
                       end_time='2018-08-25T01:00', dry_run=False)
model.import_boundary()
model.plot_grid()
breakpoint()

model.boundary_to_spectra()
model.spectra_to_waveseries()

model.import_forcing()

model.export_boundary()
model.export_spectra()
model.export_waveseries()

model.export_forcing()
