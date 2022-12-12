# =============================================================================
# IMPORT dnora
# =============================================================================
from dnora import grd, mdl, bnd, wnd
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
grid = grd.Grid(lon=(4.00, 5.73), lat=(60.53, 61.25), name='Skjerjehamn')
#grid = grd.Grid(lon=(5.724778-0.8, 5.724778+0.8), lat=(59.480911-0.4, 59.480911+0.4), name='sula_trondheim')
#grid = grd.Grid(lon=(105, 116), lat=(14, 23), name='Hanoitest')

grid.set_spacing(dm=1000)
#grid.import_topo(grd.read.EMODNET2020(tile='*5'))
#grid.import_topo(grd.read.KartverketNo50m(folder='/home/janvb/Documents/Kartverket50m'))
#grid.mesh_grid()

# Create a ModelRun-object
model = mdl.WW3_WAM4km(grid, start_time='2022-11-30T00:00',
                        end_time='2022-11-30T01:00', dry_run=False)
model.import_boundary(expansion_factor=1.2, read_cache=True)
model.boundary_to_waveseries(parameters=['hs'])
#model.cache_boundary()
# #model.export_boundary()
# #model.export_grid(grd.write.REEF3D())
#model.import_forcing(expansion_factor=1.1, write_cache=True)
#model.export_boundary()
# #model.plot_forcing()
# #model.import_boundary(read_cache=True)
# #model.cache_boundary()
# #model.import_forcing(wnd.read_metno.MEPS(), write_cache=True, read_cache=True)
# breakpoint()
# model.boundary_to_spectra(write_cache=True)
# model.plot_grid()
#
#
# model.boundary_to_spectra()
# model.spectra_to_waveseries()
#
# model.import_forcing()
#
# model.export_boundary()
# model.export_spectra()
# model.export_waveseries()
#
# model.export_forcing()
