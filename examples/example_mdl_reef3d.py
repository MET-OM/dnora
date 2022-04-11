# =============================================================================
# IMPORT dnora
# =============================================================================
import sys
dnora_directory = '../'
sys.path.insert(0, dnora_directory)
from dnora import grd, mdl, bnd
import dnora.wave_parameters as wp

grid = grd.Grid(lon=(5.95, 6.15), lat=(62.35, 62.47), name='Sula')
grid.set_spacing(dm=50)
grid.import_topo(grd.read.KartverketNo50m(tile='*'),
                                    folder='/home/konstac/bathy/'))

grid.mesh_grid()
grid.set_boundary(boundary_setter=grd.boundary.MidPointAsBoundary(edges=['W']))


model = mdl.REEF3D(grid, start_time='2020-01-19T05:00', end_time='2020-01-20T05:00')
model.import_boundary(bnd.read_metno.NORA3())
model.boundary_to_spectra()

#model.plot_topo(save_fig=True)
model.export_spectra()
model.export_grid(grd.write.REEF3D(use_raw=True))  # Use when not meshed
#model.export_grid(grd.write.REEF3D())  #S Use when meshed

h = wp.Hs()(model.spectra().data)
t = wp.Tm_10()(model.spectra().data)

print(f'Hs={h.hs.values[0][0]}, Tm_10={t.tm_10.values[0][0]}')
