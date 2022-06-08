# =============================================================================
# IMPORT dnora
# =============================================================================
from dnora import grd, mdl, bnd, inp, run
import dnora.wave_parameters as wp

nproc=12 # number of processors for mpirun
#grid = grd.Grid(lon=(5.95, 6.15), lat=(62.35, 62.47), name='Sula20')
#grid = grd.Grid(lon=(4.93, 5.34), lat=(62.10, 62.24), name='StadA20')
#grid = grd.Grid(lon=(3.150, 3.319), lat=(56.514, 56.577), name='Ekofisk20')
grid = grd.Grid(lon=(4.84, 4.92), lat=(59.27, 59.33), name='Utsira20')


grid.set_spacing(dm=20)
#grid.import_topo(grd.read.EMODNET2020(tile='*',folder='/home/konstac/bathy/'))
grid.import_topo(grd.read.KartverketNo50m(tile='*',folder='/home/konstac/bathy/'))

grid.mesh_grid()
edges=['W']
grid.set_boundary(boundary_setter=grd.boundary.MidPointAsBoundary(edges=edges))


model = mdl.REEF3D(grid, start_time='2007-11-09T00:00', end_time='2007-11-09T01:00')
model.import_boundary(bnd.read_metno.NORA3())
model.boundary_to_spectra()

model.plot_grid(save_fig=True)
model.export_spectra()
model.export_grid(grd.write.REEF3D(use_raw=True))  # Use when not meshed
#model.export_grid(grd.write.REEF3D())  #S Use when meshed

model.write_input_file(input_file_writer=inp.REEF3D(option='DiveMESH', edges=edges, nproc=nproc))
model.write_input_file(input_file_writer=inp.REEF3D(option='REEF3D',nproc=nproc))
h = wp.Hs()(model.spectra().data)
t = wp.Tm_10()(model.spectra().data)

print(f'Hs={h.hs.values[0][0]}, Tm_10={t.tm_10.values[0][0]}')
model.run_model(model_executer = run.REEF3D(nproc=nproc))
