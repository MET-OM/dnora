# REEF3D
from dnora import grd, bnd, mdl

grid = grd.Grid(lon=(5.0, 5.2), lat=(58., 58.1))

grid.import_topo(grd.read.EMODNET2020(tile='*'))

#grid.set_spacing(dm=1000)
#grid.mesh_grid()

grid.set_boundary(grd.boundary.SetAll())
grid.set_boundary(boundary_setter=grd.boundary.MidPointAsBoundary(edges=['W']))

model = mdl.REEF3D(grid, start_time='2019-03-05T06:00', end_time='2019-03-05T06:00')
model.import_boundary(bnd.read_metno.NORA3())
model.boundary_to_spectra()
model.export_spectra()

model.export_grid(grd.write.REEF3D(use_raw=True))  # Use when not meshed
#model.export_grid(grd.write.REEF3D())  # Use when meshed
