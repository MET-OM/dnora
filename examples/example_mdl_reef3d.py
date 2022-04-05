# REEF3D
from dnora import grd, bnd, mdl

grid = grd.Grid(lon=(3.0, 4.0), lat=(56., 57.))
grid.set_spacing(nx=5, ny=5)
grid.set_boundary(boundary_setter=grd.boundary.MidPointAsBoundary(edges=['N']))

#grid = grd.Grid(lon=(3.0, 3.0), lat=(56., 56.))
#grid.set_boundary(grd.boundary.SetAll())

model = mdl.REEF3D(grid, start_time='2019-03-05T06:00', end_time='2019-03-05T06:00')
model.import_boundary(bnd.read_metno.NORA3())

model.boundary_to_spectra()
model.export_spectra()
