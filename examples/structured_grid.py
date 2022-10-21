from dnora import grd, mdl

grid = grd.Grid(lon=(4.00, 5.73), lat=(60.53, 61.25), name='Skjerjehamn')
grid.set_spacing(dm=500)

grid.import_topo(grd.read.EMODNET2020(tile='D5'))
grid.mesh_grid()

grid.set_mask(grd.boundary.EdgesAsBoundary(edges=['N', 'W', 'S']))

model = mdl.WW3(grid)
model.plot_grid()
#model.export_grid()
