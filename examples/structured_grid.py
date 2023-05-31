from dnora import grd, bnd, mdl

grid = grd.Grid(lon=[4.00, 5.73], lat=[60.53, 61.25], name='Skjerjehamn')
grid.set_spacing(dm=500)
#grid.set_utm(33,'N')

grid.set_boundary(grd.mask.LonLat(lon=4, lat=60.53))
model = mdl.NORA3(grid, start_time='2018-01-01', end_time='2018-01-02')

model.import_boundary(point_picker=bnd.pick.NearestGridPoint())