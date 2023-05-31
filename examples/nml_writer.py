from dnora import grd, mdl, inp

grid = grd.Grid(lon=[4.00, 5.73], lat=[60.53, 61.25], name='Skjerjehamn')
grid.set_spacing(dm=1000)
grid.set_boundary(grd.mask.LonLat(lon=4.0, lat=61.25))

model = mdl.WW3_NORA3(grid, start_time='2018-01-01 00:00:00', end_time='2018-01-01 02:00:00')
model.set_spectral_grid(freq0=0.02, nfreq=40, ndir=36, finc=1.1)
model.import_boundary()
model.export_grid(filename='test',dry_run=True)
model.write_input_file(inp.inp.WW3())
