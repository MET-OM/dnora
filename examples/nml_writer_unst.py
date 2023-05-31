from dnora import grd, mdl, inp

grid = grd.TriGrid.from_msh(filename='ww3_Msfibluesv004.msh', name='SulaTrondheim')

model = mdl.WW3_NORA3(grid, start_time='2018-01-01 00:00:00', end_time='2018-01-01 02:00:00')
model.import_boundary()
model.write_input_file(inp.inp.WW3())
