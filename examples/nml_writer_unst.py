from dnora import grd, mdl, exp

grid = grd.TriGrid.from_msh(filename='ww3_Msfibluesv004.msh', name='SulaTrondheim')

model = mdl.ModelRun(grid, start_time='2018-01-01 00:00:00', end_time='2018-01-01 02:00:00')
exporter = exp.WW3(model)
#exporter.export_grid()