# =============================================================================
# IMPORT dnora
# =============================================================================
import sys
dnora_directory = '../'
sys.path.insert(0, dnora_directory)
from dnora import trg, grd, mdl

import numpy as np

grid = trg.Grid(name='test')


grid.import_triang(trg.read_tr.SmsReader(filename='WW3_grid_01.2dm'))

bnd_set = trg.boundary.SetArray([105, 157, 213, 271, 334, 335, 399, 469, 470, 541, 618])
grid.append_boundary(bnd_set)

topo_reader = grd.EMODNET2018()
#topo, topo_lon, topo_lat = topo_reader(lon_min = min(grid.lon()), lon_max = max(grid.lon()), lat_min = min(grid.lat()), lat_max = max(grid.lat()))
#topo_reader = grd.read.EmptyTopo(nx=10, ny=10)
grid.import_topo(topo_reader)
grid.mesh_grid()
#grid.plot_grid()

#print(grid)
#grid.export_grid(trg.write.WW3())
#grid.export_grid(grd.write.BoundaryPoints())
model = mdl.ModelRun(grid, start_time='2018-08-25T00:00', end_time='2018-08-25T03:00')
model.export_grid(trg.write.WW3())
