# =============================================================================
# IMPORT dnora
# =============================================================================
import sys
dnora_directory = '../'
sys.path.insert(0, dnora_directory)
from dnora import trg, grd

import numpy as np

grid = trg.TrGrid(name='test')
print(grid)

grid.import_triang(trg.read_tr.SmsReader(filename='WW3_grid_01.2dm'))

bnd_set = trg.boundary.SetArray([70,72])
grid.append_boundary(bnd_set)
print(grid)

#topo_reader = grd.EMODNET2018()
topo_reader = grd.read.EmptyTopo(nx=10, ny=10)
grid.import_topo(topo_reader)
print(grid)
grid.mesh_grid()
#grid.plot_grid()

print(grid)
