# =============================================================================
# IMPORT dnora
# =============================================================================
import sys
dnora_directory = '../'
sys.path.insert(0, dnora_directory)
from dnora import trg, grd


import numpy as np
grid = trg.TrGrid(name='test')

grid.import_triang(trg.read_tr.SmsReader(filename='WW3_grid_01.2dm'))


#print(grid.boundary())

#bnd_set = trg.boundary.SetArray(np.arange(101))
bnd_set = trg.boundary.SetArray([70,72])

grid.append_boundary(bnd_set)

#print(grid.boundary())


#bnd_set = trg.boundary.SetArray([1,2,3])

#grid.set_boundary(bnd_set)

topo_reader = grd.EMODNET2018()
grid.import_topo(topo_reader)


grid.mesh_grid()

grid.plot_grid()
