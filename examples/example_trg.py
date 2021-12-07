# =============================================================================
# IMPORT dnora
# =============================================================================
import sys
dnora_directory = '../'
sys.path.insert(0, dnora_directory)
from dnora import trg


import numpy as np
trgrid = trg.TrGrid(name='test')

trgrid.import_triang(trg.read_tr.SmsReader(filename='WW3_grid_01.2dm'))


print(trgrid.boundary())

#bnd_set = trg.boundary.SetArray(np.arange(101))
bnd_set = trg.boundary.SetArray([70,72])

trgrid.append_boundary(bnd_set)

print(trgrid.boundary())


bnd_set = trg.boundary.SetArray([1,2,3])

trgrid.set_boundary(bnd_set)
