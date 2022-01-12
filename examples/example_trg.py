# =============================================================================
# IMPORT dnora
# =============================================================================
import sys
dnora_directory = '../'
sys.path.insert(0, dnora_directory)
from dnora import trg, grd, mdl

import numpy as np

# =============================================================================
# CREATE triangular grid object
# =============================================================================
grid = trg.Grid(name='test')

grid.import_triang(trg.read_tr.SmsReader(filename='WW3_grid_01.2dm'))

# Set additional boundary points
bnd_set = trg.boundary.SetArray([105, 157, 213, 271, 334, 335, 399, 469, 470, 541, 618])
grid.append_boundary(bnd_set)

# =============================================================================
# READ AND MESH topography
# =============================================================================
topo_reader = grd.EMODNET2018()
grid.import_topo(topo_reader)

grid.mesh_grid()

# =============================================================================
# EXPORT grid in WW3 format
# =============================================================================
model = mdl.ModelRun(grid, start_time='2018-08-25T00:00', end_time='2018-08-25T03:00')
model.export_grid(trg.write.WW3())
