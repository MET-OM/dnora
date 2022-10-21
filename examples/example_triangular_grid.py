# =============================================================================
# IMPORT dnora
# =============================================================================
from dnora import grd, mdl
import numpy as np
# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
grid = grd.TriGrid(name='sula_trondheim_v5')
#grid.import_triang(grd.read_tr.TxtReader('sula_trondheim_v5_bathy.txt', boundary_points=np.linspace(0,129,130)))
#grid.import_triang(grd.read_tr.MshReader('sula_trondheim_v5_bathy_clean.msh'))
#grid.arange_triangulation(grd.tri_arangers.ReorganizeBoundary(two_first_nodes=[20235-1, 40504-1], number_of_nodes=200))
grid.import_triang(grd.read_tr.MshReader('output/ww3_Msfibluesv003_bathy.msh'))

#grid.arange_triangulation(grd.tri_arangers.RemoveTriangle([[22440-1, 48412-1,46818-1]]))
grid.arange_triangulation(grd.tri_arangers.RemoveTriangle([[51594-1, 26499-1,56968-1]]))
#grid.arange_triangulation(grd.tri_arangers.ClearBoundary())
grid.import_topo(grd.read.MshFile('sula_trondheim_v5_bathy_clean.msh'))
#grid.import_topo(grd.read.EMODNET2020(tile='*'))
grid.mesh_grid()
#grid.process_grid(grd.process.SetMinDepth(10, ignore_land_mask=True))
model = mdl.WW3_NORA3(grid, start_time='2018-08-25T00:00', end_time='2018-08-25T01:00')
model.export_grid(grd.write_trg.WW3(), filename='ww3_Msfibluesv003.msh')
model.import_boundary()
model.plot_grid()
#model.export_grid(grd.write_trg.WW3())
# Create a ModelRun-object
#model = mdl.OnePoint_NORA3(grid, start_time='2018-08-25T00:00',
#                       end_time='2018-08-25T12:00', dry_run=False)

#model.import_boundary()
#model.boundary_to_waveseries()

#model.export_boundary()
#model.export_spectra()
#model.export_waveseries()
