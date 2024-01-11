# =============================================================================
# IMPORT dnora
# =============================================================================
from dnora import grid, modelrun
from dnora.readers.generic_readers import ConstantGrid

# =============================================================================
# DEFINE GRID OBJECT
# =============================================================================
# Set grid definitions
# grid = grd.Grid(lon=(4.00, 5.73), lat=(60.53, 61.25), name='Skjerjehamn')
gridd = grid.Grid(lon=(4.00, 11.0), lat=(62.0, 65.0), name="sula_trondheim")


gridd.set_spacing(dm=1000)
gridd.set_boundary_points(grid.mask.Edges(edges=["N", "W", "S"]))
# grid.import_topo(grd.read.EMODNET2020(tile='C5'))
# grid.import_topo(grd.read.KartverketNo50m(folder='/home/janvb/Documents/Kartverket50m'))
# grid.mesh_grid()
# Create a ModelRun-object
model = modelrun.ModelRun(
    gridd, start_time="2018-08-25T00:00", end_time="2018-08-25T01:00", dry_run=False
)
model.import_wind(ConstantGrid())
breakpoint()
