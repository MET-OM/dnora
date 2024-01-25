# =============================================================================
# Read in a constant spectrum for one point
# =============================================================================
from dnora import grid, modelrun
from dnora.readers.generic_readers import ConstantPointData
import numpy as np
from dnora.spectral_conventions import SpectralConvention

point = grid.Grid(lon=4.00, lat=65.0)
point.set_boundary_points(grid.mask.All())
model = modelrun.ModelRun(
    point, start_time="2021-08-25T00:00", end_time="2021-08-25T05:00"
)
model.import_spectra(ConstantPointData(fp=0.1), spec=1.1)
model.import_spectra1d(
    ConstantPointData(freq=np.linspace(0.1, 0.5, 5), convention=SpectralConvention.MET),
    spec=1.1,
    dirm=45,
    spr=15,
)
print(model.spectra())
print(model.spectra1d())
