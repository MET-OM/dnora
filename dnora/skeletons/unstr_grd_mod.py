from point_skeleton import PointSkeleton
from topography import Topography
import numpy as np
import xarray as xr
from mask_factory import add_mask
from datavar_factory import add_datavar

@add_datavar(name='topo', default_value=999.)
@add_mask(name='boundary', coords='grid', default_value=0)
@add_mask(name='sea', coords='grid', default_value=1)
class UnstrGrid(PointSkeleton):
    def __init__(self, grid=None, x=None, y=None, lon=None, lat=None, name='AnonymousGrid'):
        self.name = name
        if grid is not None:
            x, y = grid.xy(strict=True)
            lon, lat = grid.lonlat(strict=True)
        self._init_structure(x, y, lon, lat)

    def _update_sea_mask(self):
        self.ds_manager.update_mask('sea', np.logical_not(np.isnan(self.topo())).astype(int))
