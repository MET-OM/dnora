from point_skeleton import PointSkeleton
from topography import Topography
import numpy as np
import xarray as xr
#from ..grd.read import TopoReader
from dnora import msg


class UnstrGrid(PointSkeleton, Topography):
    def __init__(self, grid=None, x=None, y=None, lon=None, lat=None, name='AnonymousGrid'):
        self.name = name
        if grid is not None:
            x, y = grid.xy(strict=True)
            lon, lat = grid.lonlat(strict=True)
        self.data = self._create_structure(x, y, lon, lat)
        self._reset_vars()
