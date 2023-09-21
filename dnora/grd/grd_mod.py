from geo_skeletons import GriddedSkeleton, PointSkeleton
import numpy as np
from geo_skeletons.decorators import add_mask, add_datavar
from .grid_methods import GridMethods

from ..aux_funcs import read_ww3_info
from pathlib import Path

@add_datavar(name='topo', default_value=999.)
@add_mask(name='boundary', coords='grid', default_value=0)
@add_mask(name='output', coords='grid', default_value=0)
@add_mask(name='sea', coords='grid', default_value=1, opposite_name='land')
class Grid(GriddedSkeleton, GridMethods):
    @classmethod
    def from_ww3_grid(cls, gridname: str, folder: str=''):
        """Recreate a WW3 grid object based on the _info, _bathy and _mapsta files"""

        filename = Path(folder) / f"{gridname}_info.txt"

        print(filename)
        lon_min, lon_max, lat_min, lat_max, dlon, dlat, NX, NY = read_ww3_info(filename)

        filename = Path(folder) / f'{gridname}_bathy.txt'
        topo=np.loadtxt(filename).reshape((NY,NX))
        filename = Path(folder) / f'{gridname}_mapsta.txt'
        mask=np.loadtxt(filename).reshape((NY,NX)) == 2 # Boundary points given as value 2

        grid = cls(lon=(lon_min, lon_max), lat=(lat_min, lat_max), name=gridname)
        grid.set_spacing(nx = NX, ny = NY)
        grid.set_topo(topo)
        grid.set_boundary_mask(mask)

        return grid
    
    def __init__(self, x=None, y=None, lon=None, lat=None, name='LonelyGrid', **kwargs):
        if np.all([a is None for a in [x,y,lon,lat]]):
            x, y = 0, 0
        super().__init__(x, y, lon, lat, name, **kwargs)

    def boundary_nx(self) -> int:
        """Return approximate number of grid points in the longitude direction
        """
        abs_diff = np.abs(np.diff(np.where(self.boundary_mask())))
        if abs_diff.size == 0:
            return 0
        abs_diff=np.median(abs_diff[abs_diff>0]).astype(int)

        return np.ceil(self.nx()/abs_diff).astype(int)

    def boundary_ny(self) -> int:
        """Return approximate number of grid points in the longitude direction
        """
        abs_diff = np.abs(np.diff(np.where(self.boundary_mask())))
        if abs_diff.size == 0:
            return 0
        abs_diff=np.median(abs_diff[abs_diff>0]).astype(int)

        return np.ceil(self.ny()/abs_diff).astype(int)