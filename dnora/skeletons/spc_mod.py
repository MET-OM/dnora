from point_skeleton import PointSkeleton
from gridded_skeleton import GriddedSkeleton
import numpy as np
from power_spectrum import PowerSpectrum
from coordinates import add_time, add_frequency, add_direction
from datetime import datetime, timedelta
from mask_generator import add_mask

@add_frequency(grid_coord=False)
@add_time(grid_coord=True)
class Spectra(PointSkeleton, PowerSpectrum):
    def __init__(self, grid=None, x=None, y=None, lon=None, lat=None, name: str="AnonymousSpectra"):
        self._name = name
        if grid is not None:
            x, y = grid.xy(strict=True)
            lon, lat = grid.lonlat(strict=True)
        self.data = self._create_structure(x, y, lon, lat, freq=np.array([0.1,0.2,0.3]))
        self._reset_vars()

@add_direction(grid_coord=False, coord_name='dirs')
@add_frequency(grid_coord=False)
@add_time(grid_coord=False)
class Boundary(PointSkeleton, PowerSpectrum):
    def __init__(self, grid=None, x=None, y=None, lon=None, lat=None, name: str="AnonymousSpectra"):
        self._name = name
        if grid is not None:
            x, y = grid.xy(strict=True)
            lon, lat = grid.lonlat(strict=True)
        t = np.arange(datetime(2020,1,1), datetime(2020,1,2), timedelta(hours=1)).astype(datetime)
        self.data = self._create_structure(x, y, lon, lat, time=t, freq=np.array([0.1,0.2,0.3]), dirs=np.array([0, 45, 90, 135, 180, 225, 270, 315]))
        self._init_masks()
        #self._reset_vars()
