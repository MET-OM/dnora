from point_skeleton import PointSkeleton
from gridded_skeleton import GriddedSkeleton
import numpy as np
from wave_spectrum import WaveSpectrum, WaveSpectrum2D
from coordinate_decorators import frequency, direction

class Spectra(PointSkeleton, WaveSpectrum):
    def __init__(self, grid=None, x=None, y=None, lon=None, lat=None, name: str="AnonymousSpectra"):
        self._name = name
        if grid is not None:
            x, y = grid.xy(strict=True)
            lon, lat = grid.lonlat(strict=True)
        self.data = self._create_structure(x, y, lon, lat, freq=np.array([0.1,0.2,0.3]))

class Boundary(PointSkeleton, WaveSpectrum2D):
    def __init__(self, grid=None, x=None, y=None, lon=None, lat=None, name: str="AnonymousSpectra"):
        self._name = name
        if grid is not None:
            x, y = grid.xy(strict=True)
            lon, lat = grid.lonlat(strict=True)
        self.data = self._create_structure(x, y, lon, lat, freq=np.array([0.1,0.2,0.3]), dirs=np.linspace(0,350,36))


class GriddedSpectra(GriddedSkeleton, WaveSpectrum):
    def __init__(self, grid=None, x=None, y=None, lon=None, lat=None, name: str="AnonymousSpectra"):
        self._name = name
        if grid is not None:
            x, y = grid.x(strict=True), grid.y(strict=True)
            lon, lat = grid.lon(strict=True), grid.lat(strict=True)
        self.data = self._create_structure(x, y, lon, lat, freq=np.array([0.1,0.2,0.3]))
