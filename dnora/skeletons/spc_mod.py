from point_skeleton import PointSkeleton
import numpy as np
from wave_spectrum import WaveSpectrum, WaveSpectrum2D
from coordinate_decorators import frequency, direction

class Spectra(PointSkeleton, WaveSpectrum):
    def __init__(self, x=None, y=None, lon=None, lat=None, name: str="AnonymousSpectra"):
        self._name = name
        self.data = self._create_structure(x, y, lon, lat, freq=np.array([0.1,0.2,0.3]))

class Boundary(PointSkeleton, WaveSpectrum2D):
    def __init__(self, x=None, y=None, lon=None, lat=None, name: str="AnonymousSpectra"):
        self._name = name
        self.data = self._create_structure(x, y, lon, lat, freq=np.array([0.1,0.2,0.3]), dirs=np.linspace(0,350,36))
