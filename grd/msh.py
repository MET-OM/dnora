from abc import ABC, abstractmethod
from copy import copy
import numpy as np
from .. import msg
from scipy.interpolate import griddata

class Mesher(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, data, lon, lat, lonQ, latQ):
        pass

class InterpolationMesher(Mesher):
    def __init__(self, method = 'linear'):
        self.method = method
        #msg.info(f"Initializing mesher with method: {self.method}")

        return

    def __call__(self, data, lon, lat, lonQ, latQ):
        lon0, lat0 = np.meshgrid(lon, lat)
        lon1, lat1 = np.meshgrid(lonQ, latQ)
        data[np.isnan(data)] = 0 # Keeping land points as nan lets the shoreline creep out
        M = np.column_stack((data.ravel(), lon0.ravel(),lat0.ravel()))
        meshed_data = griddata(M[:,1:], M[:,0], (lon1, lat1), method=self.method)
        meshed_data[meshed_data>=0] = 32767

        return meshed_data

class TrivialMesher(Mesher):
    def __init__(self):
        pass

    def __call__(self, data, lon, lat, lonQ, latQ):
        return copy(data)


