from abc import ABC, abstractmethod
from copy import copy
import numpy as np
from scipy.interpolate import griddata

class Mesher(ABC):
    """Abstract class for meshing the bathymetrical data to the grid."""
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, data, lon, lat, lonQ, latQ):
        """Gets the bathymetrical information and returns a version that is
        meshed to the area and resolution of the grid.

        This method is called from within the Grid-object
        """
        pass

    @abstractmethod
    def __str__(self):
        """Describes how the data is meshed.

        This is called by the Grid-objeect to provide output to the user.
        """
        pass


class Interpolate(Mesher):
    """Interpolates data to grid. A wrapper for scipy.interpolate's griddate."""
    def __init__(self, method = 'linear'):
        self.method = method

        return

    def __call__(self, data, lon, lat, lonQ, latQ):
        lon0, lat0 = np.meshgrid(lon, lat)
        lon1, lat1 = np.meshgrid(lonQ, latQ)
        data[np.isnan(data)] = 0 # Keeping land points as nan lets the shoreline creep out
        M = np.column_stack((data.ravel(), lon0.ravel(),lat0.ravel()))
        meshed_data = griddata(M[:,1:], M[:,0], (lon1, lat1), method=self.method)
        meshed_data[meshed_data>=0] = 32767

        return meshed_data

    def __str__(self):
        return(f"Meshing using {self.method} interpolation.")


class TrivialMesher(Mesher):
    """Passes along data. NB! This might not fit the grid, and is only used
    for e.g. recreating a Grid-object from an ouput file.
    """
    def __init__(self):
        pass

    def __call__(self, data, lon, lat, lonQ, latQ):
        return copy(data)

    def __str__(self):
        return("Passing input data along as final meshed grid.")
