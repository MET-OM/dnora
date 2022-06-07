import numpy as np
from abc import ABC, abstractmethod
from copy import copy
import scipy.ndimage as ndimage
# Import aux_funcsility functions
from .. import msg

class GridProcessor(ABC):
    """Abstract class for modifying bathymetrical data of the Grid-object."""
    @abstractmethod
    def __init__(self):
        pass

    def topo(self, data, lon, lat, land_sea_mask) -> np.ndarray:
        """Gets raw bathymetrical information in xyz-format and returns modified version.

        This method is called from within the Grid-object
        """
        return None

    def grid(self, data, lon, lat, land_sea_mask, boundary_mask) -> np.ndarray:
        """Gets meshed bathymetrical information and returns a modified version.

        This method is called from within the Grid-object
        """
        return None

    @abstractmethod
    def __str__(self):
        """
        Describes how the data is processed.

        This is called by the Grid-objeect to provide output to the user.
        """
        pass


class TrivialFilter(GridProcessor):
    """Returns the identical data it is passed. Used as default option."""

    def __init__(self):
        pass

    def topo(self, data, lon, lat, land_sea_mask):
        return copy(data)

    def grid(self, data, lon, lat, land_sea_mask, boundary_mask):
        return copy(data)

    def __str__(self):
        return("Doing nothing to the data, just passing it along.")

class SetMinDepth(GridProcessor):
    """Modify depth points shallower than a certain threshold.

    If to_land is True (default=False), then points shallower than min_depth
    are set to land. Otherwise the shallow points are set to min_depth.
    """

    def __init__(self, depth: float, to_land: bool=False) -> None:
        self.to_land = to_land
        self.depth = depth

    def topo(self, data, lon, lat, sea_mask):
        return self.grid(data, lon, lat, land_sea_mask)

    def grid(self, data, lon, lat, sea_mask, boundary_mask=None):
        shallow_points = data < self.depth
        if self.to_land:
            new_data = copy(data)
            new_data[np.logical_and(shallow_points,sea_mask)] = np.nan # Don't touch land points
            msg.plain(f"Affected {np.count_nonzero(np.logical_and(shallow_points, sea_mask))} points")
        else:
            # Set points to the limiter
            new_data = copy(data)
            new_data[np.logical_and(shallow_points, sea_mask)] = self.depth # Don't touch land points
            msg.plain(f"Affected {np.count_nonzero(np.logical_and(shallow_points, sea_mask))} points")
        return new_data

    def __str__(self):
        if self.to_land:
            return(f"Setting points shallower than {self.depth} to land (NaN)")
        else:
            return(f"Setting points shallower than {self.depth} to {self.depth}")

class SetMaxDepth(GridProcessor):
    """Modify depth points deeper than a certain threshold.

    If to_land is True (default=False), then points deeper than depth
    are set to land. Otherwise the shallow points are set to depth.
    """

    def __init__(self, depth: float, to_land: bool=False) -> None:
        self.to_land = to_land
        self.depth = depth

    def topo(self, data, lon, lat, sea_mask):
        return self.grid(data, lon, lat, land_sea_mask)

    def grid(self, data, lon, lat, sea_mask, boundary_mask=None):
        deep_points = data > self.depth
        if self.to_land:
            new_data = copy(data)
            new_data[np.logical_and(deep_points,sea_mask)] = np.nan # Don't touch land points
            msg.plain(f"Affected {np.count_nonzero(np.logical_and(deep_points, sea_mask))} points")
        else:
            # Set points to the limiter
            new_data = copy(data)
            new_data[np.logical_and(deep_points, sea_mask)] = self.depth # Don't touch land points
            msg.plain(f"Affected {np.count_nonzero(np.logical_and(deep_points, sea_mask))} points")
        return new_data

    def __str__(self):
        if self.to_land:
            return(f"Setting points deeper than {self.depth} to land (nan)")
        else:
            return(f"Setting points deeper than {self.depth} to {self.depth}")

class SetConstantDepth(GridProcessor):
    """Set a constant depth for the grid. """

    def __init__(self, depth: float) -> None:
        self.depth = depth

    def grid(self, data, lon, lat, sea_mask, boundary_mask=None):
        land_mask = np.logical_not(sea_mask)
        new_data = np.ones((len(lat), len(lon)))*self.depth
        new_data[land_mask] = np.nan # Don't touch land points
        msg.plain(f"Affected {np.count_nonzero(sea_mask)} points")
        return new_data

    def __str__(self):
        return(f"Creating a grid with constant depth {self.depth}")

class GaussianFilter(GridProcessor):
    """Set a constant depth for the grid. """

    def __init__(self, sigma: float=0.5) -> None:
        self.sigma = sigma

    def grid(self, data, lon, lat, sea_mask, boundary_mask=None):
        land_mask = np.logical_not(sea_mask)
        data[land_mask] = 0
        new_data = ndimage.filters.gaussian_filter(data, sigma=self.sigma)
        new_data[land_mask] = np.nan # Don't touch land points
        return new_data

    def __str__(self):
        return(f"Applying a gaussian filter with standard deviation {self.sigma} grid points")
