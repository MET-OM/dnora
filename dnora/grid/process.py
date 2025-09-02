from __future__ import annotations
import numpy as np
from abc import ABC, abstractmethod
from copy import copy
import scipy.ndimage as ndimage

from dnora import msg
from typing import Union, TYPE_CHECKING

if TYPE_CHECKING:
    from .grid import Grid, UnstrGrid


class GridProcessor(ABC):
    """Abstract class for modifying bathymetrical data of the Grid-object."""

    @abstractmethod
    def __call__(self, grid: Union[Grid, UnstrGrid]) -> np.ndarray:
        pass

    @abstractmethod
    def __str__(self):
        """
        Describes how the data is processed.

        This is called by the Grid-objeect to provide output to the user.
        """
        pass


class TrivialFilter(GridProcessor):
    """Returns the identical data it is passed. Used as default option."""

    def __call__(self, grid: Union[Grid, UnstrGrid], **kwargs) -> np.ndarray:
        return grid.topo()

    def __str__(self):
        return "Doing nothing to the data, just passing it along."


class SetMinDepth(GridProcessor):
    """Modify depth points shallower than a certain threshold.

    If to_land is True (default=False), then points shallower than min_depth
    are set to land. Otherwise the shallow points are set to min_depth.
    """

    def __init__(self, min_depth: float = 2.0, to_land: bool = False):
        self.min_depth = float(min_depth)
        self.to_land = to_land

    def __call__(
        self,
        grid: Union[Grid, UnstrGrid],
        min_depth: float = None,
        to_land: bool = None,
        ignore_land_mask: bool = False,
        **kwargs,
    ):

        if min_depth is None:
            min_depth = self.min_depth

        if to_land is None:
            to_land = self.to_land

        data = grid.topo()

        shallow_points = data < min_depth
        if ignore_land_mask:
            mask = shallow_points
        else:
            mask = np.logical_and(shallow_points, grid.sea_mask())

        if to_land:
            new_value = np.nan
        else:
            new_value = min_depth

        new_data = copy(data)
        new_data[mask] = new_value

        msg.plain(f"Affected {np.count_nonzero(mask)} points")
        return new_data

    def __str__(self):
        if self.to_land:
            return f"Setting shallow points to land (NaN)"
        else:
            return f"Setting minimum depth"


class SetMaxDepth(GridProcessor):
    """Modify depth points deeper than a certain threshold.

    If to_land is True (default=False), then points deeper than depth
    are set to land. Otherwise the shallow points are set to depth.
    """

    def __init__(self, max_depth: float = 999.0):
        self.max_depth = float(max_depth)

    def __call__(
        self,
        grid: Union[Grid, UnstrGrid],
        max_depth: float = None,
        **kwargs,
    ):
        data = grid.topo()

        if max_depth is None:
            max_depth = self.max_depth

        deep_points = data > max_depth

        new_data = copy(data)
        new_data[np.logical_and(deep_points, grid.sea_mask())] = (
            max_depth  # Don't touch land points
        )
        msg.plain(
            f"Affected {np.count_nonzero(np.logical_and(deep_points, grid.sea_mask()))} points"
        )
        return new_data

    def __str__(self):
        return f"Capping to a max depth"


class SetConstantDepth(GridProcessor):
    """Set a constant depth for the grid."""

    def __init__(self, depth: float = 999.0):
        self.depth = float(depth)

    def __call__(
        self,
        grid: Union[Grid, UnstrGrid],
        **kwargs,
    ):
        new_data = np.full(grid.size(), self.depth)

        new_data[grid.land_mask()] = np.nan  # Don't touch land points
        msg.plain(f"Affected {np.count_nonzero(grid.sea_mask())} points")
        return new_data

    def __str__(self):
        return f"Creating a grid with constant depth {self.depth}"


class GaussianFilter(GridProcessor):
    """Gaussian filter. Only works for structured grid."""

    def __call__(
        self,
        grid: Union[Grid, UnstrGrid],
        sigma: float = 0.5,
        **kwargs,
    ):
        self.sigma = sigma
        if not grid.is_gridded():
            msg.info(
                f"GaussianFilter only implemented for structured grids. Doing nothing."
            )
            return grid.topo()

        new_data = grid.topo()
        new_data[grid.land_mask()] = 0
        new_data = ndimage.filters.gaussian_filter(new_data, sigma=sigma)
        new_data[grid.land_mask()] = np.nan  # Don't touch land points

        return new_data

    def __str__(self):
        return f"Applying a gaussian filter with standard deviation {self.sigma} grid points"
