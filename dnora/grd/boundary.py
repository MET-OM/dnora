import numpy as np
from typing import List
from abc import ABC, abstractmethod

# Import aux_funcsiliry functions
from .. import msg

class BoundarySetter(ABC):
    """Set boundary points in the grid.

    The dimensions and orientation of the boolean array [True = boundary point]
    that is returned to the object should be:

    rows = latitude and colums = longitude (i.e.) shape = (nr_lat, nr_lon).

    North = [-1,:]
    South = [0,:]
    East = [:,-1]
    West = [:,0]
    """
    @abstractmethod
    def __call__(self, sea_mask: np.ndarray) -> np.ndarray:
        """This method is called from within the Grid-object."""
        return boundary_mask

    @abstractmethod
    def __str__(self):
        """Describes how the boundary points are set.

        This is called by the Grid-objeect to provide output to the user.
        """
        pass


class ClearBoundary(BoundarySetter):
    """Clears all boundary points by setting a mask with False values."""

    def __init__(self):
        pass

    def __call__(self, sea_mask: np.ndarray):
        mask_size = sea_mask.shape
        return np.full(mask_size, False)

    def __str__(self):
        return("Clearing all possible boundary points and setting an empty mask.")


class EdgesAsBoundary(BoundarySetter):
    """Set the grid edges as boundary points.

    Any combination of North, South, East, West ['N', 'S', 'E', 'W'] edges
    can be set.

    If step is e.g. 5, then only every fifth point of the edges are set. This
    is useful if the boundary spectra are coarse and we want to let the wave
    model interpolate the spectra.
    """

    def __init__(self, edges: List[str]=['N', 'S', 'E', 'W'], step: int=1) -> None:
        self.edges = edges
        if step < 1:
            raise ValueError('step cannot be smaller than 1')
        else:
            self.step = int(step)
        return

    def __call__(self, sea_mask: np.ndarray):
        mask_size = sea_mask.shape

        if mask_size == (1, 1):
            return np.full(mask_size, True)

        boundary_mask = np.full(mask_size, False)
        # --------- North boundary ----------
        if 'N' in self.edges:
            boundary_mask[-1,::self.step] = True
        ## --------- South boundary ----------
        if 'S' in self.edges:
            boundary_mask[0,::self.step] = True
        ## --------- East boundary ----------
        if 'E' in self.edges:
            boundary_mask[::self.step,-1] = True
        ## --------- West boundary ----------
        if 'W' in self.edges:
            boundary_mask[::self.step,0] = True

        return boundary_mask

    def __str__(self):
        return(f"Setting all edges {self.edges} to boundary points using step {self.step}.")


class MidPointAsBoundary(BoundarySetter):
    """Set the middle point of grid edges as boundary points.

    Any combination of North, South, East, West ['N', 'S', 'E', 'W'] edges
    can be set.
    """

    def __init__(self, edges: List[str]=['N', 'S', 'E', 'W']) -> None:
        self.edges = edges


    def __call__(self, sea_mask: np.ndarray):
        mask_size = sea_mask.shape
        if mask_size == (1, 1):
            return np.full(mask_size, True)

        boundary_mask = np.full(mask_size, False)
        ny = np.round(mask_size[0]/2).astype(int)
        nx = np.round(mask_size[1]/2).astype(int)

        # --------- North boundary ----------
        if 'N' in self.edges:
            edge = sea_mask[-1,:]
            nx = np.round(np.median(np.where(edge))).astype(int)
            boundary_mask[-1,nx] = True
        ## --------- South boundary ----------
        if 'S' in self.edges:
            edge = sea_mask[0,:]
            nx = np.round(np.median(np.where(edge))).astype(int)
            boundary_mask[0,nx] = True
        ## --------- East boundary ----------
        if 'E' in self.edges:
            edge = sea_mask[:,-1]
            ny = np.round(np.median(np.where(edge))).astype(int)
            boundary_mask[ny,-1] = True
        ## --------- West boundary ----------
        if 'W' in self.edges:
            edge = sea_mask[:,0]
            ny = np.round(np.median(np.where(edge))).astype(int)
            boundary_mask[ny,0] = True

        return boundary_mask

    def __str__(self):
        return(f"Setting mid point of edges {self.edges} to boundary points.")


class SetMatrix(BoundarySetter):
    """Set boundary points by providing a boolean array [True = boundary point].

    The dimensions and orientation of the array should be:

    rows = latitude and colums = longitude (i.e.) shape = (nr_lat, nr_lon).

    North = [-1,:]
    South = [0,:]
    East = [:,-1]
    West = [:,0]
    """

    def __init__(self, matrix):
        self.matrix = matrix
        return

    def __call__(self, mask_size: tuple):
        if self.matrix.shape == mask_size:
            return self.matrix
        else:
            raise Exception(f'Given mask for boundary points does not match the dimensions of the grid ({self.matrix.shape[0]}x{self.matrix.shape[1]} vs {mask_size[0]}x{mask_size[1]})')

    def __str__(self):
        return(f"Setting boundary points using the boolean matrix I was initialized with.")

class SetAll(BoundarySetter):
    """Set all points to boundary points. """

    def __call__(self, sea_mask: np.ndarray):
        mask_size = sea_mask.shape
        return np.full(mask_size, True)

    def __str__(self):
        return(f"Setting all points to boundary points.")
