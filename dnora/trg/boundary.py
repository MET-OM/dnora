import numpy as np
from typing import List
from abc import ABC, abstractmethod
from copy import copy
import numpy as np
# Import auxiliry functions
from .. import msg

class BoundarySetter(ABC):
    """Set boundary points in the grid."""

    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, nodes):
        """This method is called from within the Grid-object."""

        return boundary_list

    @abstractmethod
    def __str__(self):
        """Describes how the boundary points are set.

        This is called by the Grid-object to provide output to the user.
        """
        pass


class ClearBoundary(BoundarySetter):
    def __init__(self):
        return

    def __call__(self, nodes):
        return []

    def __str__(self):
        return "Clearing all boundary points"


class SetArray(BoundarySetter):
    def __init__(self, boundary_array):
        self._boundary_array = np.array(boundary_array)
        return

    def __call__(self, nodes):
        bnd_points = self._boundary_array
        nodes = np.array(nodes)
        mask = np.logical_and(bnd_points <= max(nodes), bnd_points >= min(nodes))
        return bnd_points[mask]

    def __str__(self):
        return "Setting new boundary points based on an provided array"
