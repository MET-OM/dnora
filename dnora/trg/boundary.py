import numpy as np
from typing import List
from abc import ABC, abstractmethod
from copy import copy
import numpy as np
from typing import Iterable
# Import aux_funcsiliry functions
from .. import msg

class BoundarySetter(ABC):
    """Set boundary points in the grid."""

    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, nodes, tri, lon, lat):
        """This method is called from within the Grid-object."""

        return boundary_list, triangulation, nodes, longitude, latitude

    @abstractmethod
    def __str__(self):
        """Describes how the boundary points are set.

        This is called by the Grid-object to provide output to the user.
        """
        pass

class ReorganizeBoundary(BoundarySetter):
    @staticmethod
    def find_boundary_nodes(tri: np.ndarray, two_first_nodes: Iterable, number_of_nodes: int) -> np.ndarray:
        """Finds all the consequtive boundary nodes when given the seed of the
        two first ones"""
        boundary_nodes = np.zeros(number_of_nodes).astype(int)
        boundary_nodes[0:2]=np.array(two_first_nodes)
        for n in range(2,number_of_nodes):
            last_node = boundary_nodes[n-1]
            previous_node = boundary_nodes[n-2]
            # Find all triangles that contains the last known boundary node
            last_node_inds = np.argwhere(tri==last_node)[:,0]
            # Find how many times nodes are found in combination with last known
            # boundary node.
            count = np.unique(tri[last_node_inds,:], return_counts=True)
            unique_nodes = count[0][count[1]==1]
            # Only three should be unique: the previous node, the last node and the next node
            # Remove the previous node
            next_node = np.setdiff1d(unique_nodes, previous_node)
            # Remove the last node if it was unique (depends on the exact triangulation)
            if last_node in unique_nodes:
                next_node = np.setdiff1d(next_node, last_node)

            if next_node.shape != (1,):
                # This should not happen
                raise Exception(f'Could not find a next boundary node after node nr {n} ({last_node})!')
            boundary_nodes[n] = next_node
        return boundary_nodes

    def __init__(self, two_first_nodes: Iterable, number_of_nodes: int):
        self._two_first_nodes = two_first_nodes
        self._number_of_nodes = number_of_nodes

    def __call__(self, nodes, tri, lon, lat):
        bnd_nodes = self.find_boundary_nodes(tri, self._two_first_nodes, self._number_of_nodes)

        other_nodes = np.setdiff1d(np.unique(tri),bnd_nodes)

        # Re-organize nodes so that boundary nodes are first
        new_tri = np.copy(tri)
        new_lon = np.copy(lon)
        new_lat = np.copy(lat)
        for ind, node in enumerate(bnd_nodes):
            new_tri[tri==node] = ind
            new_lon[ind] = lon[node]
            new_lat[ind] = lat[node]
        for ind, node in enumerate(other_nodes):
            new_tri[tri==node] = ind + len(bnd_nodes)
            new_lon[ind + len(bnd_nodes)] = lon[node]
            new_lat[ind + len(bnd_nodes)] = lat[node]

        return range(len(bnd_nodes)), new_tri, range(len(nodes)), new_lon, new_lat

    def __str__(self):
        return f"Setting {self._number_of_nodes} boundary starting from nodes {self._two_first_nodes[0]} and {self._two_first_nodes[1]}"


class ClearBoundary(BoundarySetter):
    def __init__(self):
        return

    def __call__(self, nodes, tri, lon, lat):
        return [], tri, nodes, lon, lat

    def __str__(self):
        return "Clearing all boundary points"


class SetArray(BoundarySetter):
    def __init__(self, boundary_array):
        self._boundary_array = np.array(boundary_array)
        return

    def __call__(self, nodes, tri, lon, lat):
        bnd_points = self._boundary_array
        nodes = np.array(nodes)
        mask = np.logical_and(bnd_points <= max(nodes), bnd_points >= min(nodes))
        return bnd_points[mask], tri, nodes, lon, lat

    def __str__(self):
        return "Setting new boundary points based on an provided array"
