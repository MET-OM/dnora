import numpy as np

from abc import ABC, abstractmethod
import numpy as np
from typing import Iterable


class TriAranger(ABC):
    """Rearange triangulation and boundary points."""

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


class RemoveTriangle(TriAranger):
    def __init__(self, list_of_triangles: Iterable):
        self._list_of_triangles = list_of_triangles

    def __call__(self, nodes, bnd_nodes, tri, lon, lat):
        for ttt in self._list_of_triangles:
            for n in range(len(tri)):
                if ttt[0] in tri[n] and ttt[1] in tri[n] and ttt[2] in tri[n]:
                    tri = np.delete(tri, n, 0)
                    break
        return bnd_nodes, tri, nodes, lon, lat

    def __str__(self):
        return f"Removing some triangles..."


class ReorganizeBoundary(TriAranger):
    @staticmethod
    def find_edge_nodes(tri: np.ndarray, two_first_nodes: Iterable) -> np.ndarray:
        """Finds all the consequtive edge nodes when given the seed of the
        two first ones"""
        # boundary_nodes = np.zeros(number_of_nodes).astype(int)
        edge_nodes = np.zeros(len(np.unique(tri))).astype(int)
        # boundary_nodes[0:2]=np.array(two_first_nodes)
        edge_nodes[0:2] = np.array(two_first_nodes)
        first_node = edge_nodes[0]
        n = 2
        next_node = -1
        while next_node != first_node:
            # for n in range(2,number_of_nodes):

            last_node = edge_nodes[n - 1]
            previous_node = edge_nodes[n - 2]
            # Find all triangles that contains the last known boundary node
            last_node_inds = np.argwhere(tri == last_node)[:, 0]
            # Find how many times nodes are found in combination with last known
            # boundary node.
            count = np.unique(tri[last_node_inds, :], return_counts=True)
            unique_nodes = count[0][count[1] == 1]
            # Only three should be unique: the previous node, the last node and the next node
            # Remove the previous node
            next_node = np.setdiff1d(unique_nodes, previous_node)
            # Remove the last node if it was unique (depends on the exact triangulation)
            if last_node in unique_nodes:
                next_node = np.setdiff1d(next_node, last_node)

            if next_node.shape != (1,):
                # This should not happen
                raise Exception(
                    f"Could not find a next boundary node after node nr {n} ({last_node})!"
                )
            edge_nodes[n] = next_node
            # if n<number_of_nodes:
            #     boundary_nodes[n] = next_node
            n += 1
        edge_nodes = edge_nodes[: n - 1]
        return edge_nodes

    def __init__(self, two_first_nodes: Iterable, number_of_nodes: int):
        self._two_first_nodes = two_first_nodes
        self._number_of_nodes = number_of_nodes

    def __call__(self, nodes, bnd_nodes, tri, lon, lat):
        edge_nodes = self.find_edge_nodes(tri, self._two_first_nodes)
        bnd_nodes = edge_nodes[: self._number_of_nodes]
        other_nodes = np.setdiff1d(np.unique(tri), edge_nodes)

        # Re-organize nodes so that boundary nodes are first
        new_tri = np.copy(tri)
        new_lon = np.copy(lon)
        new_lat = np.copy(lat)
        for ind, node in enumerate(edge_nodes):
            new_tri[tri == node] = ind
            new_lon[ind] = lon[node]
            new_lat[ind] = lat[node]
        for ind, node in enumerate(other_nodes):
            new_tri[tri == node] = ind + len(edge_nodes)
            new_lon[ind + len(edge_nodes)] = lon[node]
            new_lat[ind + len(edge_nodes)] = lat[node]
        return range(len(bnd_nodes)), new_tri, range(len(nodes)), new_lon, new_lat

    def __str__(self):
        return f"Setting {self._number_of_nodes} boundary starting from nodes {self._two_first_nodes[0]} and {self._two_first_nodes[1]}"


class ClearBoundary(TriAranger):
    def __init__(self):
        return

    def __call__(self, nodes, tri, lon, lat):
        return [], tri, nodes, lon, lat

    def __str__(self):
        return "Clearing all boundary points"


#
#
# class SetArray(TriAranger):
#     def __init__(self, boundary_array):
#         self._boundary_array = np.array(boundary_array)
#         return
#
#     def __call__(self, nodes, tri, lon, lat):
#         bnd_points = self._boundary_array
#         nodes = np.array(nodes)
#         mask = np.logical_and(bnd_points <= max(nodes), bnd_points >= min(nodes))
#         return bnd_points[mask], tri, nodes, lon, lat
#
#     def __str__(self):
#         return "Setting new boundary points based on an provided array"
