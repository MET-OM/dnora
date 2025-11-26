import numpy as np

from abc import ABC, abstractmethod
import numpy as np
from typing import Iterable

from dnora import msg
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
def reorganize_triangulation_with_edges_first(tri, lon, lat, edge_nodes):
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
    return new_tri, new_lon, new_lat

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
        
        new_tri, new_lon, new_lat = reorganize_triangulation_with_edges_first(tri, lon, lat, edge_nodes)

        return range(len(bnd_nodes)), new_tri, range(len(nodes)), new_lon, new_lat

    def __str__(self):
        return f"Setting {self._number_of_nodes} boundary starting from nodes {self._two_first_nodes[0]} and {self._two_first_nodes[1]}"


class ClearBoundary(TriAranger):
    def __init__(self):
        return

    def __call__(self, nodes, bnd_nodes,tri, lon, lat):
        return [], tri, nodes, lon, lat

    def __str__(self):
        return "Clearing all boundary points"

def reorder_triangulation_index(lon, lat, tri, old_inds):
    """Reorders the index values in old_inds to new_inds and keeps consistens lon, lat and trianulation values"""
    
    new_lon, new_lat = np.zeros(len(lon)), np.zeros(len(lat))-1
    new_tri = np.zeros(tri.shape).astype(int)
    for new, old in enumerate(old_inds):
        new_lon[new] = lon[old]
        new_lat[new] = lat[old]
        mask = tri == old
        new_tri[mask] = new
    all_old_inds = range(len(lon))
    other_old_inds = list(set(all_old_inds) - set(old_inds))

    for new, old in enumerate(other_old_inds):
        new_lon[new+len(old_inds)] = lon[old]
        new_lat[new+len(old_inds)] = lat[old]
        mask = tri == old

        new_tri[mask] = new+len(old_inds)

    return new_lon, new_lat, new_tri
def get_boundary_edges(tri):
    """Finds boundary edges in triangulation"""
    # Extract all edges from triangles
    edges = np.vstack([
        tri[:, [0, 1]],  # First edge of each triangle
        tri[:, [1, 2]],  # Second edge of each triangle
        tri[:, [2, 0]]   # Third edge of each triangle
    ])

    # Sort the edges to ensure consistent ordering (e.g., [1, 2] is the same as [2, 1])
    edges = np.sort(edges, axis=1)

    # Count the occurrences of each edge
    from collections import Counter

    edge_counts = Counter(map(tuple, edges))

    # Boundary edges are those that occur only once
    boundary_edges = [edge for edge, count in edge_counts.items() if count == 1]
    return boundary_edges

def get_boundary_nodes(edges):
    """Gets all the nodes that form the edges in a consequtive order"""
    ind = 0
    shorelines = [] # Shorelines of boundary and all islands
    bnd_nodes = []
    edge_array = np.array(edges)
    n0, n1 = edges[ind]
    bnd_nodes.append(n0)
    bnd_nodes.append(n1)
    edge_array[ind,:] = -1
    next_node = n1
    # print(f"Node {n0} connects to {n1}")
    # print(f"Next node: {next_node}")
    while np.max(edge_array) > -1:
        try:
            ind = np.argwhere(edge_array == next_node)[:,0][0]
        except IndexError:
            shorelines.append(bnd_nodes)
            bnd_nodes = []
            ind = np.where(np.max(edge_array,axis=1)>-1)[0][0]
            
        n0, n1 = edges[ind]
        
        if not bnd_nodes or n0 == bnd_nodes[-1]:
            next_node = n1
        else:
            next_node = n0
        #print(f"Node {n0} connects to {n1}")
        #print(f"Next node: {next_node}")
        bnd_nodes.append(next_node)
        edge_array[ind,:] = -1

    # Assume that the longest consequitve shoreline is boundary/land and rest are islands
    shore_ind = np.argmax([len(bnd) for bnd in shorelines])
    return shorelines[shore_ind][:-1]
class ReorderEdges(TriAranger):
    """Reorganizes the triangulation so that the points that make up the edges (both ocean and land) are consequtive and starting from 0"""
    def __init__(self):
        return

    def __call__(self, nodes, bnd_nodes,tri, lon, lat):

        msg.plain('Getting edges of mesh...')
        bnd_edges = get_boundary_edges(tri)
        msg.plain('Getting nodes of edges og mesh...')
        bnd_nodes = get_boundary_nodes(bnd_edges)

        msg.plain("Reordering triangulation...")
        new_lon, new_lat, new_tri = reorder_triangulation_index(lon, lat, tri, bnd_nodes)
        return [], new_tri, nodes, new_lon, new_lat

    def __str__(self):
        return "Rearganizing boundary to be consequtive starting from 0"


class ReorderBoundary(TriAranger):
    """Reorganizes the triangulation so that the points that make up the boundary poits start from 0"""
    def __init__(self):
        return

    def __call__(self, nodes, bnd_nodes,tri, lon, lat):
        msg.plain("Reordering triangulation...")
        new_lon, new_lat, new_tri = reorder_triangulation_index(lon, lat, tri, bnd_nodes)
        return list(range(len(bnd_nodes))), new_tri, nodes, new_lon, new_lat

    def __str__(self):
        return "Rearganizing boundary to be consequtive starting from 0"

def find_nodes_with_three_neighbours(tri) -> list[int]:
    """Finds nodes that have exacylt three neighbours"""
    nodes = range(np.max(tri))
    bad_nodes = []
    for n in nodes:
        mask = tri == n
        num_of_neighbours = np.sum(np.sum(mask, axis=1)>0)
        if num_of_neighbours == 3:
            bad_nodes.append(n)
    return bad_nodes
class Remove3NeighbourEdgeTriangles(TriAranger):
    """Reorganizes the triangulation so that the points that make up the boundary poits start from 0"""
    def __init__(self):
        return

    def __call__(self, nodes, bnd_nodes,tri, lon, lat):
        msg.plain('Getting edges of mesh...')
        bnd_edges = get_boundary_edges(tri)
        msg.plain('Getting nodes of edges og mesh...')
        bnd_nodes = get_boundary_nodes(bnd_edges)
        bad_nodes = find_nodes_with_three_neighbours(tri)
        bad_nodes = list(set.intersect(set(bnd_nodes), set(bad_nodes)))
        lon[bad_nodes] = []
        lat[bad_nodes] = []
        for n in range(tri.shape[0]):
            if tri[n,0] in bad_nodes or tri[n,1] in bad_nodes or tri[n,2] in bad_nodes:
                tri[n,:] = []
        new_bnd_nodes = list(set(bnd_nodes)- set(bad_nodes))
        new_lon, new_lat, new_tri = reorder_triangulation_index(lon, lat, tri, new_bnd_nodes)
        return list(range(len(bnd_nodes))), new_tri, nodes, new_lon, new_lat
        
        #
    def __str__(self):
        return "Rearganizing boundary to be consequtive starting from 0"

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
