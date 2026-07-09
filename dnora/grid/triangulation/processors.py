import numpy as np

from abc import ABC, abstractmethod
import numpy as np
from typing import Iterable
from .funcs import find_nodes_with_three_neighbours, get_boundary_edges, get_boundary_nodes, force_nodes_as_first_in_triangulation, triangle_angles_xy, outlier_angle_indices_median
from dnora import msg
class TriangulationProcessor(ABC):
    """Rearange triangulation and boundary points."""

    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, nodes, bnd_nodes, tri, lon, lat, flag_points: bool):
        """This method is called from within the Grid-object."""

        return nodes, bnd_nodes, tri, lon, lat

    @abstractmethod
    def __str__(self):
        """Describes how the boundary points are set.

        This is called by the Grid-object to provide output to the user.
        """
        pass




class RemoveNodes(TriangulationProcessor):
    def __init__(self, list_of_nodes: Iterable):
        self._list_of_nodes = list_of_nodes

    def __call__(self, nodes, bnd_nodes, tri, lon, lat, flag_points: bool):
        if flag_points:
            return self._list_of_nodes
        
        msg.plain('\tRemoving nodes...')
        keep = ~np.isin(nodes, self._list_of_nodes)
        nodes, lon, lat = nodes[keep], lon[keep], lat[keep]
        msg.plain('\tRemoving triangles containing the nodes...')
        
        mask_bad = np.isin(tri, self._list_of_nodes).any(axis=1)
        tri = tri[~mask_bad]
        
        
        return nodes, bnd_nodes, tri, lon, lat

    def __str__(self):
        return f"Removing {len(self._list_of_nodes)} nodes..."


# class RemoveTriangle(TriangulationProcessor):
#     def __init__(self, list_of_triangles: Iterable):
#         self._list_of_triangles = list_of_triangles

#     def __call__(self, nodes, bnd_nodes, tri, lon, lat):
#         for ttt in self._list_of_triangles:
#             for n in range(len(tri)):
#                 if ttt[0] in tri[n] and ttt[1] in tri[n] and ttt[2] in tri[n]:
#                     tri = np.delete(tri, n, 0)
#                     break
#         return bnd_nodes, tri, nodes, lon, lat

#     def __str__(self):
#         return f"Removing {len(self._list_of_triangles)} triangles..."


# class ReorganizeBoundary(TriangulationProcessor):
#     @staticmethod
#     def find_edge_nodes(tri: np.ndarray, two_first_nodes: Iterable) -> np.ndarray:
#         """Finds all the consequtive edge nodes when given the seed of the
#         two first ones"""
#         # boundary_nodes = np.zeros(number_of_nodes).astype(int)
#         edge_nodes = np.zeros(len(np.unique(tri))).astype(int)
#         # boundary_nodes[0:2]=np.array(two_first_nodes)
#         edge_nodes[0:2] = np.array(two_first_nodes)
#         first_node = edge_nodes[0]
#         n = 2
#         next_node = -1
#         while next_node != first_node:
#             # for n in range(2,number_of_nodes):

#             last_node = edge_nodes[n - 1]
#             previous_node = edge_nodes[n - 2]
#             # Find all triangles that contains the last known boundary node
#             last_node_inds = np.argwhere(tri == last_node)[:, 0]
#             # Find how many times nodes are found in combination with last known
#             # boundary node.
#             count = np.unique(tri[last_node_inds, :], return_counts=True)
#             unique_nodes = count[0][count[1] == 1]
#             # Only three should be unique: the previous node, the last node and the next node
#             # Remove the previous node
#             next_node = np.setdiff1d(unique_nodes, previous_node)
#             # Remove the last node if it was unique (depends on the exact triangulation)
#             if last_node in unique_nodes:
#                 next_node = np.setdiff1d(next_node, last_node)

#             if next_node.shape != (1,):
#                 # This should not happen
#                 raise Exception(
#                     f"Could not find a next boundary node after node nr {n} ({last_node})!"
#                 )
#             edge_nodes[n] = next_node
#             # if n<number_of_nodes:
#             #     boundary_nodes[n] = next_node
#             n += 1
#         edge_nodes = edge_nodes[: n - 1]
#         return edge_nodes

#     def __init__(self, two_first_nodes: Iterable, number_of_nodes: int):
#         self._two_first_nodes = two_first_nodes
#         self._number_of_nodes = number_of_nodes

#     def __call__(self, nodes, bnd_nodes, tri, lon, lat):
#         edge_nodes = self.find_edge_nodes(tri, self._two_first_nodes)
#         bnd_nodes = edge_nodes[: self._number_of_nodes]
        
#         new_tri, new_lon, new_lat = reorganize_triangulation_with_edges_first(tri, lon, lat, edge_nodes)

#         return range(len(bnd_nodes)), new_tri, range(len(nodes)), new_lon, new_lat

#     def __str__(self):
#         return f"Setting {self._number_of_nodes} boundary starting from nodes {self._two_first_nodes[0]} and {self._two_first_nodes[1]}"


# class ClearBoundary(TriangulationProcessor):
#     def __init__(self):
#         return

#     def __call__(self, nodes, bnd_nodes,tri, lon, lat):
#         return [], tri, nodes, lon, lat

#     def __str__(self):
#         return "Clearing all boundary points"



class ReorderEdges(TriangulationProcessor):
    """Reorganizes the triangulation so that the points that make up the edges (both ocean and land) are consequtive and starting from 0"""
    def __init__(self):
        return

    def __call__(self, nodes, bnd_nodes,tri, lon, lat, flag_points:bool):

        msg.plain('Getting edges of mesh...')
        edges = get_boundary_edges(tri)
        msg.plain('Getting nodes of edges og mesh...')
        edge_nodes = get_boundary_nodes(edges)

        msg.plain("Reordering triangulation and possible boundary nodes...")
        lon, lat, tri, bnd_nodes = force_nodes_as_first_in_triangulation(lon, lat, tri, bnd_nodes, nodes_to_be_first=edge_nodes)

        return nodes, bnd_nodes, tri, lon, lat

    def __str__(self):
        return "Reorganizing edges of grid to be consequtive starting from 0"


class ReorderBoundary(TriangulationProcessor):
    """Reorganizes the triangulation so that the points that make up the boundary poits start from 0"""
    def __init__(self):
        return

    def __call__(self, nodes, bnd_nodes,tri, lon, lat, flag_points: bool):
        msg.plain("Reordering triangulation...")
        lon, lat, tri, bnd_nodes = force_nodes_as_first_in_triangulation(lon, lat, tri, bnd_nodes, nodes_to_be_first=bnd_nodes)

        return nodes, bnd_nodes, tri, lon, lat

    def __str__(self):
        return "Rearganizing boundary points to be consequtive starting from 0"


class Remove3NeighbourTriangles(TriangulationProcessor):
    """Removes triangles containing nodes that have only three neighbours"""
    def __init__(self):
        return

    def __call__(self, nodes, bnd_nodes,tri, lon, lat, flag_points):
        # msg.plain('Getting edges of mesh...')
        # bnd_edges = get_boundary_edges(tri)
        # msg.plain('Getting nodes of edges og mesh...')
        # bnd_nodes = get_boundary_nodes(bnd_edges)
        msg.plain('Finding nodes with three neighbours...')
        bad_nodes = find_nodes_with_three_neighbours(tri)
        
        if flag_points:
            return bad_nodes

        msg.plain('Fixing triangulation...')
        for bad_node in bad_nodes:
            # Two triangles contain this node and needs to be joined
            mask = np.isin(tri, bad_node).any(axis=1)
            tri1 = tri[mask][0,:]
            tri2 = tri[mask][1,:]
            # Modify first triangle by replacing the badd node with the not-common node from the other triangle
            new_node = tri2[~np.isin(tri2, tri1)]
            tri[mask.argmax(),: ] = np.where(tri1 == bad_node, new_node, tri1)
            # Keep old bad node in other triangle thus flagging it for removal in the last step
            


        
        nodes, bnd_nodes, tri, lon, lat = RemoveNodes(bad_nodes)(nodes, bnd_nodes, tri, lon, lat, flag_points=False)

        return nodes, bnd_nodes, tri, lon, lat
        
        #
    def __str__(self):
        return "Removing triangles with only three neighbours"


class RemoveSmallAngleTriangles(TriangulationProcessor):
    """Removes triangles containing a too small angle"""
    def __init__(self, min_angle: float):
        self._min_angle = min_angle

    def __call__(self, nodes, bnd_nodes,tri, lon, lat, flag_points: bool):
        # msg.plain('Getting edges of mesh...')
        # bnd_edges = get_boundary_edges(tri)
        # msg.plain('Getting nodes of edges og mesh...')
        # bnd_nodes = get_boundary_nodes(bnd_edges)
        msg.plain('Calculating angles for triangles...')
        tri_angles = triangle_angles_xy(lon[tri], lat[tri])
        bad_mask = np.any(tri_angles<self._min_angle, axis=1)
        bad_angles = tri_angles[bad_mask]
        msg.plain(f'Identified {sum(bad_mask)} bad triangles...')
        
        
        msg.plain('Determining which node in triangle to remove...')
        cols = outlier_angle_indices_median(bad_angles)
        
        if flag_points:
            tri_to_remove = tri[bad_mask]
            rows = np.arange(tri_to_remove.shape[0])
            bad_nodes = tri_to_remove[rows, cols]
            return bad_nodes
        
        tri = tri[~bad_mask]
        


        # nodes, bnd_nodes, tri, lon, lat = RemoveNodes(bad_nodes)(nodes, bnd_nodes, tri, lon, lat)

        return nodes, bnd_nodes, tri, lon, lat
        
        #
    def __str__(self):
        return f"Removing triangles with an angle smaller than {self._min_angle}"


class EvenOutSpacing(TriangulationProcessor):
    """Removes triangles containing a too small angle"""
    def __init__(self, nodes: Iterable):
        self._nodes = nodes

    def __call__(self, nodes, bnd_nodes,tri, lon, lat, flag_points: bool):
        if flag_points:
            msg.plain(f'Flagging {len(self._nodes)} affected nodes')
            return self._nodes
        
        x, y = lon[self._nodes], lat[self._nodes]

        s = np.r_[0.0, np.cumsum(np.hypot(np.diff(x), np.diff(y)))]
        # Target arc-length positions for N points
        s_target = np.linspace(0.0, s[-1], len(x))

        # Interpolate coordinates at target arc lengths
        lon[self._nodes] = np.interp(s_target, s, x)
        lat[self._nodes] = np.interp(s_target, s, y)

        return nodes, bnd_nodes, tri, lon, lat
        
        #
    def __str__(self):
        return f"Setting given nodes to have an even spacing"
