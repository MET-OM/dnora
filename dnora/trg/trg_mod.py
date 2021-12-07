from copy import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from .read_tr import TriangReader
from .boundary import BoundarySetter, ClearBoundary


class TrGrid:
    def __init__(self, name='AnonymousTrianGrid'):
        self._name = copy(name)
        return

    def import_triang(self, triang_reader: TriangReader):
        """Reads a triangular mesh."""
        tri, nodes, lon, lat, types, edge_nodes = triang_reader()

        self._tri = tri
        self._nodes = nodes
        self._lon = lon
        self._lat = lat
        self._types = types #???
        edge_nodes = np.array(edge_nodes)
        edge_nodes = edge_nodes.astype(int)
        self._boundary = edge_nodes

    def plot_grid(self) -> None:
        plt.triplot(self.lon(), self.lat(), triangles = self.tri(), linewidth = 0.2, color='black')
        plt.plot(self.lon()[self._edge_nodes],self.lat()[self._edge_nodes],'rx')
        plt.show()
        return

    def append_boundary(self, boundary_setter: BoundarySetter) -> None:
        """Set new boundary points but keep the old ones."""

        new_points = boundary_setter(self.nodes())
        new_points = np.array(new_points)
        new_points = new_points.astype(int)
        self._boundary = np.union1d(new_points, self.boundary())

        return

    def set_boundary(self, boundary_setter: BoundarySetter) -> None:
        """Set boundary points. Old ones not kept."""

        new_points = boundary_setter(self.nodes())
        new_points = np.array(new_points)
        new_points = new_points.astype(int)
        self._boundary = new_points

        return

    def name(self):
        return copy(self._name)

    def tri(self):
        return copy(self._tri)

    def nodes(self):
        return copy(self._nodes)

    def lon(self):
        return copy(self._lon)

    def lat(self):
        return copy(self._lat)

    def boundary(self):
        return copy(self._boundary)
