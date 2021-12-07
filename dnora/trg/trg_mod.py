
from copy import copy
import numpy as np

from .read_tr import TriangReader
from .boundary import BoundarySetter, ClearBoundary
from .plot import TrGridPlotter, TriPlotter

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

    def plot_grid(self, grid_plotter: TrGridPlotter=None) -> None:
        if grid_plotter is None:
            self._grid_plotter = self._get_grid_plotter()
        else:
            self._grid_plotter = grid_plotter

        if self._grid_plotter is None:
            raise Exception('Define a TrGridPlotter!')

        fig, filename = self._grid_plotter(self)
        fig.show()

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

    def _get_grid_plotter(self) -> TrGridPlotter:
        return TriPlotter()
