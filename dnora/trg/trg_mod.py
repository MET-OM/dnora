from copy import copy
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from .read_tr import TriangReader


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
        self._edge_nodes = edge_nodes

    def plot_grid(self) -> None:
        plt.triplot(self.lon(), self.lat(), triangles = self.tri(), linewidth = 0.2, color='black')
        plt.plot(self.lon()[self._edge_nodes],self.lat()[self._edge_nodes],'rx')
        plt.show()
        return


    def name(self):
        return self._name

    def tri(self):
        return self._tri

    def nodes(self):
        return self._nodes

    def lon(self):
        return self._lon

    def lat(self):
        return self._lat

    def edge_nodes(self):
        return self._edge_nodes
