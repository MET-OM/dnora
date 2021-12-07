from copy import copy
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
