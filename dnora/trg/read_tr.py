from copy import copy
from functools import reduce
from abc import ABC, abstractmethod
from .fvgrid import read_sms_mesh
#import utm
from typing import Tuple, Iterable
import numpy as np

class TriangReader(ABC):
    """Abstract class for reading the triangular object."""

    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, filename: str) -> Tuple:
        """Reads the triangular grid.

        This method is called from within the TrGrid-object
        """
        pass

    @abstractmethod
    def __str__(self):
        """Describes what triangular grid is read and from wher.

        This is called by the TrGrid-object to provide output to the user.
        """
        pass

class TxtReader(TriangReader):
    def __init__(self, filename: str, boundary_points=None):
        self.filename = copy(filename)
        self.boundary_points = boundary_points
        if self.boundary_points is None:
            self.boundary_points = []
        return

    def __call__(self) -> Tuple:

        import utm
        with open(self.filename, 'r') as f:
            nr_of_nodes=int(f.readline())
            nodes = np.array(range(nr_of_nodes)).astype(int)
            x = np.full(nodes.shape,0)
            y = np.full(nodes.shape,0)
            for n in nodes:
                x[n], y[n] = np.array(f.readline().split(' ')).astype(float)

            nr_of_triangs=int(f.readline())
            tri = np.full((nr_of_triangs, 3),0)

            for n in range(nr_of_triangs):
                tri[n,:] = np.fliplr(np.array([f.readline().split(' ')]).astype(int))[0]
                #tri[n,:] = np.array([f.readline().split(' ')]).astype(int)[0]

        #tri, nodes, x, y, Z, types, nodeStrings = read_sms_mesh(self.filename, nodestrings=True)
        lat, lon = utm.to_latlon(x, y, 33, zone_letter = 'W', strict = False)
        nodeStrings = np.array(self.boundary_points)
        types = None
        #if nodestring_subset is not None:
        #    nodeStrings = [nodeStrings[i] for i in nodestring_subset]
        #nodeStrings = reduce(lambda x, y: x+y, nodeStrings)
        return tri, nodes, lon, lat, types, nodeStrings

    def __str__(self):
        return "Reading triangular grid from SMS-file."

class SmsReader(TriangReader):
    def __init__(self, filename: str):
        self.filename = copy(filename)
        return

    def __call__(self, nodestring_subset: Iterable=None) -> Tuple:
        """
        Parameters
        ----------
        nodestring_subset : 0-based indices of individual nodestrings to include
            as open boundary nodes, as ordered in the .2dm file.
            Default: include all

        The read_sms_mesh-function is taken directly from the PyFVCOM package
        https://github.com/pwcazenave/pyfvcom"""

        import utm

        tri, nodes, X, Y, Z, types, nodeStrings = read_sms_mesh(self.filename, nodestrings=True)
        lat, lon = utm.to_latlon(X, Y, 33, zone_letter = 'W', strict = False)

        if nodestring_subset is not None:
            nodeStrings = [nodeStrings[i] for i in nodestring_subset]
        nodeStrings = reduce(lambda x, y: x+y, nodeStrings)
        return tri, nodes, lon, lat, types, nodeStrings

    def __str__(self):
        return "Reading triangular grid from SMS-file."

class MshReader(TriangReader):
    def __init__(self, filename: str):
        self.filename = copy(filename)
        return

    def __call__(self) -> Tuple:

        import meshio
        mesh = meshio.read(self.filename)

        for cell in mesh.cells:
            if cell.type == 'vertex': # Boundary points
                nodeStrings = cell.data[:,0]
            elif cell.type == 'triangle':
                tri = cell.data

        lon = mesh.points[:,0]
        lat = mesh.points[:,1]
        #Z = mesh.points[:,2]

        types = None # This is not used in DNORA and I have no idea what it is
        nodes = np.array((range(len(mesh.points[:,0]))))

        return tri, nodes, lon, lat, types, nodeStrings

    def __str__(self):
        return "Reading triangular grid from Msh-file."
