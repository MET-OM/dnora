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
        nodes = np.array((range(len(mesh.points[:,0]))))+1

        return tri, nodes, lon, lat, types, nodeStrings

    def __str__(self):
        return "Reading triangular grid from Msh-file."
