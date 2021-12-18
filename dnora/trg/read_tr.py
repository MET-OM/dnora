from copy import copy
from abc import ABC, abstractmethod
from .fvgrid import read_sms_mesh
#import utm
from typing import Tuple

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

    def __call__(self) -> Tuple:
        """The read_sms_mesh-function is taken directly from the PyFVCOM package
        https://github.com/pwcazenave/pyfvcom"""

        import utm

        tri, nodes, X, Y, Z, types, nodeStrings = read_sms_mesh(self.filename, nodestrings=True)

        lat, lon = utm.to_latlon(X, Y, 33, zone_letter = 'W', strict = False)
        nodeStrings=nodeStrings[0]
        return tri, nodes, lon, lat, types, nodeStrings

    def __str__(self):
        return "Reading triangular grid from SMS-file."
