from copy import copy
from functools import reduce
from abc import ABC, abstractmethod
from .fvgrid import read_sms_mesh
import meshio

from typing import Iterable
import numpy as np
from dnora.type_manager.data_sources import DataSource
import xarray as xr


class TriangReader(ABC):
    """Abstract class for reading the triangular object."""

    @abstractmethod
    def __call__(self, filename: str) -> tuple:
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

    def default_data_source(self) -> DataSource:
        return DataSource.UNDEFINED


# class NetcdfTriangReader(TriangReader):
#     def __call__(
#         self,
#         source: DataSource,
#         folder: str,
#         filename: str = None,
#         utm: tuple[int, str] = (None, None),
#     ) -> tuple:
#         self.filename = filename

#         ds = xr.open_dataset(get_url(folder, filename))
#         coord_dict = get_coordinates_from_ds(ds, return_dict=True)
#         tri = ds.triangles.values
#         edge_nodes = ds.inds.values[ds.boundary_mask.values.astype(bool)]
#         utm = (None, None)
#         return (
#             tri,
#             coord_dict,
#             edge_nodes,
#             utm[0],
#             utm[1],
#         )

#     def __str__(self):
#         return f"Reading triangular grid from Msh-file {self.filename}."


# class TxtReader(TriangReader):
#     def __init__(self, filename: str, boundary_points=None):
#         self.filename = copy(filename)
#         self.boundary_points = boundary_points
#         if self.boundary_points is None:
#             self.boundary_points = []
#         return

#     def __call__(self) -> tuple:
#         with open(self.filename, "r") as f:
#             nr_of_nodes = int(f.readline())
#             nodes = np.array(range(nr_of_nodes)).astype(int)
#             x = np.full(nodes.shape, 0).astype(float)
#             y = np.full(nodes.shape, 0).astype(float)
#             for n in nodes:
#                 x[n], y[n] = np.array(f.readline().split(" ")).astype(float)
#             nr_of_triangs = int(f.readline())
#             tri = np.full((nr_of_triangs, 3), 0)

#             for n in range(nr_of_triangs):
#                 tri[n, :] = np.fliplr(np.array([f.readline().split(" ")]).astype(int))[
#                     0
#                 ]
#                 # tri[n,:] = np.array([f.readline().split(' ')]).astype(int)[0]

#         # tri, nodes, x, y, Z, types, nodeStrings = read_sms_mesh(self.filename, nodestrings=True)
#         # lat, lon = utm.to_latlon(x, y, 33, zone_letter = 'W', strict = False)
#         nodeStrings = np.array(self.boundary_points)
#         types = None
#         # if nodestring_subset is not None:
#         #    nodeStrings = [nodeStrings[i] for i in nodestring_subset]
#         # nodeStrings = reduce(lambda x, y: x+y, nodeStrings)
#         return tri, nodes, None, None, x, y, types, nodeStrings, 33, "W"

#     def __str__(self):
#         return "Reading triangular grid from Txt-file."


class SmsFileTri(TriangReader):
    def __call__(
        self,
        filename: str,
        utm: tuple[int, str] = (33, "W"),
        nodestring_subset: Iterable = None,
    ) -> tuple:
        """
        Parameters
        ----------
        nodestring_subset : 0-based indices of individual nodestrings to include
            as open boundary nodes, as ordered in the .2dm file.
            Default: include all

        The read_sms_mesh-function is taken directly from the PyFVCOM package
        https://github.com/pwcazenave/pyfvcom"""
        self.filename = filename
        self.umt = utm
        tri, nodes, X, Y, Z, types, nodeStrings = read_sms_mesh(
            self.filename, nodestrings=True
        )
        # lat, lon = utm.to_latlon(X, Y, utm[0], zone_letter=utm[1], strict=False)

        if nodestring_subset is not None:
            nodeStrings = [nodeStrings[i] for i in nodestring_subset]
        nodeStrings = reduce(lambda x, y: x + y, nodeStrings)

        return (
            tri,
            {"x": X, "y": Y},
            nodeStrings,
            utm[0],
            utm[1],
        )

    def __str__(self):
        return f"Reading triangular grid from SMS-file {self.filename} assuming UTM {self.utm}"


class MshReader(TriangReader):
    def __call__(
        self,
        source: DataSource,
        folder: str,
        filename: str = None,
        utm: tuple[int, str] = (None, None),
    ) -> tuple:
        self.filename = filename
        mesh = meshio.read(self.filename)

        for cell in mesh.cells:
            if cell.type == "vertex":  # Boundary points
                nodeStrings = cell.data[:, 0]
            elif cell.type == "triangle":
                tri = cell.data

        x = mesh.points[:, 0]
        y = mesh.points[:, 1]
        # Z = mesh.points[:,2]

        # types = None  # This is not used in DNORA and I have no idea what it is
        # nodes = np.array((range(len(mesh.points[:, 0]))))

        if utm[0] is None:
            coord_dict = {"lon": x, "lat": y}
        else:
            coord_dict = {"x": x, "y": y}

        return (
            tri,
            coord_dict,
            nodeStrings,
            utm[0],
            utm[1],
        )

    def __str__(self):
        return f"Reading triangular grid from Msh-file {self.filename}."
