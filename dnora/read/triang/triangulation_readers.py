from copy import copy
from functools import reduce
from .fvgrid import read_sms_mesh
import meshio

from typing import Iterable
import numpy as np
from dnora.type_manager.data_sources import DataSource
import xarray as xr
from dnora.read.abstract_readers import DataReader


class SmsReader(DataReader):
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


class MshReader(DataReader):
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
