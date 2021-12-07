import xarray as xr
from copy import copy
import numpy as np

from ..grd import TopoReader
from .read_tr import TriangReader
from .boundary import BoundarySetter, ClearBoundary
from .plot import TrGridPlotter, TriTopoPlotter
from ..grd.mesh import Mesher, Interpolate
from .. import msg

class TrGrid:
    def __init__(self, name='AnonymousTriangGrid'):
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

    def import_topo(self, topo_reader: TopoReader) -> None:
        """Reads the raw bathymetrical data."""

        msg.header(topo_reader, "Importing topography...")
        print(topo_reader)
        topo, lon, lat = topo_reader(min(self.lon()), max(self.lon()), min(self.lat()), max(self.lat()))
        topo[topo < 0] = 0
        coords_dict = {'lon': lon, 'lat': lat}
        vars_dict = {'topo': (['lat', 'lon'], topo)}
        self.rawdata = xr.Dataset(
                    coords=(coords_dict
                    ),
                    data_vars=(vars_dict
                    ),
                    )
        return

    def mesh_grid(self, mesher: Mesher=Interpolate(method = 'linear')) -> None:
        """Meshes the raw data down to the grid definitions."""

        if self.tri() is not None:

            msg.header(mesher, "Meshing grid bathymetry...")
            print(mesher)
            topo = mesher(self.raw_topo(), self.raw_lon(), self.raw_lat(), self.lon(), self.lat())
            topo[topo < 0] = 0
            self.data = topo
            # coords_dict = {'lon': self.lon(), 'lat': self.lat()}
            # vars_dict = {'topo': (['lat', 'lon'], topo)}
            # self.data = xr.Dataset(
            #             coords=(coords_dict
            #             ),
            #             data_vars=(vars_dict
            #             ),
            #             )
            return

        else:
            msg.plain('No triangular grid imported!')
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


    def name(self):
        return copy(self._name)

    def tri(self):
        if hasattr(self, '_tri'):
            return copy(self._tri)
        else:
            return None

    def raw_topo(self):
        """Returns an array containing the unmeshed imported topography."""
        return copy(self.rawdata.topo.values)

    def raw_lon(self):
        """Returns a longitude vector of the unmeshed imported topography."""
        return copy(self.rawdata.lon.values)

    def raw_lat(self):
        """Returns a latitude vector of the unmeshed imported topography."""
        return copy(self.rawdata.lat.values)

    def nodes(self):
        return copy(self._nodes)

    def topo(self):
        """Returns an array containing the unmeshed imported topography."""
        return copy(self.data)

    def lon(self):
        return copy(self._lon)

    def lat(self):
        return copy(self._lat)

    def boundary(self):
        return copy(self._boundary)

    def _get_grid_plotter(self) -> TrGridPlotter:
        return TriTopoPlotter()
