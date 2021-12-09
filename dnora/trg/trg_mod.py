import xarray as xr
from copy import copy
import numpy as np

from ..grd import TopoReader
from .read_tr import TriangReader
from .write import TrGridWriter
from .boundary import BoundarySetter, ClearBoundary
from .plot import TrGridPlotter, TriTopoPlotter
from ..grd.mesh import Mesher, Interpolate
from .. import msg
from ..aux import create_filename_obj, check_if_folder, add_folder_to_filename
from ..defaults import dflt_grd
import sys

class Grid:
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

    def append_boundary(self, boundary_setter: BoundarySetter) -> None:
        """Set new boundary points but keep the old ones."""

        new_points = boundary_setter(self.nodes())
        new_points = np.array(new_points)
        new_points = new_points.astype(int)
        self._boundary = np.union1d(new_points, self.boundary_inds())

        return

    def set_boundary(self, boundary_setter: BoundarySetter) -> None:
        """Set boundary points. Old ones not kept."""

        new_points = boundary_setter(self.nodes())
        new_points = np.array(new_points)
        new_points = new_points.astype(int)
        self._boundary = new_points

        return

    def mesh_grid(self, mesher: Mesher=Interpolate(method = 'linear')) -> None:
        """Meshes the raw data down to the grid definitions."""

        if self.tri() is not None:

            msg.header(mesher, "Meshing grid bathymetry...")
            print(mesher)
            topo = mesher(self.raw_topo(), self.raw_lon(), self.raw_lat(), self.lon(), self.lat())
            topo[topo < 0] = 0
            self.data = topo
            coords_dict = {'nodes': self.nodes()}
            #coords_dict = {'lon': self.lon(), 'lat': self.lat()}
            vars_dict = {'lon': ('nodes', self.lon()), 'lat': ('nodes', self.lat()), 'topo': ('nodes', topo)}
            self.data = xr.Dataset(
                     coords=(coords_dict
                     ),
                     data_vars=(vars_dict
                     ),
                     )
            return

        else:
            msg.plain('No triangular grid imported!')
            return

    def export_grid(self, grid_writer: TrGridWriter, out_format: str=None, filestring: str=None, infofilestring: str=None, folder: str=None) -> None:
        """Exports the boundary spectra to a file.

        The grid_writer defines the file format.
        """

        msg.header(grid_writer, f"Writing grid topography from {self.name()}")

        ### Formats of file names etc.
        out_format = out_format or grid_writer._preferred_format()

        filestring = filestring or dflt_grd['fs'][out_format]
        infofilestring = infofilestring or dflt_grd['info'][out_format]
        folderstring = folder or dflt_grd['fldr'][out_format]

        ### File names
        filename = create_filename_obj(filestring=filestring, objects=[self])
        infofilename = create_filename_obj(filestring=infofilestring, objects=[self])
        folder = create_filename_obj(filestring=folderstring, objects=[self])

        existed = check_if_folder(folder=folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {folder}")

        output_files, output_folder = grid_writer(self, filename=filename, infofilename=infofilename, folder=folder)

        return output_files, output_folder

    def plot_grid(self, grid_plotter: TrGridPlotter=None) -> None:
        self._grid_plotter = grid_plotter or self._get_grid_plotter()

        if self._grid_plotter is None:
            raise Exception('Define a TrGridPlotter!')

        fig, filename = self._grid_plotter(self)
        fig.show()

        return


    def name(self):
        if hasattr(self, '_name'):
            return copy(self._name)
        else:
            return ''

    def structured(self):
        return False

    def tri(self):
        if hasattr(self, '_tri'):
            return copy(self._tri)
        else:
            return None

    def raw_topo(self):
        """Returns an array containing the unmeshed imported topography."""
        if hasattr(self, 'rawdata') and hasattr(self.rawdata, 'topo'):
            return copy(self.rawdata.topo.values)
        else:
            return None

    def raw_lon(self):
        """Returns a longitude vector of the unmeshed imported topography."""
        if hasattr(self, 'rawdata') and hasattr(self.rawdata, 'lon'):
            return copy(self.rawdata.lon.values)
        else:
            return None

    def raw_lat(self):
        """Returns a latitude vector of the unmeshed imported topography."""
        if hasattr(self, 'rawdata') and hasattr(self.rawdata, 'lat'):
            return copy(self.rawdata.lat.values)
        else:
            return None

    def nodes(self):
        if hasattr(self, '_nodes'):
            return copy(self._nodes)
        else:
            return None

    def topo(self):
        """Returns an array containing the meshed topography."""
        if hasattr(self, 'data') and hasattr(self.data, 'topo'):
            return copy(self.data.topo.values)
        else:
            return np.array([])

    def lon(self):
        """Return longitude vector."""
        if hasattr(self, '_lon'):
            return copy(self._lon)
        else:
            return None

    def lat(self):
        """Return latitude vector."""
        if hasattr(self, '_lat'):
            return copy(self._lat)
        else:
            return None

    def boundary_inds(self):
        if hasattr(self, '_boundary'):
            return copy(self._boundary)
        else:
            return np.array([])

    def boundary_points(self):
        mask = self.boundary_inds()
        return np.transpose(np.array([self.lon()[mask], self.lat()[mask]]))

    def write_status(self, filename='', folder='') -> None:
        """Writes out the status of the grid to a file."""

        if not filename:
            filename = f"{self.name()}_info.txt"

        filename = add_folder_to_filename(filename, folder)
        msg.to_file(filename)

        stdout = sys.stdout
        sys.stdout = open(filename, 'w')
        print(self)
        sys.stdout.close()
        sys.stdout = stdout

        return

    def __str__(self) -> str:
        """Prints status of the grid."""

        empty_topo = np.mean(self.topo()) == 9999
        #empty_topo = False # Testing
        msg.header(self, f"Status of grid {self.name()}")
        if self.lon() is not None:
            msg.plain(f'lon: {min(self.lon())} - {max(self.lon())}, lat: {min(self.lat())} - {max(self.lat())}')

        if self.topo().size > 0 and not empty_topo:
            msg.plain(f"Mean depth: {np.mean(self.topo()):.1f} m")
            msg.plain(f"Max depth: {np.max(self.topo()):.1f} m")
            msg.plain(f"Min depth: {np.min(self.topo()):.1f} m")

        if self.nodes() is not None:
            msg.print_line()
            msg.plain('Grid contains:')
            msg.plain(f'{len(self.nodes())} nodes')
            msg.plain(f'{self.tri().shape[0]} triangle objects.')
        if len(self.boundary_inds())>0:
            msg.plain(f'{len(self.boundary_inds()):d} boundary points')

        msg.print_line()

        return ''

    def _get_grid_plotter(self) -> TrGridPlotter:
        return TriTopoPlotter()
