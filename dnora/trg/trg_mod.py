import xarray as xr
from copy import copy
import numpy as np
from typing import Iterable
from ..grd import TopoReader
from .read_tr import TriangReader
from .write import TrGridWriter
from .boundary import BoundarySetter, ClearBoundary
from .plot import TrGridPlotter, TriTopoPlotter
from ..grd.mesh import Mesher, Interpolate
from ..grd.process import GridProcessor
from ..grd.grd_mod import force_to_xyz
from .. import msg
from .. import file_module
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

        topo, lon, lat = force_to_xyz(topo, lon, lat)

        # Depth is positive, so set everything that is not positive to nan
        topo[topo<=0]=np.nan

        # This was used for structured topography
        #coords_dict = {'lon': lon, 'lat': lat}
        #vars_dict = {'topo': (['lat', 'lon'], topo)}

        points = [x for x in range(len(lon))]
        coords_dict = {'points': points}
        vars_dict = {'topo': (['points'], topo), 'lon': (['points'], lon), 'lat': (['points'], lat)}
        self.rawdata = xr.Dataset(
                    coords=(coords_dict
                    ),
                    data_vars=(vars_dict
                    ),
                    )


        return

    def append_boundary(self, boundary_setter: BoundarySetter) -> None:
        """Set new boundary points but keep the old ones."""

        new_points, new_tri, new_nodes, new_lon, new_lat = boundary_setter(self.nodes(), self.tri(), self.lon(), self.lat())

        self._boundary = np.union1d(np.array(new_points).astype(int), self.boundary_inds())
        self._tri = np.array(new_tri).astype(int)
        self._nodes = np.array(new_nodes).astype(int)
        self._lon = np.array(new_lon)
        self._lat = np.array(new_lat)

    def set_boundary(self, boundary_setter: BoundarySetter) -> None:
        """Set boundary points. Old ones not kept."""

        new_points, new_tri, new_nodes, new_lon, new_lat = boundary_setter(self.nodes(), self.tri(), self.lon(), self.lat())
        self._boundary = np.array(new_points).astype(int)
        self._tri = np.array(new_tri).astype(int)
        self._nodes = np.array(new_nodes).astype(int)
        self._lon = np.array(new_lon)
        self._lat = np.array(new_lat)

    def mesh_grid(self, mesher: Mesher=Interpolate(method = 'linear')) -> None:
        """Meshes the raw data down to the grid definitions."""

        if self.tri() is not None:

            msg.header(mesher, "Meshing grid bathymetry...")
            print(mesher)
            topo = mesher(self.raw_topo(), self.raw_lon(), self.raw_lat(), self.lon(), self.lat())
            topo[topo < 0] = 0
            self.data = topo
            coords_dict = {'nodes': self.nodes()}
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

    def process_grid(self, filt: GridProcessor=None) -> None:
        """Processes the gridded bathymetrical data, e.g. with a filter."""

        if filt is None:
            msg.warning('Provide a GridProcessor!')
            return

        msg.header(filt, "Processing meshed grid...")
        print(filt)
        true_mask = np.full(self.topo().shape, True)
        topo = filt.grid(self.topo(), self.lon(), self.lat(),
                        sea_mask = np.full(self.topo().shape, True),
                        boundary_mask = np.full(self.topo().shape, False))

        if topo is None:
            msg.warning('Processing of mesed topography is not implemented')
        else:
            vars_dict = {'topo': ('nodes', topo)}
            self.data = self.data.assign(vars_dict)


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
            return np.array([])

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
            return np.array([])

    def ny(self) -> int:
        """Return the number of points in longitude direction."""
        return len(self.lat())


    def nx(self) -> int:
        """Return the number of points in latitude direction."""
        return 1

    def size(self) -> tuple:
        """Returns the size (nx, ny) of the grid."""
        #return self.land_sea_mask().shape
        return (self.ny(), self.nx())


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

        filename = file_module.add_folder_to_filename(filename, folder)
        msg.to_file(filename)

        stdout = sys.stdout
        sys.stdout = open(filename, 'w')
        print(self)
        sys.stdout.close()
        sys.stdout = stdout

        return

    def lon_edges(self) -> tuple[float, float]:
        return (np.min(self.lon()), np.max(self.lon()))

    def lat_edges(self) -> tuple[float, float]:
        return (np.min(self.lat()), np.max(self.lat()))

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

    def __repr__(self):
        lines = [f"<dnora Grid object> (unstructured)", f"  Name: {self.name()}"]

        if self.topo().shape==(0,):
            empty_topo = True
        elif np.mean(self.topo()) == 9999:
            empty_topo = True
        else:
            empty_topo = False

        if self.raw_topo().shape==(0,):
            empty_raw_topo = True
        elif np.mean(self.raw_topo()) == 9999:
            empty_raw_topo = True
        else:
            empty_raw_topo = False

        if len(self.nodes())>0:
            lines.append(f"  Number of points: {self.nodes().shape[0]}")
            lines.append(f"  Number of triangles: {self.tri().shape[0]}")
        else:
            lines.append(f"  Number of points: Use method .import_triang() to set structure.")
        lines.append(f"  Data:")
        if not empty_raw_topo:
            lines.append(f'\traw_topo: {self.raw_topo().shape}')
        else:
            lines.append(f'\traw_topo: import using .import_topo()')
        if not empty_topo:
            lines.append(f'\ttopo {self.topo().shape}')
        else:
            lines.append(f'\ttopo: mesh using .mesh_grid()')

        lines.append('\n  Use print() for grid details.')

        return "\n".join(lines)
