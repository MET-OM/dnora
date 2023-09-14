from skeletons import GriddedSkeleton, PointSkeleton
import numpy as np
from .. import msg
from skeletons.mask_factory import add_mask
from skeletons.datavar_factory import add_datavar
from .read import TopoReader
from .read import MshFile as topo_MshFile
from .read_tr import TriangReader
from .read_tr import MshFile as triang_MshFile
from .tri_arangers import TriAranger
from copy import copy
from typing import Union
from .. import aux_funcs
from .process import GridProcessor
from .mesh import Mesher, Interpolate


class Topography:
    def import_topo(self, topo_reader: TopoReader) -> None:
        """Reads the raw bathymetrical data."""

        msg.header(topo_reader, "Importing topography...")
        print(topo_reader)
        topo, lon, lat, x, y, zone_number, zone_letter = topo_reader(self.edges('lon'), self.edges('lat'), self.edges('x'), self.edges('y'))

        if 0 in topo.shape:
            msg.warning('Imported topography seems to be empty. Maybe using wrong tile?')
            return
        
        if aux_funcs.is_gridded(topo, lon, lat) or aux_funcs.is_gridded(topo, x, y):
            self._raw = Grid(lon=lon, lat=lat, x=x, y=y)
            x = x or lon
            y = y or lat
            self._raw.set_spacing(nx=len(x), ny=len(y))
        else:
            self._raw = UnstrGrid(lon=lon, lat=lat, x=x, y=y)

        if self.edges('lon', native=True)[0] < self.raw().edges(self.x_str)[0] or self.edges('lon', native=True)[1] > self.raw().edges(self.x_str)[1]:
            msg.warning(f"The data gotten from the TopoReader doesn't cover the grid in the {self.x_str} direction. Grid: {self.edges('lon')}, imported topo: {self.raw().edges(self.x_str)}")

        if self.edges('lat', native=True)[0] < self.raw().edges(self.y_str)[0] or self.edges('lat', native=True)[1] > self.raw().edges(self.y_str)[1]:
            msg.warning(f"The data gotten from the TopoReader doesn't cover the grid in the {self.y_str} direction. Grid: {self.edges('lat')}, imported topo: {self.raw().edges(self.y_str)}")

        if zone_number is not None:
            self._raw.set_utm(zone_number, zone_letter)
        self.raw().set_topo(topo)
        self.raw().set_sea_mask(self.raw().topo()>0)

    def mesh_grid(self, mesher: Mesher=Interpolate(method = 'nearest')) -> None:
        """Meshes the raw data down to the grid definitions."""

        if self.raw() is None:
            msg.warning('Import topography using .import_topo() before meshing!')
            return

        msg.header(mesher, "Meshing grid bathymetry...")
        print(mesher)

        if self.is_gridded():
            xQ, yQ = np.meshgrid(self.x(native=True), self.y(native=True))
        else:
            xQ, yQ = self.xy(native=True)

        if self.is_cartesian():
            x, y = self.raw().xy()
        else:
            x, y = self.raw().lonlat()

        topo = mesher(self.raw().topo().ravel(), x, y, xQ, yQ)

        self.set_topo(topo)
        self.set_sea_mask(self.topo()>0)

    def process_grid(self, grid_processor: GridProcessor=None) -> None:
        """Processes the gridded bathymetrical data, e.g. with a filter."""
        if grid_processor is None:
            return

        msg.header(grid_processor, "Processing meshed grid...")
        print(grid_processor)
        if self.is_gridded():
            topo = grid_processor.grid(self.topo(), self.lon(), self.lat(), self.sea_mask(), self.boundary_mask())
            if topo is None:
                msg.warning('Filtering of gridded topography is not implemented in this GridProcessor.')
                return
        else:
            topo = grid_processor.topo(self.topo(), self.lon(), self.lat(), self.sea_mask())
            if topo is None:
                msg.warning('Filtering of unstructured topography is not implemented in this GridProcessor.')
                return

        self.set_topo(topo)
        self.set_sea_mask(self.topo()>0)

    def set_boundary_points(self, mask_setter) -> None:
        boundary_mask = mask_setter(self)
        self.set_boundary_mask(boundary_mask)

    def time(self) -> tuple:
        return (None, None)

    def tri(self):
        if hasattr(self, '_tri'):
            return copy(self._tri)
        else:
            return None

    def raw(self):
        if hasattr(self, '_raw'):
            return self._raw
        return None
    
    def cfl(self, dx=None, f0=0.041180):
        """Calculates approximate time step [s].
        Based on grid resolution and given lowest frequency [Hz] (default=0.041180)
        """
        if dx is None:
            dx = min(self.dx(), self.dy()) # Grid spacing [m]

        cg = 1.56/f0*0.5 # Deep water group velocity [m/s]
        dt = dx/cg

        print(f'Grid spacing dx = {dx:.0f} m and f0 = {f0:.8f} Hz')
        print(f'Approximate minimum time step: dt = dx/cg = {dx:.0f}/{cg:.1f} = {dt:.1f} s')

        return dt

@add_datavar(name='topo', default_value=999.)
@add_mask(name='boundary', coords='grid', default_value=0)
@add_mask(name='sea', coords='grid', default_value=1, opposite_name='land')
class Grid(GriddedSkeleton, Topography):
    def boundary_nx(self) -> int:
        """Return approximate number of grid points in the longitude direction
        """
        abs_diff = np.abs(np.diff(np.where(self.boundary_mask())))
        if abs_diff.size == 0:
            return 0
        abs_diff=np.median(abs_diff[abs_diff>0]).astype(int)

        return np.ceil(self.nx()/abs_diff).astype(int)

    def boundary_ny(self) -> int:
        """Return approximate number of grid points in the longitude direction
        """
        abs_diff = np.abs(np.diff(np.where(self.boundary_mask())))
        if abs_diff.size == 0:
            return 0
        abs_diff=np.median(abs_diff[abs_diff>0]).astype(int)

        return np.ceil(self.ny()/abs_diff).astype(int)


@add_datavar(name='topo', default_value=999.)
@add_mask(name='boundary', coords='grid', default_value=0)
@add_mask(name='sea', coords='grid', default_value=1, opposite_name='land')
class UnstrGrid(PointSkeleton, Topography):
    pass

class TriGrid(UnstrGrid):
    @classmethod
    def from_msh(cls, filename: str, name: str='AnonymousGrid'):
        tri_grid = cls(name=name)
        tri_grid.import_triang(triang_MshFile(filename))
        tri_grid.import_topo(topo_MshFile(filename))
       
        return tri_grid

    def __init__(self, x=None, y=None, lon=None, lat=None, name='AnonymousGrid'):
        self.name = name
        # Only initialize if x, y, lon, lat given
        if [a for a in (x, y, lon, lat) if a is not None]:
            self._init_structure(x, y, lon, lat, name)

    def import_triang(self, triang_reader: TriangReader):
        """Reads a triangular mesh."""
        tri, nodes, lon, lat, x, y, types, edge_nodes, zone_number, zone_letter = triang_reader()

        self._init_structure(x, y, lon, lat)

        self.set_utm(zone_number, zone_letter)
        edge_nodes = np.array(edge_nodes)
        edge_nodes = edge_nodes.astype(int)
        self._update_boundary(edge_nodes)
        self._tri = tri
        #self._nodes = nodes # These are now in self.inds()
        self._types = types #???

    def arange_triangulation(self, tri_aranger: TriAranger) -> None:
        print(tri_aranger)
        bnd_nodes, tri, nodes, x, y = tri_aranger(self.inds(),  np.where(self.boundary_mask())[0], self.tri(), self.x(native=True), self.y(native=True))

        x, y = self.xy(strict=True)
        lon, lat = self.lonlat(strict=True)
        self._init_structure(x=x, y=y, lon=lon, lat=lat)
        
        self._update_boundary(bnd_nodes)
        self._tri = tri

    def _update_boundary(self, boundary_inds):
        mask = np.array([ind in boundary_inds for ind in self.inds()])
        self.set_boundary_mask(mask)



# if self.x_str == 'x':
#             self._init_structure(x=x, y=y, lon=None, lat=None)
#         else:
#             self._init_structure(x=None, y=None, lon=x, lat=y)