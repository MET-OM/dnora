from ..skeletons.gridded_skeleton import GriddedSkeleton
from ..skeletons.point_skeleton import PointSkeleton
from ..skeletons.topography import topography_methods
import numpy as np
import xarray as xr
from .. import msg
from ..skeletons.coordinate_factory import add_time
from ..skeletons.mask_factory import add_mask
from ..skeletons.datavar_factory import add_datavar
from .. import aux_funcs
from .read import TopoReader
from .read import MshFile as topo_MshFile
from .read_tr import TriangReader
from .read_tr import MshFile as triang_MshFile
from .tri_arangers import TriAranger
from copy import copy
from typing import Union

@topography_methods
@add_datavar(name='topo', default_value=999., stash_get=True)
@add_mask(name='waveseries', coords='grid', default_value=0)
@add_mask(name='spectra', coords='grid', default_value=0)
@add_mask(name='boundary', coords='grid', default_value=0)
@add_mask(name='sea', coords='grid', default_value=1, opposite_name='land')
class Grid(GriddedSkeleton):
    def __init__(self, x=None, y=None, lon=None, lat=None, name='AnonymousGrid'):
        self.name = name
        self._init_structure(x, y, lon, lat)

    def set_spacing(self, dlon: float=0., dlat: float=0.,
                    dx: float=0., dy: float=0., dm: float=0,
                    nx: int=0, ny: int=0, floating_edge: bool=False) -> None:
        """Defines longitude and latitude vectors based on desired spacing.

        Options (priority in this order)
        nx, ny [grid pnts]: Grid resolution is set to have nx points in
                            longitude and ny points in latitude direction.

        dlon, dlat [deg]:   Grid spacing is set as close to the given resolution
                            as possible (edges are fixed).

        dm [m]:             Grid spacing is set as close as dm metres as
                            possible.

        dx, dy [m]:         Grid spacing is set as close as dx and xy metres as
                            possible.

        Set floating_edge=True to force exact dlon, dlat
        and instead possibly move lon_max, lat_max slightly
        to make it work (only compatibel with native coordinates).

        """
        def determine_nx_ny(nx, ny, dx, dy, dlon, dlat):
            x_end = self.edges('x', native=True)[1]
            y_end = self.edges('y', native=True)[1]

            if nx and ny:
                return int(nx), int(ny), x_end, y_end

            if dlon and dlat:
                nx = np.round((self.edges('lon')[1]-self.edges('lon')[0])/dlon) + 1
                ny = np.round((self.edges('lat')[1]-self.edges('lat')[0])/dlat) + 1
                if floating_edge:
                    if self.is_cartesian():
                        raise Exception('Grid is cartesian, so cant set exact dlon/dlat using floating_edge!')
                    x_end = self.edges('lon')[0]+(nx-1)*dlon
                    y_end = self.edges('lat')[0]+(ny-1)*dlat
                return int(nx), int(ny), x_end, y_end

            if dm:
                dx = dm
                dy = dm

            if dx and dy:
                nx = np.round((self.edges('x')[1]-self.edges('x')[0])/dx) + 1
                ny = np.round((self.edges('y')[1]-self.edges('y')[0])/dy) + 1
                if floating_edge:
                    if not self.is_cartesian():
                        raise Exception('Grid is spherical, so cant set exact dx/dy using floating_edge!')
                    x_end = self.edges('x')[0]+(nx-1)*dx
                    y_end = self.edges('y')[0]+(ny-1)*dy
                return int(nx), int(ny), x_end, y_end

            raise ValueError('Give a combination of nx/xy, dlon/dlat, dx/dy or dm')

        nx, ny, native_x_end, native_y_end = determine_nx_ny(nx,ny, dx, dy, dlon, dlat)
        x_native = np.linspace(self.native_x()[0], native_x_end, nx)
        y_native = np.linspace(self.native_y()[0], native_y_end, ny)
        if self.is_cartesian():
            x = x_native
            y = y_native
            lon = None
            lat = None
        else:
            lon = x_native
            lat = y_native
            x = None
            y = None
        self._init_structure(x, y, lon, lat)
        print(self)

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
        else:
            self._raw = UnstrGrid(lon=lon, lat=lat, x=x, y=y)

        if self.edges('lon', native=True)[0] < self.raw().edges(self.x_str)[0] or self.edges('lon', native=True)[1] > self.raw().edges(self.x_str)[1]:
            msg.warning(f"The data gotten from the TopoReader doesn't cover the grid in the {self.x_str} direction. Grid: {self.edges('lon')}, imported topo: {self.raw().edges(self.x_str)}")

        if self.edges('lat', native=True)[0] < self.raw().edges(self.y_str)[0] or self.edges('lat', native=True)[1] > self.raw().edges(self.y_str)[1]:
            msg.warning(f"The data gotten from the TopoReader doesn't cover the grid in the {self.y_str} direction. Grid: {self.edges('lat')}, imported topo: {self.raw().edges(self.y_str)}")



        if zone_number is not None:
            self._raw.set_utm(zone_number, zone_letter)

        self._raw._update_datavar('topo', topo)
        self._raw._update_sea_mask()

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

    def time(self) -> tuple:
        return (None, None)

    def raw(self):
        if hasattr(self, '_raw'):
            return self._raw
        return None

    def __str__(self) -> str:
        """Prints status of the grid."""

        msg.header(self, f"Status of grid {self.name}")
        msg.plain(f'lon: {self.lon()[0]} - {self.lon()[-1]}, lat: {self.lat()[0]} - {self.lat()[-1]}')

        if self.dlon() is not None:
            msg.plain(f'dlon, dlat = {self.dlon()}, {self.dlat()} deg')
            if self.dlon() > 0 and self.dlat() > 0:
                msg.plain(f'Inverse: dlon, dlat = 1/{1/self.dlon()}, 1/{1/self.dlat()} deg')

        if self.dx() is not None:
            msg.plain(f'dx, dy approximately {self.dx()}, {self.dy()} metres')

        if self.nx() is not None:
            msg.plain(f'nx, ny = {self.nx()} x {self.ny()} grid points')

        if not self.is_empty('topo'):
            msg.plain(f"Mean depth: {np.mean(self.topo()[self.sea_mask()]):.1f} m")
            msg.plain(f"Max depth: {np.max(self.topo()[self.sea_mask()]):.1f} m")
            msg.plain(f"Min depth: {np.min(self.topo()[self.sea_mask()]):.1f} m")

        if not self.is_empty('topo'):
            msg.print_line()
            msg.plain('Grid contains:')
            msg.plain(f'{sum(sum(self.sea_mask())):d} sea points')
            msg.plain(f'{sum(sum(np.logical_not(self.sea_mask()))):d} land points')
            msg.plain(f'{sum(sum(np.logical_and(self.boundary_mask(), self.sea_mask()))):d} boundary points')

        msg.print_line()

        return ''


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

@topography_methods
@add_datavar(name='topo', default_value=999., stash_get=True)
@add_mask(name='buoy', coords='grid', default_value=0)
@add_mask(name='spec', coords='grid', default_value=0)
@add_mask(name='boundary', coords='grid', default_value=0)
@add_mask(name='sea', coords='grid', default_value=1, opposite_name='land')
class UnstrGrid(PointSkeleton):

    @classmethod
    def from_grid(cls, grid: Union[GriddedSkeleton, PointSkeleton], name: str=None, mask: str=None):
        if name is None:
            name = grid.name
        if mask is None:
            lon, lat = grid.lonlat(strict=True)
            x, y = grid.xy(strict=True)
        else:
            if f"{mask}_mask" not in grid.masks():
                raise KeyError(f"Grid doesn't have a {mask} mask!")
            lon, lat = eval(f"grid.{mask}_points(type='lon', strict=True)")
            x, y = eval(f"grid.{mask}_points(type='x', strict=True)")
        unstr_grid = cls(lon=lon, lat=lat, x=x, y=y, name=name)
        utm_number, utm_zone = grid.utm()
        unstr_grid.set_utm(utm_number, utm_zone)
        return unstr_grid
    
    def __init__(self, grid=None, x=None, y=None, lon=None, lat=None, name='AnonymousGrid'):
        self.name = name
        if grid is not None:
            x, y = grid.xy(strict=True)
            lon, lat = grid.lonlat(strict=True)
        self._init_structure(x, y, lon, lat)

    def import_topo(self, topo_reader: TopoReader) -> None:
        """Reads the raw bathymetrical data."""

        msg.header(topo_reader, "Importing topography...")
        print(topo_reader)
        topo, lon, lat, x, y, zone_number, zone_letter = topo_reader(self.edges('lon'), self.edges('lat'), self.edges('x'), self.edges('y'))
        if 0 in topo.shape:
            msg.warning('Imported topography seems to be empty. Maybe using wrong tile?')

        if aux_funcs.is_gridded(topo, lon, lat) or aux_funcs.is_gridded(topo, x, y):
            self._raw = Grid(lon=lon, lat=lat, x=x, y=y)
        else:
            self._raw = UnstrGrid(lon=lon, lat=lat, x=x, y=y)

        if self.edges('lon', native=True)[0] < self.raw().edges(self.x_str)[0] or self.edges('lon', native=True)[1] > self.raw().edges(self.x_str)[1]:
            msg.warning(f"The data gotten from the TopoReader doesn't cover the grid in the {self.x_str} direction. Grid: {self.edges('lon')}, imported topo: {self.raw().edges(self.x_str)}")

        if self.edges('lat', native=True)[0] < self.raw().edges(self.y_str)[0] or self.edges('lat', native=True)[1] > self.raw().edges(self.y_str)[1]:
            msg.warning(f"The data gotten from the TopoReader doesn't cover the grid in the {self.y_str} direction. Grid: {self.edges('lat')}, imported topo: {self.raw().edges(self.y_str)}")


        if zone_number is not None:
            self.raw().set_utm(zone_number, zone_letter)

        self.raw()._update_datavar('topo', topo)
        self.raw()._update_sea_mask()

    def tri(self):
        return None

    def raw(self):
        if hasattr(self, '_raw'):
            return self._raw
        return None

    def boundary_nx(self) -> int:
        """Return approximate number of grid points in the longitude direction
        """
        return np.round(len(self.inds())/len(self.boundary_points()[0])/2).astype(int)

    def boundary_ny(self) -> int:
        return np.round(len(self.inds())/len(self.boundary_points()[0])/2).astype(int)

class TriGrid(UnstrGrid):
    @classmethod
    def from_msh(cls, filename: str, name: str=None):
        if name is None:
            tri_grid = cls()
        else:
            tri_grid = cls(name=name)
        tri_grid.import_triang(triang_MshFile(filename))
        tri_grid.import_topo(topo_MshFile(filename))
       
        return tri_grid



    def __init__(self, x=None, y=None, lon=None, lat=None, name='AnonymousGrid'):
        self.name = name
        # Only initialize if x, y, lon, lat given
        if [a for a in (x, y, lon, lat) if a is not None]:
            self._init_structure(x, y, lon, lat)

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

    def tri(self):
        if hasattr(self, '_tri'):
            return copy(self._tri)
        else:
            return None

    def arange_triangulation(self, tri_aranger: TriAranger) -> None:
        print(tri_aranger)
        bnd_nodes, tri, nodes, x, y = tri_aranger(self.inds(),  np.where(self.boundary_mask())[0], self.tri(), self.native_x(), self.native_y())
        if self.x_str == 'x':
            self._init_structure(x=x, y=y, lon=None, lat=None)
        else:
            self._init_structure(x=None, y=None, lon=x, lat=y)
        self._update_boundary(bnd_nodes)
        self._tri = tri

    def _update_boundary(self, boundary_inds):
        mask = np.array([ind in boundary_inds for ind in self.inds()])
        self._update_mask('boundary', mask)
