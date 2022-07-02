from gridded_skeleton import GriddedSkeleton
from unstr_grd_mod import UnstrGrid
from topography import Topography
import numpy as np
import xarray as xr
from dnora.grd.read import TopoReader
from dnora import msg
from dnora.grd.mesh import Mesher, Interpolate
from coordinate_factory import add_time
from mask_factory import add_mask
from datavar_factory import add_datavar


def is_gridded(data: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> bool:
    if data.shape == (len(lat), len(lon)):
        return True

    if len(data.shape) == 1 and len(lat) == data.shape[0] and len(lon) == data.shape[0]:
        return False

    raise Exception(f"Size of data is {data.shape} but len(lat) = {len(lat)} and len(lon) = {len(lon)}. I don't know what is going on!")

@add_datavar(name='topo', default_value=999.)
@add_mask(name='boundary', coords='grid', default_value=0)
@add_mask(name='sea', coords='grid', default_value=1)
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
                return int(nx), int(ny),

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

        x = np.linspace(self.native_x()[0], native_x_end, nx)
        y = np.linspace(self.native_y()[0], native_y_end, ny)
        self._init_structure(x, y)
        print(self)
    # def set_boundary(self, boundary_setter: BoundarySetter) -> None:
    #     """Marks the points that should be treated as boundary points in the
    #     grid.
    #
    #     The boundary points are stored in a boolean array where True values
    #     mark a boundary point.
    #
    #     NB! No check for land points are done, so it is possible that a land
    #     point is marked as a boundary point. Possibly accounting for this is
    #     the responsibility of the GridWriter.
    #     """
    #
    #     msg.header(boundary_setter, "Setting boundary points...")
    #     print(boundary_setter)
    #
    #     boundary_mask = boundary_setter(self.sea_mask()).astype(int)
    #
    #     self.ds_manager.update_mask('boundary', boundary_mask)


    def import_topo(self, topo_reader: TopoReader) -> None:
        """Reads the raw bathymetrical data."""

        if isinstance(topo_reader, Grid) or isinstance(topo_reader, UnstrGrid):
            msg.header(topo_reader, "Getting topography from Grid-object...")
            lon, lat = topo_reader.lon(), topo_reader.lat()
            topo = topo_reader.topo()
        else:
            msg.header(topo_reader, "Importing topography...")
            print(topo_reader)
            topo, lon, lat = topo_reader(self.lon()[0], self.lon()[-1], self.lat()[0], self.lat()[-1])

        if is_gridded(topo, lon, lat):
            self.raw = Grid(lon=lon, lat=lat)
        else:
            self.raw = UnstrGrid(lon=lon, lat=lat)

        self.raw.ds_manager.update_datavar('topo', topo)
        self.raw._update_sea_mask()

    def mesh_grid(self, mesher: Mesher=Interpolate(method = 'linear')) -> None:
        """Meshes the raw data down to the grid definitions."""

        msg.header(mesher, "Meshing grid bathymetry...")
        print(mesher)
        lonQ, latQ = np.meshgrid(self.lon(), self.lat())
        lon, lat = self.raw.lonlat()

        topo = mesher(self.raw.topo().ravel(), lon, lat, lonQ, latQ)
        self.ds_manager.update_datavar('topo', topo)
        self._update_sea_mask()
#            print(self)

    def _update_sea_mask(self):
        self.ds_manager.update_mask('sea', (self.topo()>0).astype(int))

    def __str__(self) -> str:
        """Prints status of the grid."""

        # if self.topo().shape==(0,):
        #     empty_topo = True
        # elif np.mean(self.topo()[self.land_sea_mask()]) == 9999:
        #     empty_topo = True
        # else:
        #     empty_topo = False

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

        # if self.topo().size > 0 and (not empty_topo):
        #     msg.plain(f"Mean depth: {np.mean(self.topo()[self.land_sea_mask()]):.1f} m")
        #     msg.plain(f"Max depth: {np.max(self.topo()[self.land_sea_mask()]):.1f} m")
        #     msg.plain(f"Min depth: {np.min(self.topo()[self.land_sea_mask()]):.1f} m")
        #
        # if self.land_sea_mask().size > 0:
        #     msg.print_line()
        #     msg.plain('Grid contains:')
        #     msg.plain(f'{sum(sum(self.land_sea_mask())):d} sea points')
        #     msg.plain(f'{sum(sum(np.logical_not(self.land_sea_mask()))):d} land points')
        #
        # if self.boundary_mask().size > 0:
        #     msg.plain(f'{sum(sum(np.logical_and(self.boundary_mask(), self.land_sea_mask()))):d} boundary points')

        msg.print_line()

        return ''
