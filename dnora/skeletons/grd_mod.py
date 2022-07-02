from gridded_skeleton import GriddedSkeleton
from unstr_grd_mod import UnstrGrid
from topography import Topography
import numpy as np
import xarray as xr
from dnora.grd.read import TopoReader
from dnora import msg
from dnora.grd.mesh import Mesher, Interpolate
from coordinates import include_time
from mask_generator import add_mask

def is_gridded(data: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> bool:
    if data.shape == (len(lat), len(lon)):
        return True

    if len(data.shape) == 1 and len(lat) == data.shape[0] and len(lon) == data.shape[0]:
        return False

    raise Exception(f"Size of data is {data.shape} but len(lat) = {len(lat)} and len(lon) = {len(lon)}. I don't know what is going on!")

@add_mask(name='boundary', coords='grid', default_value=0)
@add_mask(name='land', coords='grid', default_value=0)
@add_mask(name='sea', coords='grid', default_value=1)
class Grid(GriddedSkeleton, Topography):
    def __init__(self, x=None, y=None, lon=None, lat=None, name='AnonymousGrid'):
        self.name = name
        self._create_structure(x, y, lon, lat)
        self._reset_vars()

    def set_spacing(self, dlon: float=0., dlat: float=0.,
                    dx: float=0., dy: float=0., dm: float=0,
                    nx: int=0, ny: int=0, floating_edge: bool=False) -> None:
        """Defines longitude and latitude vectors based on desired spacing.

        Options (priority in this order)
        nx, ny [grid pnts]: Grid resolution is set to have nx points in
                            longitude and ny points in latitude direction.

        dlon, dlat [deg]:   Grid spacing is set as close to the given resolution
                            as possible (edges are fixed).

                            Set floating_edge=True to force exact dlon, dlat
                            and instead possibly move lon_max, lat_max slightly
                            to make it work.

        dm [m]:             Grid spacing is set as close as dm metres as
                            possible.

        dx, dy [m]:         Grid spacing is set as close as dx and xy metres as
                            possible.
        """
        def determine_nx_ny(dlon, dlat, dx, dy, dm, nx, ny):
            if nx and ny:
                return int(nx), int(ny)

            if dlon and dlat:
                nx = np.round((self.lon_edges()[1]-self.lon_edges()[0])/dlon) + 1
                ny = np.round((self.lat_edges()[1]-self.lat_edges()[0])/dlat) + 1
                return int(nx), int(ny)

            if dm:
                dx = dm
                dy = dm

            if dx and dy:
                nx = np.round((self.x_edges()[1]-self.x_edges()[0])/dx) + 1
                ny = np.round((self.y_edges()[1]-self.y_edges()[0])/dy) + 1
                return int(nx), int(ny)

            raise ValueError('Give a combination of nx/xy, dlon/dlat, dx/dy or dm')

        nx, ny = determine_nx_ny(dlon, dlat, dx, dy, dm, nx, ny)

        x = np.linspace(self.native_x_edges()[0], self.native_x_edges()[1], nx)
        y = np.linspace(self.native_y_edges()[0], self.native_y_edges()[1], ny)
        self.data = self._init_ds(x, y)
        self._reset_vars()

    def import_topo(self, topo_reader: TopoReader) -> None:
        """Reads the raw bathymetrical data."""

        msg.header(topo_reader, "Importing topography...")
        print(topo_reader)
        topo, lon, lat = topo_reader(self.lon()[0], self.lon()[-1], self.lat()[0], self.lat()[-1])

        if is_gridded(topo, lon, lat):
            self.raw = Grid(lon=lon, lat=lat)
        else:
            self.raw = UnstrGrid(lon=lon, lat=lat)
        topo[topo<=0]=np.nan

        self.raw._set(topo, 'topo')
        self.raw._update_sea_mask()

    def mesh_grid(self, mesher: Mesher=Interpolate(method = 'linear')) -> None:
        """Meshes the raw data down to the grid definitions."""

        msg.header(mesher, "Meshing grid bathymetry...")
        print(mesher)
        lonQ, latQ = np.meshgrid(self.lon(), self.lat())
        lon, lat = self.raw.lonlat()

        topo = mesher(self.raw.topo().ravel(), lon, lat, lonQ, latQ)
        topo[topo<=0]=np.nan
        self._set(topo, 'topo')
        self._update_sea_mask()
#            print(self)
