from __future__ import annotations # For TYPE_CHECKING

import xarray as xr
import numpy as np
from copy import copy
import sys
import re

# Import abstract classes and needed instances of them
from .read import TopoReader, EmptyTopo
from .boundary import BoundarySetter, ClearBoundary
from .mesh import Mesher, Interpolate
from .process import GridProcessor, TrivialFilter
from typing import TYPE_CHECKING, Tuple
if TYPE_CHECKING:
    from .write import GridWriter

# Import default values and aux_funcsiliry functions
from .. import msg
from .. import file_module
from ..aux_funcs import force_to_xyz, distance_2points, domain_size_in_km, set_spacing_dlon_dlat_fixed_edges, set_spacing_dlon_dlat_floating_edges, set_spacing_dx_dy, set_spacing_nx_ny


class Grid:
    def __init__(self, lon: Tuple[float, float]=(0.,0.), lat: Tuple[float, float]=(0.,0.), name: str="AnonymousGrid") -> None:
        """Initializes a new grid by setting the bounding box and name"""

        coords_dict = {'lon': np.array(np.unique(lon)), 'lat': np.array(np.unique(lat))}
        attr_dict = {'name': name}
        self.data = xr.Dataset(
                    coords=coords_dict,
                    attrs=attr_dict
                    )

    def import_topo(self, topo_reader: TopoReader) -> None:
        """Reads the raw bathymetrical data."""

        msg.header(topo_reader, "Importing topography...")
        print(topo_reader)
        topo, lon, lat = topo_reader(self.lon()[0], self.lon()[-1], self.lat()[0], self.lat()[-1])

        topo, lon, lat = force_to_xyz(topo, lon, lat)

        # Depth is positive, so set everything that is not positive to nan
        #topo[topo<=0]=np.nan

        points = [x for x in range(len(lon))]
        coords_dict = {'points': points}
        vars_dict = {'topo': (['points'], topo), 'lon': (['points'], lon), 'lat': (['points'], lat)}
        self.rawdata = xr.Dataset(
                    coords=coords_dict,
                    data_vars=vars_dict
                    )

    def process_topo(self, filt: GridProcessor=TrivialFilter()) -> None:
        """Processes the raw bathymetrical data, e.g. with a filter."""

        msg.header(filt, "Filtering topography...")
        land_sea_mask = self.raw_topo() > 0 # Sea points set to true

        print(filt)
        topo = filt.topo(self.raw_topo(), self.raw_lon(), self.raw_lat(), land_sea_mask)

        if topo is None:
            msg.warning('Filtering of raw topography is not implemented')
        else:
            vars_dict = {'topo': (['points'], topo)}
            self.rawdata = self.rawdata.assign(vars_dict)

    def mesh_grid(self, mesher: Mesher=Interpolate(method = 'linear')) -> None:
        """Meshes the raw data down to the grid definitions."""

        if hasattr(self, 'lon') and hasattr(self, 'lon'):

            msg.header(mesher, "Meshing grid bathymetry...")
            print(mesher)
            lonQ, latQ = np.meshgrid(self.lon(), self.lat())
            topo = mesher(self.raw_topo(), self.raw_lon(), self.raw_lat(), lonQ, latQ)

            vars_dict = {'topo': (['lat', 'lon'], topo)}
            self.data = self.data.assign(vars_dict)

            self._update_masks()
            print(self)
        else:
            msg.templates('no_spacing')

    def process_grid(self, filt: GridProcessor=TrivialFilter()) -> None:
        """Processes the gridded bathymetrical data, e.g. with a filter."""

        msg.header(filt, "Filtering meshed grid...")
        print(filt)
        topo = filt.grid(self.topo(), self.lon(), self.lat(), self.land_sea_mask(), self.boundary_mask())

        if topo is None:
            msg.warning('Filtering of mesed topography is not implemented')
        else:
            vars_dict = {'topo': (['lat', 'lon'], topo)}
            self.data = self.data.assign(vars_dict)

        msg.info('Upodating land-sea mask and boundary mask')
        self._update_masks()

    def _reset_grid(self, dlon: float, dlat: float, dx: float, dy: float, lon: np.ndarray, lat: np.ndarray) -> None:
        # Old topography conflicts in size, so drop them first
        self._drop_topo_and_masks()

        # Set the new values to the xarray
        attr_dict = {'dx': dx, 'dy': dy, 'dlon': dlon, 'dlat': dlat}
        self.data = self.data.assign_attrs(attr_dict)
        coords_dict = {'lon': lon, 'lat': lat}
        self.data = self.data.assign_coords(coords_dict)

        # Initialize the grid with an empty topography
        msg.info("Initializing with an empty topography")
        self.import_topo(topo_reader=EmptyTopo(self))
        self.mesh_grid(mesher=Interpolate(method='nearest'))

    def set_spacing(self, dlon: float=0, dlat: float=0, dm: float=0, nx: int=0, ny: int=0, floating_edge: bool=False) -> None:
        """Defines longitude and latitude vectors based on desired spacing.

        Options (priority in this order)
        dlon, dlat [deg]:   Grid spacing is set as close to the given resolution
                            as possible (edges are fixed).

                            Set floating_edge=True to force exact dlon, dlat
                            and instead possibly move lon_max, lat_max slightly
                            to make it work.

        dm [m]:             Grid spacing is set as close as dm metres as
                            possible.

        nx, ny [grid pnts]: Grid resolution is set to have nx points in
                            longitude and ny points in latitude direction.
        """

        msg.header(self, "Setting grid spacing...")

        if dlon and dlat:
            msg.plain(f"Setting spacing based on dlon = {dlon} and dlat = {dlat}")

            if floating_edge: #Use exactly given dlon/dlat and change lon_max/lat_max accordingly
                msg.plain("floating_edge = True. Making sure dlon/dlat are keep exactly fixed")
                dlon, dlat, dx, dy, lon_array, lat_array = set_spacing_dlon_dlat_floating_edges(dlon=dlon, dlat=dlat, lon=self.lon_edges(), lat=self.lat_edges())
            else: # Keeping edges fixed and rounding dlon/dlat to something suitable
                # Number of points
                dlon, dlat, dx, dy, lon_array, lat_array = set_spacing_dlon_dlat_fixed_edges(dlon=dlon, dlat=dlat, lon=self.lon_edges(), lat=self.lat_edges())
        elif dm:
            msg.plain(f"Setting spacing based on (approximately) dm={dm} metres")
            dlon, dlat, dx, dy, lon_array, lat_array = set_spacing_dx_dy(dx=dm, dy=dm, lon=self.lon_edges(), lat=self.lat_edges())
        elif nx and ny:
            # Cant expand one point grid
            if self.nx()<2:
                nx = 1
            if self.ny()<2:
                ny=1

            msg.plain(f"Setting spacing to have nx = {nx}, ny = {ny} points.")
            dlon, dlat, dx, dy, lon_array, lat_array = set_spacing_nx_ny(nx=nx, ny=ny, lon=self.lon_edges(), lat=self.lat_edges())
        else:
            print(f'Can not set grid spacing to zero since no non-zero values given!')
            return

        self._reset_grid(dlon, dlat, dx, dy, lon_array, lat_array)

    def set_boundary(self, boundary_setter: BoundarySetter) -> None:
        """Marks the points that should be treated as boundary points in the
        grid.

        The boundary points are stored in a boolean array where True values
        mark a boundary point.

        NB! No check for land points are done, so it is possible that a land
        point is marked as a boundary point. Possibly accounting for this is
        the responsibility of the GridWriter.
        """

        msg.header(boundary_setter, "Setting boundary points...")
        print(boundary_setter)

        boundary_mask = boundary_setter(self.land_sea_mask())

        vars_dict = {'boundary_mask': (['lat', 'lon'], boundary_mask)}
        self.data = self.data.assign(vars_dict)

    def land_sea_mask(self) -> np.ndarray:
        """Returns bool array of the land-sea mask (True = sea point)"""

        if hasattr(self.data, 'land_sea_mask'):
            return copy(self.data.land_sea_mask.values)
        else:
            return np.full((self.ny(), self.nx()), True)

    def boundary_mask(self) -> np.ndarray:
        """Returns bool array of boundary points (True = boundary point)"""

        if hasattr(self.data, 'boundary_mask'):
            return copy(self.data.boundary_mask.values)
        else:
            return np.full((self.ny(), self.nx()), False)

    def boundary_points(self) -> np.ndarray:
        """Returns a lon, lat list of the set boundary points."""

        if self.boundary_mask().size > 0:
            mask = np.logical_and(self.boundary_mask(), self.land_sea_mask())
            BOUND = self._point_list(mask)
            return BOUND
        else:
            return np.array([])

    def land_points(self) -> np.ndarray:
        """Returns a lon, lat list of land points."""

        if self.boundary_mask().size > 0:
            mask = np.logical_not(self.land_sea_mask())
            LAND = self._point_list(mask)
            return LAND
        else:
            return np.array([])

    def sea_points(self) -> np.ndarray:
        """Returns a lon, lat list of sea points."""

        if self.boundary_mask().size > 0:
            mask = self.land_sea_mask()
            land = self._point_list(mask)
            return land
        else:
            return np.array([])

    def name(self) -> str:
        """Return the name of the grid (set at initialization)."""
        return copy(self.data.name)

    def structured(self) -> bool:
        return True

    def size(self) -> Tuple[int, int]:
        """Returns the size (nx, ny) of the grid."""
        #return self.land_sea_mask().shape
        return (self.ny(), self.nx())

    def topo(self, land: float=-999) -> np.ndarray:
        """Returns an array containing the meshed topography of the grid."""
        if hasattr(self.data, 'topo'):
            topo = copy(self.data.topo.values)
            topo[np.logical_not(self.land_sea_mask())] = land
            return topo
        else:
            return np.array([])

    def nx(self) -> int:
        """Return the number of points in longitude direction."""
        return len(self.lon())

    def ny(self) -> int:
        """Return the number of points in latitude direction."""
        return len(self.lat())

    def boundary_nx(self) -> int:
        """Return approximate number of grid points in the longitude direction
        """
        abs_diff = np.abs(np.diff(np.where(self.boundary_mask())))
        if abs_diff.size == 0:
            return 0
        abs_diff=np.median(abs_diff[abs_diff>0]).astype(int)

        return np.ceil(self.nx()/abs_diff+1).astype(int)

    def boundary_ny(self) -> int:
        """Return approximate number of grid points in the longitude direction
        """
        abs_diff = np.abs(np.diff(np.where(self.boundary_mask())))
        if abs_diff.size == 0:
            return 0
        abs_diff=np.median(abs_diff[abs_diff>0]).astype(int)

        return np.ceil(self.ny()/abs_diff+1).astype(int)


    def lon(self) -> np.ndarray:
        """Returns a longitude vector of the grid."""
        return copy(self.data.lon.values)

    def lat(self) -> np.ndarray:
        """Returns a latitude vector of the grid."""
        return copy(self.data.lat.values)

    def lon_edges(self) -> Tuple[float, float]:
        return (np.min(self.lon()), np.max(self.lon()))

    def lat_edges(self) -> Tuple[float, float]:
        return (np.min(self.lat()), np.max(self.lat()))

    def dlon(self) -> float:
        if hasattr(self.data, 'dlon'):
            return copy(self.data.dlon)
        else:
            return None

    def dlat(self) -> float:
        if hasattr(self.data, 'dlat'):
            return copy(self.data.dlat)
        else:
            return None

    def dx(self) -> float:
        if hasattr(self.data, 'dx'):
            return copy(self.data.dx)
        else:
            return None

    def dy(self) -> float:
        if hasattr(self.data, 'dy'):
            return copy(self.data.dy)
        else:
            return None

    def raw_topo(self) -> np.ndarray:
        """Returns an array containing the unmeshed imported topography."""
        if hasattr(self, 'rawdata'):
            return copy(self.rawdata.topo.values)
        else:
            return np.array([])

    def raw_lon(self) -> np.ndarray:
        """Returns a longitude vector of the unmeshed imported topography."""
        if hasattr(self, 'rawdata'):
            return copy(self.rawdata.lon.values)
        else:
            return np.array([])

    def raw_lat(self) -> np.ndarray:
        """Returns a latitude vector of the unmeshed imported topography."""
        if hasattr(self, 'rawdata'):
            return copy(self.rawdata.lat.values)
        else:
            return np.array([])

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

    def _update_masks(self) -> None:
        """Sets land-sea mask and boundary point mask.

        This is called after the data is meshed to the grid or the gridded data
        is processed in order to make sure that everything is consistent.
        """

        self._set_land_sea_mask(land_sea_mask = self.data.topo.values > -20) # Only positive values

        # Create empty (no boundary points) if doesn't exist
        if self.boundary_mask().size == 0:
            self.set_boundary(boundary_setter = ClearBoundary())
        return

    def _drop_topo_and_masks(self) -> None:
        """Drops the gridded data and masks."""

        if hasattr(self.data, 'topo'):
            self.data = self.data.drop('topo')
        if hasattr(self.data, 'land_sea_mask'):
            self.data = self.data.drop('land_sea_mask')
        if hasattr(self.data, 'boundary_mask'):
            self.data = self.data.drop('boundary_mask')
        return

    def _set_land_sea_mask(self, land_sea_mask) -> None:
        """Sets the land-sea mask based on land points in the gridded
        bathymetrical data.

        A certain land-sea mask can also be forced by providing a boolean
        matrix by the user (as a bit of hack). Note, that this will be
        overwritten automatically by the method _update_masks() if any further
        changes are made to the topography.
        """

        vars_dict = {'land_sea_mask': (['lat', 'lon'], land_sea_mask)}
        self.data = self.data.assign(vars_dict)

        return

    def _point_list(self, mask=None):
        """Provides a list on longitudes and latitudes with a given mask.

        Used to e.g. generate list of boundary points or land points.
        """

        if mask is None:
            mask = np.ones(self.size(), dtype=bool)
        meshlon, meshlat=np.meshgrid(self.lon(),self.lat())
        lonlat_flat = np.column_stack((meshlon.ravel(),meshlat.ravel()))
        mask_flat = mask.ravel()

        return lonlat_flat[mask_flat]

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

    def __str__(self) -> str:
        """Prints status of the grid."""

        if self.topo().shape==(0,):
            empty_topo = True
        elif np.mean(self.topo()[self.land_sea_mask()]) == 9999:
            empty_topo = True
        else:
            empty_topo = False

        msg.header(self, f"Status of grid {self.name()}")
        msg.plain(f'lon: {self.lon()[0]} - {self.lon()[-1]}, lat: {self.lat()[0]} - {self.lat()[-1]}')

        if self.dlon() is not None:
            msg.plain(f'dlon, dlat = {self.dlon()}, {self.dlat()} deg')
            if self.dlon() > 0 and self.dlat() > 0:
                msg.plain(f'Inverse: dlon, dlat = 1/{1/self.dlon()}, 1/{1/self.dlat()} deg')

        if self.dx() is not None:
            msg.plain(f'dx, dy approximately {self.dx()}, {self.dy()} metres')

        if self.nx() is not None:
            msg.plain(f'nx, ny = {self.nx()} x {self.ny()} grid points')

        if self.topo().size > 0 and (not empty_topo):
            msg.plain(f"Mean depth: {np.mean(self.topo()[self.land_sea_mask()]):.1f} m")
            msg.plain(f"Max depth: {np.max(self.topo()[self.land_sea_mask()]):.1f} m")
            msg.plain(f"Min depth: {np.min(self.topo()[self.land_sea_mask()]):.1f} m")

        if self.land_sea_mask().size > 0:
            msg.print_line()
            msg.plain('Grid contains:')
            msg.plain(f'{sum(sum(self.land_sea_mask())):d} sea points')
            msg.plain(f'{sum(sum(np.logical_not(self.land_sea_mask()))):d} land points')

        if self.boundary_mask().size > 0:
            msg.plain(f'{sum(sum(np.logical_and(self.boundary_mask(), self.land_sea_mask()))):d} boundary points')

        msg.print_line()

        return ''

    def __repr__(self):
        lines = [f"<dnora Grid object> (structured)", f"  Name: {self.name()}"]

        if self.topo().shape==(0,):
            empty_topo = True
        elif np.mean(self.topo()[self.land_sea_mask()]) == 9999:
            empty_topo = True
        else:
            empty_topo = False

        if self.raw_topo().shape==(0,):
            empty_raw_topo = True
        elif np.mean(self.raw_topo()) == 9999:
            empty_raw_topo = True
        else:
            empty_raw_topo = False

        if self.nx()>2 or self.ny()>2:
            lines.append(f"  Number of points: (nx={self.nx()}, ny={self.ny()})")
        else:
            lines.append(f"  Number of points: Use method .set_spacing() to set structure.")
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
