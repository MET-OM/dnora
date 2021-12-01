from __future__ import annotations # For TYPE_CHECKING
import xarray as xr
import numpy as np
from copy import copy

import sys
import re
from .. import msg
from ..aux import distance_2points, add_prefix, add_suffix, add_file_extension, create_filename_obj, add_folder_to_filename, clean_filename, check_if_folder
from ..defaults import dflt_grd, list_of_placeholders
from .read import TopoReader, EmptyTopo
from .boundary import BoundarySetter, ClearBoundary
from .mesh import Mesher, TrivialMesher, Interpolate
from .process import GridProcessor, TrivialFilter

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .write import GridWriter # Abstract class

class Grid:
    def __init__(self, lon_min: float = 0., lon_max: float = 0., lat_min: float = 0., lat_max: float = 0., name: str = "AnonymousGrid"):
        """Initializes a new grid by setting the bounding box and name"""

        data_dict = {'lon_min': lon_min, 'lon_max': lon_max, 'lat_min': lat_min, 'lat_max': lat_max, 'name': name}
        self.data = xr.Dataset(
                    attrs=(data_dict
                    ),
                    )

        return

    def import_topo(self, topo_reader: TopoReader) -> None:
        """Reads the raw bathymetrical data."""

        msg.header(topo_reader, "Importing topography...")
        print(topo_reader)
        topo, lon, lat = topo_reader(self.data.lon_min, self.data.lon_max, self.data.lat_min, self.data.lat_max)

        coords_dict = {'lon': lon, 'lat': lat}
        vars_dict = {'topo': (['lat', 'lon'], topo)}
        self.rawdata = xr.Dataset(
                    coords=(coords_dict
                    ),
                    data_vars=(vars_dict
                    ),
                    )
        return

    def process_topo(self, filt: GridProcessor = TrivialFilter()) -> None:
        """Processes the raw bathymetrical data, e.g. with a filter."""

        msg.header(filt, "Filtering topography...")

        empty_mask = np.full(self.raw_topo().shape, False)
        land_sea_mask = self.raw_topo() < 0 # Sea points set to true

        print(filt)
        topo = filt(self.raw_topo(), self.raw_lon(), self.raw_lat(), land_sea_mask, empty_mask)

        vars_dict = {'topo': (['lat', 'lon'], topo)}
        self.rawdata = self.rawdata.assign(vars_dict)

        return

    def mesh_grid(self, mesher: Mesher=Interpolate(method = 'linear')) -> None:
        """Meshes the raw data down to the grid definitions."""

        if hasattr(self, 'lon') and hasattr(self, 'lon'):

            msg.header(mesher, "Meshing grid bathymetry...")
            print(mesher)
            topo = mesher(self.raw_topo(), self.raw_lon(), self.raw_lat(), self.lon(), self.lat())
            vars_dict = {'topo': (['lat', 'lon'], topo)}
            self.data = self.data.assign(vars_dict)

            self._update_masks()
            print(self)
            return

        else:
            msg.templates('no_spacing')
            return

    def process_grid(self, filt: GridProcessor = TrivialFilter()) -> None:
        """Processes the gridded bathymetrical data, e.g. with a filter."""

        msg.header(filt, "Filtering meshed grid...")
        print(filt)
        topo = filt(self.topo(), self.lon(), self.lat(), self.land_sea_mask(), self.boundary_mask())

        vars_dict = {'topo': (['lat', 'lon'], topo)}
        self.data = self.data.assign(vars_dict)

        msg.info('Upodating land-sea mask and boundary mask')
        self._update_masks()

        return

    def export_grid(self, grid_writer: GridWriter, out_format: str=None, filestring: str=None, infofilestring: str=None, folder: str=None) -> None:
        """Exports the boundary spectra to a file.

        The grid_writer defines the file format.
        """

        msg.header(grid_writer, f"Writing grid topography from {self.name()}")

        if out_format is None:
            out_format = grid_writer._preferred_format()

        if filestring is not None:
            filename = self.filename(filestring=filestring)
        else:
            filename = self.filename(filestring=dflt_grd['fs'][out_format])

        if infofilestring is not None:
            infofilename = self.filename(filestring=infofilestring)
        else:
            infofilename = self.filename(filestring=dflt_grd['info'][out_format])

        if folder is not None:
            folder = self.filename(filestring=folder)
        else:
            folder = self.filename(filestring=dflt_grd['fldr'][out_format])

        existed = check_if_folder(folder=folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {folder}")

        output_files, output_folder = grid_writer(self, filename=filename, infofilename=infofilename, folder=folder)

        return output_files, output_folder

    def filename(self, filestring: str=dflt_grd['fs']['General']) -> str:
        """Creates a filename for the object.

        The filename can be based on e.g. the name of the Grid object itself,
        or the start and end times.

        This is typically called by a GridWriter object when using
        the .export_grid() method.
        """

        # Substitute placeholders for $Grid
        filename = create_filename_obj(filestring=filestring, objects=[self])

        return filename

    def folder(self, folderstring: str=dflt_grd['fldr']['General']) -> str:
        # Substitute placeholders for $Grid
        folder = create_filename_obj(filestring=folderstring, objects=[self])
        folder = clean_filename(folder, list_of_placeholders)

        return folder

    def set_spacing(self, dlon: float=0, dlat: float=0, dm: float=0, nx: int=0, ny: int=0, floating_edge: bool=False) -> None:
        """Defines longitude and latitude vectors based on desired spacing.

        Options
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

        # Resetting spacing will lose all information about the old grid.
        # We therefore check that
        reset_grid = False

        if dlon and dlat:
            msg.plain(f"Setting spacing based on dlon = {dlon} and dlat = {dlat}")

            if dm:
                msg.plain(f"dlon and dlat given. Ignoring value dm={dm} metres!")

            if floating_edge: #Use exactly given dlon/dlat and change lon_max/lat_max accordingly
                msg.plain("floating_edge = True. Making sure dlon/dlat are keep exactly fixed")

                lon=np.arange(self.data.lon_min,self.data.lon_max+dlon/2,dlon)
                lat=np.arange(self.data.lat_min,self.data.lat_max+dlon/2,dlat)

                msg.plain(f"Setting lon_max ({self.lon()[-1]} >> {lon[-1]}), lat_max ({self.lat()[-1]} >> {lat[-1]})")

                distance_x = distance_2points((lat[0]+lat[-1])/2, lon[0], (lat[0]+lat[-1])/2, lon[-1])
                distance_y = distance_2points(lat[0], lon[0], lat[-1], lon[0])

                # Number of points
                nx = len(lon)
                ny = len(lat)

                # dx, dy in metres
                dx = distance_x*1000/nx
                dy = distance_y*1000/ny


                attr_dict = {'lon_max': lon[-1], 'lat_max': lat[-1]}
                self.data = self.data.assign_attrs(attr_dict)

                reset_grid = True

            else: # Keeping edges fixed and rounding dlon/dlat to something suitable
                # Number of points
                nx = int((self.data.lon_max-self.data.lon_min)/dlon + 1)
                ny = int((self.data.lat_max-self.data.lat_min)/dlat + 1)

                # Define longitudes and latitudes
                lon = np.linspace(self.data.lon_min, self.data.lon_max, nx)
                lat = np.linspace(self.data.lat_min, self.data.lat_max, ny)
                dlon = (self.data.lon_max-self.data.lon_min)/(nx-1)
                dlat = (self.data.lat_max-self.data.lat_min)/(ny-1)

                distance_x = distance_2points((self.data.lat_min+self.data.lat_max)/2, self.data.lon_min, (self.data.lat_min+self.data.lat_max)/2, self.data.lon_max)
                distance_y = distance_2points(self.data.lat_min, self.data.lon_min, self.data.lat_max, self.data.lon_min)

                # dx, dy in metres
                dx = distance_x*1000/nx
                dy = distance_y*1000/ny

                reset_grid = True

        elif dm:
            msg.plain(f"Setting spacing based on (approximately) dm={dm} metres")
            distance_x = distance_2points((self.data.lat_min+self.data.lat_max)/2, self.data.lon_min, (self.data.lat_min+self.data.lat_max)/2, self.data.lon_max)
            distance_y = distance_2points(self.data.lat_min, self.data.lon_min, self.data.lat_max, self.data.lon_min)

            # Number of points
            nx = int(np.round(distance_x*1000/dm))
            ny = int(np.round(distance_y*1000/dm))

            # dx, dy in metres
            dx = distance_x*1000/nx
            dy = distance_y*1000/ny

            # Define longitudes and latitudes
            lon = np.linspace(self.data.lon_min, self.data.lon_max, nx)
            lat = np.linspace(self.data.lat_min, self.data.lat_max, ny)
            dlon = (self.data.lon_max-self.data.lon_min)/(nx-1)
            dlat = (self.data.lat_max-self.data.lat_min)/(ny-1)

            reset_grid = True

        elif nx and ny:
            msg.plain(f"Setting spacing to have nx = {nx}, ny = {ny} points.")

            # Define longitudes and latitudes
            lon = np.linspace(self.data.lon_min, self.data.lon_max, nx)
            lat = np.linspace(self.data.lat_min, self.data.lat_max, ny)
            dlon = (self.data.lon_max-self.data.lon_min)/(nx-1)
            dlat = (self.data.lat_max-self.data.lat_min)/(ny-1)

            distance_x = distance_2points((self.data.lat_min+self.data.lat_max)/2, self.data.lon_min, (self.data.lat_min+self.data.lat_max)/2, self.data.lon_max)
            distance_y = distance_2points(self.data.lat_min, self.data.lon_min, self.data.lat_max, self.data.lon_min)

            # dx, dy in metres
            dx = distance_x*1000/nx
            dy = distance_y*1000/ny

            reset_grid = True

        if reset_grid:
            # Old topography conflicts in size, so drop them first
            self._drop_topo_and_masks()

            # Set the new values to the xarray
            attr_dict = {'nx': nx, 'ny': ny, 'dx': dx, 'dy': dy, 'dlon': dlon, 'dlat': dlat}
            self.data = self.data.assign_attrs(attr_dict)

            coords_dict = {'lon': lon, 'lat': lat}
            self.data = self.data.assign_coords(coords_dict)

            # Initialize the grid with an empty topography
            msg.info("Initializing with an empty topography")
            self.import_topo(topo_reader=EmptyTopo(self))

            self.mesh_grid(mesher=TrivialMesher())

        else:
            msg.advice("Doing nothing. Run set_spacing with either dlon AND dlat (in degrees), nx AND ny (in grid points), or dm (in metres).")

        return

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

        boundary_mask = boundary_setter(self.land_sea_mask().shape)

        vars_dict = {'boundary_mask': (['lat', 'lon'], boundary_mask)}
        self.data = self.data.assign(vars_dict)

        return

    def land_sea_mask(self):
        """Returns bool array of the land-sea mask (True = sea point)"""

        if hasattr(self.data, 'land_sea_mask'):
            return copy(self.data.land_sea_mask.values)
        else:
            return np.array([])

    def boundary_mask(self):
        """Returns bool array of boundary points (True = boundary point)"""

        if hasattr(self.data, 'boundary_mask'):
            return copy(self.data.boundary_mask.values)
        else:
            return np.array([])

    def boundary_points(self):
        """Returns a lon, lat list of the set boundary points."""

        if self.boundary_mask().size > 0:
            mask = np.logical_and(self.boundary_mask(), self.land_sea_mask())
            BOUND = self._point_list(mask)
            return BOUND
        else:
            return np.array([])

    def land_points(self):
        """Returns a lon, lat list of land points."""

        if self.boundary_mask().size > 0:
            mask = np.logical_not(self.land_sea_mask())
            LAND = self._point_list(mask)
            return LAND
        else:
            return np.array([])

    def sea_points(self):
        """Returns a lon, lat list of sea points."""

        if self.boundary_mask().size > 0:
            mask = self.land_sea_mask()
            LAND = self._point_list(mask)
            return LAND
        else:
            return np.array([])

    def name(self) -> str:
        """Return the name of the grid (set at initialization)."""
        return copy(self.data.name)

    def size(self) -> tuple:
        """Returns the size (nx, ny) of the grid."""
        return self.land_sea_mask().shape

    def topo(self):
        """Returns an array containing the meshed topography of the grid."""
        if hasattr(self.data, 'topo'):
            return copy(self.data.topo.values)
        else:
            return np.array([])

    def nx(self) -> int:
        """Return the number of points in longitude direction."""
        if hasattr(self.data, 'nx'):
            return int(copy(self.data.nx))
        else:
            return 0

    def ny(self) -> int:
        """Return the number of points in latitude direction."""
        if hasattr(self.data, 'ny'):
            return int(copy(self.data.ny))
        else:
            return 0

    def lon(self):
        """Returns a longitude vector of the grid."""
        if hasattr(self.data, 'lon'):
            lon = copy(self.data.lon.values)
        else:
            lon = np.array([self.data.lon_min, self.data.lon_max])
        return lon

    def lat(self):
        """Returns a latitude vector of the grid."""
        if hasattr(self.data, 'lat'):
            lat = copy(self.data.lat.values)
        else:
            lat = np.array([self.data.lat_min, self.data.lat_max])
        return lat

    def raw_topo(self):
        """Returns an array containing the unmeshed imported topography."""
        return copy(self.rawdata.topo.values)

    def raw_lon(self):
        """Returns a longitude vector of the unmeshed imported topography."""
        return copy(self.rawdata.lon.values)

    def raw_lat(self):
        """Returns a latitude vector of the unmeshed imported topography."""
        return copy(self.rawdata.lat.values)

    def _update_masks(self) -> None:
        """Sets land-sea mask and boundary point mask.

        This is called after the data is meshed to the grid or the gridded data
        is processed in order to make sure that everything is consistent.
        """

        self._set_land_sea_mask(land_sea_mask = self.topo() > 0) # Land points -999

        # Create empty (no boundary points) if doesn't exist
        if self.boundary_mask().size == 0:
            self.set_boundary(boundary_setter = ClearBoundary())
        return

    def _drop_topo_and_masks(self) -> None:
        """Drops the gridded data and masks."""

        if self.topo().size > 0:
            self.data =self.data.drop('topo')
        if self.land_sea_mask().size > 0:
            self.data =self.data.drop('land_sea_mask')
        if self.boundary_mask().size > 0:
            self.data =self.data.drop('boundary_mask')
        return

    def _set_land_sea_mask(self, land_sea_mask) -> None:
        """Sets the land-sea mask based on land points in the gridded
        bathymetrical data.

        A certain land-sea mask can also be forced by providing a boolean
        matrix by the user (as a bit of hack). Note, that this will be
        overwritten automatically by the method _update_masks() if any further
        changes are made to the topography.
        """

        # if matrix is None:
        #     land_sea_mask = np.full(self.topo().shape, False)
        #     land_sea_mask = self.topo() > 0 # Sea points set to true
        # else:
        #     land_sea_mask = matrix

        vars_dict = {'land_sea_mask': (['lat', 'lon'], land_sea_mask)}
        self.data = self.data.assign(vars_dict)

        return

    def _point_list(self, mask):
        """Provides a list on longitudes and latitudes with a given mask.

        Used to e.g. generate list of boundary points or land points.
        """
        meshlon, meshlat=np.meshgrid(self.lon(),self.lat())
        lonlat_flat = np.column_stack((meshlon.ravel(),meshlat.ravel()))
        mask_flat = mask.ravel()

        return lonlat_flat[mask_flat]

    def write_status(self, filename='', folder='') -> None:
        """Writes out the status of the grid to a file."""

        if not filename:
            filename = f"{self.data.name}_info.txt"

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

        empty_topo = np.mean(self.topo()[self.land_sea_mask()]) == 9999

        msg.header(self, f"Status of grid {self.data.name}")
        msg.plain(f'lon: {self.data.lon_min} - {self.data.lon_max}, lat: {self.data.lat_min} - {self.data.lat_max}')

        if hasattr(self.data, 'dlon') and hasattr(self.data, 'dlon'):
            msg.plain(f'dlon, dlat = {self.data.dlon}, {self.data.dlat} deg')
            msg.plain(f'Inverse: dlon, dlat = 1/{1/self.data.dlon}, 1/{1/self.data.dlat} deg')

        if hasattr(self.data, 'dx') and hasattr(self.data, 'dy'):
            msg.plain(f'dx, dy approximately {self.data.dx}, {self.data.dy} metres')

        if hasattr(self.data, 'nx') and hasattr(self.data, 'ny'):
            msg.plain(f'nx, ny = {self.data.nx} x {self.data.ny} grid points')

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
