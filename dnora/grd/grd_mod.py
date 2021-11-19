from abc import ABC, abstractmethod
import xarray as xr
import numpy as np
from copy import copy
from scipy.interpolate import griddata
import sys
import re
from .. import msg
from ..aux import distance_2points, add_prefix, add_suffix, add_file_extension, create_filename_obj, add_folder_to_filename
from ..defaults import dflt_grd

#from ..bnd.bnd_mod import Boundary # Boundary object
## -----------------------------------------------------------------------------
## READING THE DATA
## -----------------------------------------------------------------------------
class TopoReader(ABC):
    """Abstract class for reading the bathymetry."""
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float):
        """Reads the bathymetrical information from a source and returns the data.

        This method is called from within the Grid-object
        """
        pass

    @abstractmethod
    def __str__(self):
        """Describes what topography is read and from where.

        This is called by the Grid-objeect to provide output to the user.
        """
        pass


class EmptyTopo(TopoReader):
    """Creates an empty topography. Called when setting initial spacing."""
    def __init__(self, grid):
        self.grid = grid
        pass

    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float):
        """Creates a trivial topography with all water points."""
        topo = np.ones((self.grid.data.ny,self.grid.data.nx))*-9999
        topo_lon = copy(self.grid.lon())
        topo_lat = copy(self.grid.lat())
        return topo, topo_lon, topo_lat

    def __str__(self):
        return("Creating an empty topography with depth values -9999.")
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## DEFINING BOUNDARY POINTS
## -----------------------------------------------------------------------------
class BoundarySetter(ABC):
    """Abstract class for defining metods for setting boundary points in the
    grid."""
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, mask_size: tuple):
        """This method is called from within the Grid-object
        """
        return boundary_mask

    @abstractmethod
    def __str__(self):
        """Describes how the boundary points are set.

        This is called by the Grid-objeect to provide output to the user.
        """
        pass


class ClearBoundary(BoundarySetter):
    """Clears all boundary points by setting a mask with False values."""
    def __init__(self):
        pass

    def __call__(self, mask_size: tuple):
        return np.full(mask_size, False)

    def __str__(self):
        return("Clearing all possible boundary points and setting an empty mask.")
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## PROCESSING THE GRID DATA
## -----------------------------------------------------------------------------
class GridProcessor(ABC):
    """Abstract class for modifying bathymetrical data of the Grid-object."""
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, data, lon, lat, land_sea_mask, boundary_mask):
        """Gets the bathymetrical information and returns a modified version.

        This method is called from within the Grid-object
        """
        pass

    @abstractmethod
    def __str__(self):
        """
        Describes how the data is processed.

        This is called by the Grid-objeect to provide output to the user.
        """
        pass


class TrivialFilter(GridProcessor):
    """Returns the identical data it is passed. Used as default option."""
    def __init__(self):
        pass

    def __call__(self, data, lon, lat, land_sea_mask, boundary_mask):
        return copy(data)

    def __str__(self):
        return("Doing nothing to the data, just passing it along.")
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## MESHING THE DATA DOWN TO THE GRID
## -----------------------------------------------------------------------------
class Mesher(ABC):
    """Abstract class for meshing the bathymetrical data to the grid."""
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, data, lon, lat, lonQ, latQ):
        """Gets the bathymetrical information and returns a version that is
        meshed to the area and resolution of the grid.

        This method is called from within the Grid-object
        """
        pass

    @abstractmethod
    def __str__(self):
        """Describes how the data is meshed.

        This is called by the Grid-objeect to provide output to the user.
        """
        pass


class Interpolate(Mesher):
    """Interpolates data to grid. A wrapper for scipy.interpolate's griddate."""
    def __init__(self, method = 'linear'):
        self.method = method

        return

    def __call__(self, data, lon, lat, lonQ, latQ):
        lon0, lat0 = np.meshgrid(lon, lat)
        lon1, lat1 = np.meshgrid(lonQ, latQ)
        data[np.isnan(data)] = 0 # Keeping land points as nan lets the shoreline creep out
        M = np.column_stack((data.ravel(), lon0.ravel(),lat0.ravel()))
        meshed_data = griddata(M[:,1:], M[:,0], (lon1, lat1), method=self.method)
        meshed_data[meshed_data>=0] = 32767

        return meshed_data

    def __str__(self):
        return(f"Meshing using {self.method} interpolation.")


class TrivialMesher(Mesher):
    """Passes along data. NB! This might not fit the grid, and is only used
    for e.g. recreating a Grid-object from an ouput file.
    """
    def __init__(self):
        pass

    def __call__(self, data, lon, lat, lonQ, latQ):
        return copy(data)

    def __str__(self):
        return("Passing input data along as final meshed grid.")
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## THE GRID OBJECT
## -----------------------------------------------------------------------------
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

        msg.header(f'Importing topography with {type(topo_reader).__name__}')
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
            msg.header(f'Filtering topography with {type(filt).__name__}')

            empty_mask = np.full(self.raw_topo().shape, False)
            land_sea_mask = self.raw_topo() < 0 # Sea points set to true

            print(filt)
            topo = filt(self.raw_topo(), self.raw_lon(), self.raw_lat(), land_sea_mask, empty_mask)

            vars_dict = {'topo': (['lat', 'lon'], topo)}
            self.rawdata = self.rawdata.assign(vars_dict)

            return

    def mesh_grid(self, mesher: Mesher = Interpolate(method = 'linear')) -> None:
        """Meshes the raw data down to the grid definitions."""

        if hasattr(self, 'lon') and hasattr(self, 'lon'):

            msg.header(f'Meshing grid bathymetry with {type(mesher).__name__}')
            print(mesher)
            topo = mesher(self.raw_topo(), self.raw_lon(), self.raw_lat(), self.lon(), self.lat())
            vars_dict = {'topo': (['lat', 'lon'], topo)}
            self.data = self.data.assign(vars_dict)

            self._update_masks()
            return

        else:
            msg.templates('no_spacing')
            return

    def process_grid(self, filt: GridProcessor = TrivialFilter()) -> None:
        """Processes the gridded bathymetrical data, e.g. with a filter."""

        msg.header(f'Filtering meshed grid with {type(filt).__name__}')
        print(filt)
        topo = filt(self.topo(), self.lon(), self.lat(), self.land_sea_mask(), self.boundary_mask())

        vars_dict = {'topo': (['lat', 'lon'], topo)}
        self.data = self.data.assign(vars_dict)

        msg.info('Upodating land-sea mask and boundary mask')
        self._update_masks()

        return

    def export_grid(self, grid_writer):
        output_file, output_folder = grid_writer(self)

        # This is set as info in case an input file needs to be generated
        self._written_as = output_file
        self._written_to = output_folder
        return

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

        msg.header("Setting grid spacing")

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

        print(self)

        return

    def set_boundary(self, boundary_setter: BoundarySetter) -> None:
        """Marks the points that should be treated as boundary points in the
        grid."""

        msg.header(f'Setting boundary points with {type(boundary_setter).__name__}')
        print(boundary_setter)

        boundary_mask = boundary_setter(self.land_sea_mask().shape)

        vars_dict = {'boundary_mask': (['lat', 'lon'], boundary_mask)}
        self.data = self.data.assign(vars_dict)

        return

    def land_sea_mask(self):
        """Returns bool array of the land-sea mask (True=sea point)"""
        if hasattr(self.data, 'land_sea_mask'):
            return copy(self.data.land_sea_mask.values)
        else:
            return np.array([])

    def boundary_mask(self):
        """Returns bool array of boundary points (True=boundary point)"""
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
        """Returns a lon, lat list of land points."""
        if self.boundary_mask().size > 0:
            mask = self.land_sea_mask()
            LAND = self._point_list(mask)
            return LAND
        else:
            return np.array([])

    def filename(self, filestring: str=dflt_grd['fs']['General'], extension: str='', prefix: str='', suffix: str=''):
        # Substitute placeholders for $Grid
        filename = create_filename_obj(filestring=filestring, objects=[self])

        filename = add_prefix(filename=filename, prefix=prefix)
        filename = add_suffix(filename=filename, suffix=suffix)

        # Possible clean up
        filename = re.sub(f"__", '_', filename)
        filename = re.sub(f"_$", '', filename)

        if extension:
            filename = add_file_extension(filename, extension=extension)

        return filename

    def written_as(self, filestring: str=dflt_grd['fs']['General']):
        if hasattr(self, '_written_as'):
            return self._written_as
        else:
            return self.filename(filestring=filestring)

    def written_to(self):
        if hasattr(self, '_written_to'):
            return self._written_to
        else:
            return ''

    def is_written(self):
        return hasattr(self, '_written_as')


    def name(self):
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

        self._set_land_sea_mask()

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

    def _set_land_sea_mask(self, matrix = None) -> None:
        """Sets the land-sea mask based on land points in the gridded
        bathymetrical data.

        A certain land-sea mask can also be forced by providing a boolean
        matrix by the user (as a bit of hack). Note, that this will be
        overwritten automatically by the method _update_masks() if any further
        changes are made to the topography.
        """

        if matrix is None:
            land_sea_mask = np.full(self.topo().shape, False)
            land_sea_mask = self.topo() < 0 # Sea points set to true
        else:
            land_sea_mask = matrix

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
            msg.header(f"Status of grid {self.data.name}")
            msg.plain(f'lon: {self.data.lon_min} - {self.data.lon_max}, lat: {self.data.lat_min} - {self.data.lat_max}')
            if hasattr(self.data, 'dlon') and hasattr(self.data, 'dlon'):
                msg.plain(f'dlon, dlat = {self.data.dlon}, {self.data.dlat} deg')
                msg.plain(f'Inverse: dlon, dlat = 1/{1/self.data.dlon}, 1/{1/self.data.dlat} deg')
            if hasattr(self.data, 'dx') and hasattr(self.data, 'dy'):
                msg.plain(f'dx, dy approximately {self.data.dx}, {self.data.dy} metres')
            if hasattr(self.data, 'nx') and hasattr(self.data, 'ny'):
                msg.plain(f'nx, ny = {self.data.nx} x {self.data.ny} grid points')
            if self.topo().size > 0:
                msg.plain(f"Mean depth: {np.mean(self.topo()[self.land_sea_mask()]):.1f} m")
                msg.plain(f"Max depth: {np.min(self.topo()[self.land_sea_mask()]):.1f} m")
                msg.plain(f"Min depth: {np.max(self.topo()[self.land_sea_mask()]):.1f} m")
            if self.land_sea_mask().size > 0:
                msg.print_line()
                msg.plain('Grid contains:')
                msg.plain(f'{sum(sum(self.land_sea_mask())):d} sea points')
                msg.plain(f'{sum(sum(np.logical_not(self.land_sea_mask()))):d} land points')
            if self.boundary_mask().size > 0:
                msg.plain(f'{sum(sum(np.logical_and(self.boundary_mask(), self.land_sea_mask()))):d} boundary points')
            msg.print_line()
            return ''
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## WRITING THE GRID DATA TO FILES
## -----------------------------------------------------------------------------
class GridWriter(ABC):
    """Abstract class for writing the Grid-object's data to files to be Used
    by the wave models.
    """
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, grid: Grid):
        pass
## -----------------------------------------------------------------------------
