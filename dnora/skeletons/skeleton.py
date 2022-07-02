import numpy as np
import pandas as pd
import xarray as xr
import utm
from copy import copy

class Skeleton:
    """Contains methods and data of the spatial x,y / lon, lat coordinates and
    makes possible conversions between them.

    Keeps track of the native structure of the grid (cartesian UTM / sperical).
    """
    _x_str = 'x'
    _y_str = 'y'
    _name = 'LonelySkeleton'
    _zone_number = 33
    _zone_letter = 'W'
    strict = False # If this is True, no coordinate conversions will be done (return None instead)

    _xy_dict = {'x': 'x', 'lon': 'x', 'y': 'y', 'lat': 'y'}
    _cartesian_strings = ['x','y','xy']
    _spherical_strings = ['lon','lat','lonlat']

    def is_cartesian(self) -> bool:
        """Checks if the grid is cartesian (True) or spherical (False)."""
        if self.x_str == 'x' and self.y_str == 'y':
            return True
        elif self.x_str == 'lon' and self.y_str == 'lat':
            return False
        raise Exception(f"Expected x- and y string to be either 'x' and 'y' or 'lon' and 'lat', but they were {x_str} and {y_str}")

    def set_utm(self, zone_number: int=33, zone_letter: str='W'):
        """Set UTM zone and number to be used for cartesian coordinates."""
        self._zone_number = copy(zone_number)
        self._zone_letter = copy(zone_letter)

    def inds(self) -> np.ndarray:
        return self._get('inds')

    def native_x(self) -> np.ndarray:
        """Returns x-vector for cartesian grids and lon-vector for sperical
        grids.
        """
        if self.is_cartesian():
            return self.x()
        return self.lon()

    def native_y(self) -> np.ndarray:
        """Returns y-vector for cartesian grids and lat-vector for sperical
        grids.
        """
        if self.is_cartesian():
            return self.y()
        return self.lat()


    def x(self, strict=False) -> np.ndarray:
        """Returns the cartesian x-coordinate. If the grid is spherical a
        conversion to UTM coordinates is made based on the medain latitude.

        If strict=True, then None is returned if grid is sperical.
        """
        if not self.is_cartesian() and (strict or self.strict):
            return None

        if not self.is_cartesian():
            lat = np.median(self.lat())
            x, __, __, __ = utm.from_latlon(lat, self.lon(), force_zone_number=self._zone_number, force_zone_letter=self._zone_letter)
            return x

        return self._get('x')

    def y(self, strict=False) -> np.ndarray:
        """Returns the cartesian y-coordinate. If the grid is spherical a
        conversion to UTM coordinates is made based on the medain longitude.

        If strict=True, then None is returned if grid is sperical.
        """
        if not self.is_cartesian() and (strict or self.strict):
            return None

        if not self.is_cartesian():
            lon = np.median(self.lon())
            __, y, __, __ = utm.from_latlon(self.lat(), lon, force_zone_number=self._zone_number, force_zone_letter=self._zone_letter)
            return y

        return self._get('y')

    def lon(self, strict=False) -> np.ndarray:
        """Returns the spherical lon-coordinate. If the grid is cartesian (UTM)
        a conversion to spherical coordinates is made based on the medain
        y-values.

        If strict=True, then None is returned if grid is cartesian.
        """
        if self.is_cartesian() and (strict or self.strict):
            return None

        if self.is_cartesian():
            y = np.median(self.y())
            __, lon = utm.to_latlon(self.x(), y, self._zone_number, zone_letter=self._zone_letter, strict = False)
            return lon

        return self._get('lon')

    def lat(self, strict=False) -> np.ndarray:
        """Returns the spherical at-coordinate. If the grid is cartesian (UTM)
        a conversion to spherical coordinates is made based on the medain
        x-values.

        If strict=True, then None is returned if grid is cartesian.
        """
        if self.is_cartesian() and (strict or self.strict):
            return None

        if self.is_cartesian():
            x = np.median(self.x())
            lat, __ = utm.to_latlon(x, self.y(), self._zone_number, zone_letter=self._zone_letter, strict = False)
            return lat

        return self._get('lat')

    def _xy(self, x: np.ndarray=None, y: np.ndarray=None, strict=False) -> tuple[np.ndarray, np.ndarray]:
        """Converts list of points to x and y (UTM) if necessary.

        Input is assumed to be in the native format (UTM for cartesian, WGS84
        for spherical).

        If strict=True, then (None, None) is returned if grid is sperical.
        """
        if not self.is_cartesian() and (strict or self.strict):
            return None, None

        if not self.is_cartesian():
            x, y, __, __ = utm.from_latlon(y, x, force_zone_number=self._zone_number, force_zone_letter=self._zone_letter)
        return x, y

    def _lonlat(self, lon: np.ndarray=None, lat: np.ndarray=None, strict=False) -> tuple[np.ndarray, np.ndarray]:
        """Converts list of points to longitude and latitude if necessary.

        Input is assumed to be in the native format (UTM for cartesian, WGS84
        for spherical).

        If strict=True, then (None, None) is returned if grid is cartesian.
        """
        if self.is_cartesian() and (strict or self.strict):
            return None, None

        if self.is_cartesian():
            lat, lon = utm.to_latlon(lon, lat, self._zone_number, zone_letter=self._zone_letter, strict = False)
        return lon, lat

    def edges(self, coord: str, native: bool=False, strict=False) -> tuple[float, float]:
        """Min and max values of x. Conversion made for sperical grids."""
        if self._xy_dict.get(coord) is None:
            print("coord need to be 'x', 'y', 'lon' or 'lat'.")
            return

        if native:
            method = f'native_{self._xy_dict[coord]}'
            val = getattr(self, method)()
        else:
            val = getattr(self, coord)(strict=(strict or self.strict))

        if val is None:
            return (None, None)

        return np.min(val), np.max(val)

    def nx(self) -> int:
        """Length of x/lon-vector."""
        return len(self.native_x())

    def ny(self):
        """Length of y/lat-vector."""
        return len(self.native_y())

    def dx(self, strict=False):
        """Mean grid spacing of the x vector. Conversion made for
        spherical grids."""
        if not self.is_cartesian() and (strict or self.strict):
            return None

        if self.nx() == 1:
            return 0.

        edges = self.edges(coord='x')
        return (edges[1]-edges[0])/(self.nx()-1)

    def dy(self, strict=False):
        """Mean grid spacing of the y vector. Conversion made for
        spherical grids."""
        if not self.is_cartesian() and (strict or self.strict):
            return None

        if self.ny() == 1:
            return 0.

        edges = self.edges(coord='y')
        return (edges[1]-edges[0])/(self.ny()-1)

    def dlon(self, strict=False):
        """Mean grid spacing of the longitude vector. Conversion made for
        cartesian grids."""
        if self.is_cartesian() and (strict or self.strict):
            return None
        if self.nx() == 1:
            return 0.

        edges = self.edges(coord='dlon')
        return (edges[1]-edges[0])/(self.nx()-1)

    def dlat(self):
        """Mean grid spacing of the latitude vector. Conversion made for
        cartesian grids."""
        if self.is_cartesian() and (strict or self.strict):
            return None
        if self.ny() == 1:
            return 0.

        edges = self.edges(coord='dlat')
        return (edges[1]-edges[0])/(self.ny()-1)

    def native_dx(self):
        """Mean grid spacing of x vector for cartesian grids.
        Mean grid spacing of lon vector for spherical grids."""
        if self.nx() == 1:
            return 0.

        edges = self.edges(coord='x', native=True)
        return (edges[1]-edges[0])/(self.nx()-1)

    def native_dy(self):
        """Mean grid spacing of y vector for cartesian grids.
        Mean grid spacing of lat vector for spherical grids."""
        if self.ny() == 1:
            return 0.

        edges = self.edges(coord='y', native=True)
        return (edges[1]-edges[0])/(self.ny()-1)

    @property
    def x_str(self) -> str:
        """Return string compatible with the type of spacing used:

        'x' for cartesian grid.
        'lon' for spherical grid.
        """
        return self._x_str

    @x_str.setter
    def x_str(self, new_str):
        if new_str in ['x', 'lon']:
            self._x_str = new_str
        else:
            raise ValueError("x_str need to be 'x' or 'lon'")

    @property
    def y_str(self) -> str:
        """Return string compatible with the type of spacing used:

        'y' for cartesian grid.
        'lat' for spherical grid.
        """
        return self._y_str

    @y_str.setter
    def y_str(self, new_str):
        if new_str in ['y', 'lat']:
            self._y_str = new_str
        else:
            raise ValueError("y_str need to be 'y' or 'lat'")

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, new_name):
        if isinstance(new_name, str):
            self._name = new_name
        else:
            raise ValueError("name needs to be a string")
