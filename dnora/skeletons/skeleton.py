import numpy as np
import pandas as pd
import xarray as xr
import utm
from copy import copy
from .dataset_manager import DatasetManager
from .. import aux_funcs

class Skeleton:
    """Contains methods and data of the spatial x,y / lon, lat coordinates and
    makes possible conversions between them.

    Keeps track of the native structure of the grid (cartesian UTM / sperical).
    """
    #_name = 'LonelySkeleton'
    #strict = False # If this is True, no coordinate conversions will be done (return None instead)

    _xy_dict = {'x': 'x', 'lon': 'x', 'y': 'y', 'lat': 'y'}


    def _init_structure(self, x=None, y=None, lon=None, lat=None, **kwargs) -> None:
        """Determines grid type (Cartesian/Spherical), generates a DatasetManager
        and initializes the Xarray dataset within the DatasetManager.

        The initial coordinates and variables are read from the method of the
        subclass (e.g. PointSkeleton)
        """

        # Migth have been already created by decorators
        if not hasattr(self, '_coord_manager'):
            self._coord_manager =  CoordinateManager()

        #if not self.is_initialized():
        x_str, y_str, xvec, yvec = will_grid_be_spherical_or_cartesian(x, y, lon, lat)

        self.x_str = x_str
        self.y_str = y_str

        # Initial values defined in subclass (e.g. GriddedSkeleton)
        self._coord_manager.set_initial_coords(self._initial_coords())
        self._coord_manager.set_initial_vars(self._initial_vars())
        # else:
        #     xvec = x
        #     yvec = y

        if not hasattr(self, 'ds_manager'):
            self.ds_manager = DatasetManager(self._coord_manager)

        # The manager contains the Xarray Dataset
        self.ds_manager.create_structure(xvec,yvec,self.x_str,self.y_str,**kwargs)
        # if self.ds() is None or not append:
        #     self.ds_manager.set_new_ds(ds)
        # else:
        #     self.ds_manager.set_new_ds(xr.concat([self.ds(), ds], dim="time").sortby('time'))
        self._reset_masks()
        self._reset_datavars()
        self._structure_initialized = True

    def _absorb_object(self, obj, dimension: str) -> None:
        """Absorb another object of same type. This is used e.g. when pathcing
        cached data and joining different Boundary etc. over time.
        """
        self.ds_manager.set_new_ds(xr.concat([self.ds(), obj.ds()], dim=dimension).sortby(dimension))

    def _reset_masks(self) -> None:
        """Resets the mask to default values."""
        for name in self._coord_manager.added_masks():
            # update-method sets empty mask when no is provided
            self._update_mask(name)

    def _reset_datavars(self) -> None:
        """Resets the data variables to default values."""
        for name in self._coord_manager.added_vars():
            # update-method sets empty mask when no is provided
            self._update_datavar(name)

    def _update_mask(self, name: str, updated_mask=None) -> None:
        coords = self._coord_manager.added_masks().get(name)
        if name is None:
            raise ValueError(f'A mask named {name} has not been defines ({list(masks.keys())})')

        if updated_mask is None:
            updated_mask = self.get(f'{name}_mask',empty=True)
        self.ds_manager.set(data=updated_mask.astype(int), data_name=f'{name}_mask', coord_type=coords)

    def _update_datavar(self, name: str, updated_var=None) -> None:
        coords = self._coord_manager.added_vars().get(name)
        if name is None:
            raise ValueError(f'A data variable named {name} has not been defines ({list(vars.keys())})')

        if updated_var is None:
            updated_var = self.get(name, empty=True)
        self.ds_manager.set(data=updated_var, data_name=name, coord_type=coords)

    def get(self, name, empty=False):
        """Gets a mask or data variable.

        The ds_manager always gets what is in the Dataset (integers for masks).
        The Skeletons get-method gives boolen masks, and you can also
        request empty masks that willb e return even if data doesn't exist."""
        if empty:
            return eval(f'self.{name}(empty=True)')

        data = self.ds_manager.get(name)
        if data is None:
            return None

        data = eval(f'self.{name}()')
        return data

    def is_empty(self, name):
        """Checks if a Dataset variable is empty.

        Empty means all initial values OR all 0 values."""
        data = self.get(name)
        if data is None:
            return False
        empty_data = self.get(name, empty=True)
        land_data = data*0
        is_empty = np.allclose(data.astype(float), empty_data.astype(float)) or np.allclose(data.astype(float), land_data.astype(float))
        return is_empty


    def is_initialized(self) -> bool:
        return hasattr(self, 'x_str') and hasattr(self, 'y_str')

    def is_cartesian(self) -> bool:
        """Checks if the grid is cartesian (True) or spherical (False)."""
        if self.x_str == 'x' and self.y_str == 'y':
            return True
        elif self.x_str == 'lon' and self.y_str == 'lat':
            return False
        raise Exception(f"Expected x- and y string to be either 'x' and 'y' or 'lon' and 'lat', but they were {x_str} and {y_str}")

    def is_gridded(self) -> bool:
        return aux_funcs.is_gridded(self.topo(), self.native_x(), self.native_y())

    def set_utm(self, zone_number: int=33, zone_letter: str='W'):
        """Set UTM zone and number to be used for cartesian coordinates."""
        self._zone_number = copy(zone_number)
        self._zone_letter = copy(zone_letter)

    def ds(self):
        return self.ds_manager.ds()

    def size(self, type: str='all') -> tuple[int]:
        """Returns the size of the Dataset.

        'all': size of entire Dataset
        'spatial': size over coordinates from the Skeleton (x, y, lon, lat, inds)
        'grid': size over coordinates for the grid (e.g. z, time)
        'gridpoint': size over coordinates for a grid point (e.g. frequency, direcion or time)
        """
        return self.ds_manager.size(type)

    def inds(self, **kwargs) -> np.ndarray:
        return self.ds_manager.get('inds', **kwargs).values.copy()

    def native_x(self, **kwargs) -> np.ndarray:
        """Returns x-vector for cartesian grids and lon-vector for sperical
        grids.
        """
        if self.is_cartesian():
            return self.x(**kwargs)
        return self.lon(**kwargs)

    def native_y(self, **kwargs) -> np.ndarray:
        """Returns y-vector for cartesian grids and lat-vector for sperical
        grids.
        """
        if self.is_cartesian():
            return self.y(**kwargs)
        return self.lat(**kwargs)


    def x(self, strict=False, **kwargs) -> np.ndarray:
        """Returns the cartesian x-coordinate. If the grid is spherical a
        conversion to UTM coordinates is made based on the medain latitude.

        If strict=True, then None is returned if grid is sperical.
        """
        if not self.is_cartesian() and (strict or self.strict):
            return None

        if not self.is_cartesian():
            lat = np.median(self.lat(**kwargs))
            number, letter = self.utm()
            x, __, __, __ = utm.from_latlon(lat, self.lon(**kwargs), force_zone_number=number, force_zone_letter=letter)
            return x

        return self.ds_manager.get('x', **kwargs).values.copy()

    def y(self, strict=False, **kwargs) -> np.ndarray:
        """Returns the cartesian y-coordinate. If the grid is spherical a
        conversion to UTM coordinates is made based on the medain longitude.

        If strict=True, then None is returned if grid is sperical.
        """
        if not self.is_cartesian() and (strict or self.strict):
            return None

        if not self.is_cartesian():
            lon = np.median(self.lon(**kwargs))
            number, letter = self.utm()
            __, y, __, __ = utm.from_latlon(self.lat(**kwargs), lon, force_zone_number=number, force_zone_letter=letter)
            return y

        return self.ds_manager.get('y', **kwargs).values.copy()

    def lon(self, strict=False, **kwargs) -> np.ndarray:
        """Returns the spherical lon-coordinate. If the grid is cartesian (UTM)
        a conversion to spherical coordinates is made based on the medain
        y-values.

        If strict=True, then None is returned if grid is cartesian.
        """
        if self.is_cartesian() and (strict or self.strict):
            return None

        if self.is_cartesian():
            if self.is_gridded(): # This will rotate the grid, but is best estimate to keep it strucutred
                y = np.median(self.y(**kwargs))
            else:
                y = self.y(**kwargs)
            number, letter = self.utm()
            __, lon = utm.to_latlon(self.x(**kwargs), y, zone_number=number, zone_letter=letter, strict = False)

            return lon
        return self.ds_manager.get('lon', **kwargs).values.copy()

    def lat(self, strict=False, **kwargs) -> np.ndarray:
        """Returns the spherical at-coordinate. If the grid is cartesian (UTM)
        a conversion to spherical coordinates is made based on the medain
        x-values.

        If strict=True, then None is returned if grid is cartesian.
        """
        if self.is_cartesian() and (strict or self.strict):
            return None

        if self.is_cartesian():
            if self.is_gridded(): # This will rotate the grid, but is best estimate to keep it strucutred
                x = np.median(self.x(**kwargs))
            else:
                x = self.x(**kwargs)
            number, letter = self.utm()
            lat, __ = utm.to_latlon(x, self.y(**kwargs), zone_number=number, zone_letter=letter, strict = False)
            return lat

        return self.ds_manager.get('lat', **kwargs).values.copy()

    def _xy(self, x: np.ndarray=None, y: np.ndarray=None, strict=False) -> tuple[np.ndarray, np.ndarray]:
        """Converts list of points to x and y (UTM) if necessary.

        Input is assumed to be in the native format (UTM for cartesian, WGS84
        for spherical).

        If strict=True, then (None, None) is returned if grid is sperical.
        """
        if not self.is_cartesian() and (strict or self.strict):
            return None, None

        if not self.is_cartesian():
            number, letter = self.utm()
            x, y, __, __ = utm.from_latlon(y, x, force_zone_number=number, force_zone_letter=letter)
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
            number, letter = self.utm()
            lat, lon = utm.to_latlon(lon, lat, zone_number=number, zone_letter=letter, strict = False)
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

        edges = self.edges('x')
        return (edges[1]-edges[0])/(self.nx()-1)

    def dy(self, strict=False):
        """Mean grid spacing of the y vector. Conversion made for
        spherical grids."""
        if not self.is_cartesian() and (strict or self.strict):
            return None

        if self.ny() == 1:
            return 0.

        edges = self.edges('y')
        return (edges[1]-edges[0])/(self.ny()-1)

    def dlon(self, strict=False):
        """Mean grid spacing of the longitude vector. Conversion made for
        cartesian grids."""
        if self.is_cartesian() and (strict or self.strict):
            return None
        if self.nx() == 1:
            return 0.

        edges = self.edges('lon')
        return (edges[1]-edges[0])/(self.nx()-1)

    def dlat(self, strict=False):
        """Mean grid spacing of the latitude vector. Conversion made for
        cartesian grids."""
        if self.is_cartesian() and (strict or self.strict):
            return None
        if self.ny() == 1:
            return 0.

        edges = self.edges('lat')
        return (edges[1]-edges[0])/(self.ny()-1)

    def native_dx(self):
        """Mean grid spacing of x vector for cartesian grids.
        Mean grid spacing of lon vector for spherical grids."""
        if self.nx() == 1:
            return 0.

        edges = self.edges('x', native=True)
        return (edges[1]-edges[0])/(self.nx()-1)

    def native_dy(self):
        """Mean grid spacing of y vector for cartesian grids.
        Mean grid spacing of lat vector for spherical grids."""
        if self.ny() == 1:
            return 0.

        edges = self.edges('y', native=True)
        return (edges[1]-edges[0])/(self.ny()-1)

    def utm(self) -> tuple[int, str]:
        """Returns UTM zone number and letter. Returns 33, 'W' as default
        value if it hasn't been set by the user."""
        number, letter = 33, 'W'
        if hasattr(self, '_zone_number'):
            number = self._zone_number
        if hasattr(self, '_zone_letter'):
            letter = self._zone_letter
        return number, letter

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
    def strict(self) -> bool:
        """If this is set to true, then no coordinate conversion will ever be
        made when requesting lon, lat, x, y etc."""
        if not hasattr(self, '_strict'):
            return False
        return self._strict

    @strict.setter
    def name(self, strict: bool) -> None:
        if isinstance(strict, bool):
            self._strict = strict
        else:
            raise ValueError("strict needs to be a bool")

    @property
    def name(self) -> str:
        if not hasattr(self,'_name'):
            return 'LonelySkeleton'
        return self._name

    @name.setter
    def name(self, new_name):
        if isinstance(new_name, str):
            self._name = new_name
        else:
            raise ValueError("name needs to be a string")


    def size(self, type: str='all', **kwargs) -> tuple[int]:
        """Returns the size of the Dataset.

        'all': size of entire Dataset
        'spatial': size over coordinates from the Skeleton (x, y, lon, lat, inds)
        'grid': size over coordinates for the grid (e.g. z, time)
        'gridpoint': size over coordinates for a grid point (e.g. frequency, direcion or time)
        """
        return self.ds_manager.coords_to_size(self.ds_manager.coords(type), **kwargs)

def will_grid_be_spherical_or_cartesian(x, y, lon, lat):
    """Determines if the grid will be spherical or cartesian based on which
    inputs are given and which are None.

    Returns the ringth vector and string to identify the native values.
    """
    xy = False
    lonlat = False
    if isinstance(x,tuple) and np.all(x==(None, None)):
        x = None
    if isinstance(y,tuple) and np.all(y==(None, None)):
        y = None
    if isinstance(lon,tuple) and np.all(lon==(None, None)):
        lon = None
    if isinstance(lat,tuple) and np.all(lat==(None, None)):
        lat = None
#    none_tuple = (None, None)
#    if (x is not None and x != none_tuple) and (y is not None and y != none_tuple):
    if (x is not None) and (y is not None):
        xy = True
        native_x = 'x'
        native_y = 'y'
        xvec = x
        yvec = y
        xvec = np.array(x)
        yvec = np.array(y)

    #if (lon is not None and lon != none_tuple) and (lat is not None and lat != none_tuple):
    if (lon is not None) and (lat is not None):
        lonlat = True
        native_x = 'lon'
        native_y = 'lat'
        xvec = np.array(lon)
        yvec = np.array(lat)

    if xy and lonlat:
        raise ValueError("Can't set both lon/lat and x/y!")

    if not xy and not lonlat:
        raise ValueError('Have to set either lon/lat or x/y!')

    if not xvec.shape:
        xvec = np.array([xvec])

    if not yvec.shape:
        yvec = np.array([yvec])

    if len(np.unique(xvec)) == 1 and len(xvec) == 2: # e.g. lon=(4.0, 4.0) should behave like lon=4.0
        xvec = np.unique(xvec)

    if len(np.unique(yvec)) == 1 and len(yvec) == 2:
        yvec = np.unique(yvec)

    return native_x, native_y, xvec, yvec
