import numpy as np
import pandas as pd
import xarray as xr
import utm
from copy import copy
from .dataset_manager import DatasetManager
from .. import aux_funcs
from .. import msg

def cap_lat_for_utm(lat):
    if isinstance(lat, float):
        lat = np.array([lat])
    if max(lat) > 84:
        msg.warning(f'Max latitude {max(lat)}>84. These points well be capped to 84 deg in UTM conversion!')
        lat[lat>84.] = 84.
    if min(lat) < -80:
        lat[lat<-80.] = -80.
        msg.warning(f'Min latitude {min(lat)}<-80. These points well be capped to -80 deg in UTM conversion!')
    return lat

class Skeleton:
    """Contains methods and data of the spatial x,y / lon, lat coordinates and
    makes possible conversions between them.

    Keeps track of the native structure of the grid (cartesian UTM / sperical).
    """
    #_name = 'LonelySkeleton'
    #strict = False # If this is True, no coordinate conversions will be done (return None instead)

    _xy_dict = {'x': 'x', 'lon': 'x', 'y': 'y', 'lat': 'y'}
    _x_to_xy_dict = {'x': 'xy', 'y': 'xy', 'lon': 'lonlat', 'lat': 'lonlat'}

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
        self.set_utm()
        # if self.ds() is None or not append:
        #     self.ds_manager.set_new_ds(ds)
        # else:
        #     self.ds_manager.set_new_ds(xr.concat([self.ds(), ds], dim="time").sortby('time'))
        self._reset_masks()
        self._reset_datavars()

    def _structure_initialized(self) -> bool:
        return hasattr(self, 'ds_manager')

    def _absorb_object(self, obj, dimension: str) -> None:
        """Absorb another object of same type. This is used e.g. when pathcing
        cached data and joining different Boundary etc. over time.
        """
        self.ds_manager.set_new_ds(xr.concat([self.ds(), obj.ds()], dim=dimension, data_vars='minimal').sortby(dimension))

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
        if not self._structure_initialized():
            return None

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
        if not self._structure_initialized():
            return False
        if self.x_str == 'x' and self.y_str == 'y':
            return True
        elif self.x_str == 'lon' and self.y_str == 'lat':
            return False
        raise Exception(f"Expected x- and y string to be either 'x' and 'y' or 'lon' and 'lat', but they were {x_str} and {y_str}")

    # def is_gridded(self) -> bool:
    #     if not self._structure_initialized():
    #         return False
    #     return aux_funcs.is_gridded(self.topo(), self.native_x(), self.native_y())

    def set_utm(self, zone_number: int=None, zone_letter: str=None):
        """Set UTM zone and number to be used for cartesian coordinates.

        If not given for a spherical grid, they will be deduced.

        If not given for a cartesian grid, will be set to 33, 'W'
        """

        if zone_number is None or zone_letter is None:
            if self.is_cartesian():
                zone_number = 33
                zone_letter = 'W'
            else:
                lon, lat = self.lonlat()
                # *** utm.error.OutOfRangeError: latitude out of range (must be between 80 deg S and 84 deg N)
                mask = np.logical_and(lat < 84, lat > -80)
                # raise OutOfRangeError('longitude out of range (must be between 180 deg W and 180 deg E)')

                lat, lon = lat[mask], lon[mask]

                # *** ValueError: latitudes must all have the same sign
                if len(lat[lat>=0]) > len(lat[lat<0]):
                    lat, lon = lat[lat>=0], lon[lat>=0]
                else:
                    lat, lon = lat[lat<0], lon[lat<0]

                __, __, zone_number, zone_letter = utm.from_latlon(lat, lon)

        if isinstance(zone_number, int) or isintance(zone_number, float):
            self._zone_number = copy(int(zone_number))
        else:
            raise ValueError("zone_number needs to be an integer")

        if isinstance(zone_letter, str):
            self._zone_letter = copy(zone_letter)
        else:
            raise ValueError("zone_number needs to be an integer")

        msg.header(self, f'Setting UTM ({zone_number}, {zone_letter})')

    def utm(self) -> tuple[int, str]:
        """Returns UTM zone number and letter. Returns 33, 'W' as default
        value if it hasn't been set by the user."""
        number, letter = 33, 'W'
        if hasattr(self, '_zone_number'):
            number = self._zone_number
        if hasattr(self, '_zone_letter'):
            letter = self._zone_letter
        return number, letter

    def ds(self):
        if not self._structure_initialized():
            return None
        return self.ds_manager.ds()

    def size(self, type: str='all') -> tuple[int]:
        """Returns the size of the Dataset.

        'all': size of entire Dataset
        'spatial': size over coordinates from the Skeleton (x, y, lon, lat, inds)
        'grid': size over coordinates for the grid (e.g. z, time)
        'gridpoint': size over coordinates for a grid point (e.g. frequency, direcion or time)
        """
        if not self._structure_initialized():
            return None
        return self.ds_manager.size(type)

    def inds(self, **kwargs) -> np.ndarray:
        if not self._structure_initialized():
            return None
        inds = self.ds_manager.get('inds', **kwargs)
        if inds is None:
            return None
        return inds.values.copy()

    def native_x(self, **kwargs) -> np.ndarray:
        """Returns x-vector for cartesian grids and lon-vector for sperical
        grids.
        """
        if not self._structure_initialized():
            return None
        if self.is_cartesian():
            return self.x(**kwargs)
        return self.lon(**kwargs)

    def native_y(self, **kwargs) -> np.ndarray:
        """Returns y-vector for cartesian grids and lat-vector for sperical
        grids.
        """
        if not self._structure_initialized():
            return None
        if self.is_cartesian():
            return self.y(**kwargs)
        return self.lat(**kwargs)


    def x(self, strict=False, **kwargs) -> np.ndarray:
        """Returns the cartesian x-coordinate. If the grid is spherical a
        conversion to UTM coordinates is made based on the medain latitude.

        If strict=True, then None is returned if grid is sperical.
        """

        if not self._structure_initialized():
            return None
        if not self.is_cartesian() and (strict or self.strict):
            return None

        if not self.is_cartesian():
            number, letter = self.utm()
            if self.is_gridded(): # This will rotate the grid, but is best estimate to keep it strucutred
                lat = np.median(self.lat(**kwargs))
                msg.warning('Regridding spherical grid to cartesian coordinates. This will cause a rotation!')
                x, __, __, __ = utm.from_latlon(lat, self.lon(**kwargs), force_zone_number=number, force_zone_letter=letter)
            else:
                lat = self.lat(**kwargs)
                lat = cap_lat_for_utm(lat)

                posmask = lat>=0
                negmask = lat<0
                x=np.zeros(len(lat))
                if np.any(posmask):
                    x[posmask], __, __, __ = utm.from_latlon(lat[posmask], self.lon(**kwargs)[posmask], force_zone_number=number, force_zone_letter=letter)
                if np.any(negmask):
                    x[negmask], __, __, __ = utm.from_latlon(-lat[negmask], self.lon(**kwargs)[negmask], force_zone_number=number, force_zone_letter=letter)
            #x, __, __, __ = utm.from_latlon(lat, self.lon(**kwargs), force_zone_number=number, force_zone_letter=letter)
            return x

        return self.ds_manager.get('x', **kwargs).values.copy()

    def y(self, strict=False, **kwargs) -> np.ndarray:
        """Returns the cartesian y-coordinate. If the grid is spherical a
        conversion to UTM coordinates is made based on the medain longitude.

        If strict=True, then None is returned if grid is sperical.
        """
        if not self._structure_initialized():
            return None
        if not self.is_cartesian() and (strict or self.strict):
            return None

        if not self.is_cartesian():
            number, letter = self.utm()
            posmask = self.lat(**kwargs)>=0
            negmask = self.lat(**kwargs)<0
            if self.is_gridded(): # This will rotate the grid, but is best estimate to keep it strucutred
                lon = np.median(self.lon(**kwargs))
                msg.warning('Regridding spherical grid to cartesian coordinates. This will cause a rotation!')
                y = np.zeros(len(self.lat(**kwargs)))
                if np.any(posmask):
                    _, y[posmask], __, __ = utm.from_latlon(self.lat(**kwargs)[posmask], lon, force_zone_number=number, force_zone_letter=letter)
                if np.any(negmask):
                    _, y[negmask], __, __ = utm.from_latlon(-self.lat(**kwargs)[negmask], lon, force_zone_number=number, force_zone_letter=letter)
                    y[negmask] = -y[negmask]
            else:
                lon = self.lon(**kwargs)
                lat = cap_lat_for_utm(self.lat(**kwargs))

                y = np.zeros(len(self.lat(**kwargs)))
                if np.any(posmask):
                    _, y[posmask], __, __ = utm.from_latlon(self.lat(**kwargs)[posmask], lon[posmask], force_zone_number=number, force_zone_letter=letter)
                if np.any(negmask):
                    _, y[negmask], __, __ = utm.from_latlon(-self.lat(**kwargs)[negmask], lon[negmask], force_zone_number=number, force_zone_letter=letter)
                    y[negmask] = -y[negmask]

            #
            return y

        return self.ds_manager.get('y', **kwargs).values.copy()

    def lon(self, strict=False, **kwargs) -> np.ndarray:
        """Returns the spherical lon-coordinate. If the grid is cartesian (UTM)
        a conversion to spherical coordinates is made based on the medain
        y-values.

        If strict=True, then None is returned if grid is cartesian.
        """
        if not self._structure_initialized():
            return None

        if self.is_cartesian() and (strict or self.strict):
            return None

        if self.is_cartesian():
            if self.is_gridded(): # This will rotate the grid, but is best estimate to keep it strucutred
                y = np.median(self.y(**kwargs))
                msg.warning('Regridding cartesian grid to spherical coordinates. This will cause a rotation!')
            else:
                y = self.y(**kwargs)
            number, letter = self.utm()
            __, lon = utm.to_latlon(self.x(**kwargs), np.mod(y, 10_000_000), zone_number=number, zone_letter=letter, strict = False)

            return lon
        return self.ds_manager.get('lon', **kwargs).values.copy()

    def lat(self, strict=False, **kwargs) -> np.ndarray:
        """Returns the spherical at-coordinate. If the grid is cartesian (UTM)
        a conversion to spherical coordinates is made based on the medain
        x-values.

        If strict=True, then None is returned if grid is cartesian.
        """
        if not self._structure_initialized():
            return None

        if self.is_cartesian() and (strict or self.strict):
            return None

        if self.is_cartesian():
            if self.is_gridded(): # This will rotate the grid, but is best estimate to keep it strucutred
                x = np.median(self.x(**kwargs))
                msg.warning('Regridding cartesian grid to spherical coordinates. This will cause a rotation!')
            else:
                x = self.x(**kwargs)
            number, letter = self.utm()
            lat, __ = utm.to_latlon(x, np.mod(self.y(**kwargs), 10_000_000), zone_number=number, zone_letter=letter, strict = False)
            return lat

        return self.ds_manager.get('lat', **kwargs).values.copy()

    def _xy(self, x: np.ndarray=None, y: np.ndarray=None, strict=False) -> tuple[np.ndarray, np.ndarray]:
        """Converts list of points to x and y (UTM) if necessary.

        Input is assumed to be in the native format (UTM for cartesian, WGS84
        for spherical).

        If strict=True, then (None, None) is returned if grid is sperical.
        """
        if not self._structure_initialized():
            return None

        if not self.is_cartesian() and (strict or self.strict):
            return None, None

        if not self.is_cartesian():
            y = cap_lat_for_utm(y)
            number, letter = self.utm()

            posmask = y>=0
            negmask = y<0

            if np.any(posmask):
                x[posmask], y[posmask], __, __ = utm.from_latlon(y[posmask], x[posmask], force_zone_number=number, force_zone_letter=letter)
            if np.any(negmask):
                x[negmask], y[negmask], __, __ = utm.from_latlon(-y[negmask], x[negmask], force_zone_number=number, force_zone_letter=letter)
                y[negmask] = -y[negmask]
        return x, y

    def _lonlat(self, lon: np.ndarray=None, lat: np.ndarray=None, strict=False) -> tuple[np.ndarray, np.ndarray]:
        """Converts list of points to longitude and latitude if necessary.

        Input is assumed to be in the native format (UTM for cartesian, WGS84
        for spherical).

        If strict=True, then (None, None) is returned if grid is cartesian.
        """
        if not self._structure_initialized():
            return None

        if self.is_cartesian() and (strict or self.strict):
            return None, None

        if self.is_cartesian():
            number, letter = self.utm()
            lat, lon = utm.to_latlon(lon, np.mod(lat,10_000_000), zone_number=number, zone_letter=letter, strict = False)
        return lon, lat

    def edges(self, coord: str, native: bool=False, strict=False) -> tuple[float, float]:
        """Min and max values of x. Conversion made for sperical grids."""
        if not self._structure_initialized():
            return (None, None)

        if self._xy_dict.get(coord) is None:
            print("coord need to be 'x', 'y', 'lon' or 'lat'.")
            return

        if native:
            method = f'native_{self._xy_dict[coord]}' # e.g. self.native_x()
            val = getattr(self, method)()
        else:
            # Need to get full list (not gridded) to not rotate the coordinates
            vals = {}
            vals['x'], vals['y'] = getattr(self, self._x_to_xy_dict.get(coord))(strict=(strict or self.strict)) # e.g. coord = x -> self.xy()
            val = vals[self._xy_dict[coord]]

        if val is None:
            return (None, None)

        return np.min(val), np.max(val)

    def nx(self) -> int:
        """Length of x/lon-vector."""
        if not self._structure_initialized():
            return 0
        return len(self.native_x())

    def ny(self):
        """Length of y/lat-vector."""
        if not self._structure_initialized():
            return 0
        return len(self.native_y())

    def dx(self, strict=False):
        """Mean grid spacing of the x vector. Conversion made for
        spherical grids."""
        if not self._structure_initialized():
            return None

        if not self.is_cartesian() and (strict or self.strict):
            return None

        if self.nx() == 1:
            return 0.

        edges = self.edges('x')
        return (edges[1]-edges[0])/(self.nx()-1)

    def dy(self, strict=False):
        """Mean grid spacing of the y vector. Conversion made for
        spherical grids."""
        if not self._structure_initialized():
            return None

        if not self.is_cartesian() and (strict or self.strict):
            return None

        if self.ny() == 1:
            return 0.

        edges = self.edges('y')
        return (edges[1]-edges[0])/(self.ny()-1)

    def dlon(self, strict=False):
        """Mean grid spacing of the longitude vector. Conversion made for
        cartesian grids."""
        if not self._structure_initialized():
            return None

        if self.is_cartesian() and (strict or self.strict):
            return None
        if self.nx() == 1:
            return 0.

        edges = self.edges('lon')
        return (edges[1]-edges[0])/(self.nx()-1)

    def dlat(self, strict=False):
        """Mean grid spacing of the latitude vector. Conversion made for
        cartesian grids."""
        if not self._structure_initialized():
            return None

        if self.is_cartesian() and (strict or self.strict):
            return None
        if self.ny() == 1:
            return 0.

        edges = self.edges('lat')
        return (edges[1]-edges[0])/(self.ny()-1)

    def native_dx(self):
        """Mean grid spacing of x vector for cartesian grids.
        Mean grid spacing of lon vector for spherical grids."""
        if not self._structure_initialized():
            return None

        if self.nx() == 1:
            return 0.

        edges = self.edges('x', native=True)
        return (edges[1]-edges[0])/(self.nx()-1)

    def native_dy(self):
        """Mean grid spacing of y vector for cartesian grids.
        Mean grid spacing of lat vector for spherical grids."""
        if not self._structure_initialized():
            return None

        if self.ny() == 1:
            return 0.

        edges = self.edges('y', native=True)
        return (edges[1]-edges[0])/(self.ny()-1)

    def metadata(self) -> dict:
        """Return metadata of the dataset:

        """
        if not self._structure_initialized():
            return None
        return self.ds().attrs

    def set_metadata(self, metadata: dict, append=False) -> None:
        if not self._structure_initialized():
            return None
        if append:
            old_metadata = self.metadata()
            old_metadata.update(metadata)
            metadata = old_metadata
        self.ds_manager.set_attrs(metadata)

    @property
    def x_str(self) -> str:
        """Return string compatible with the type of spacing used:

        'x' for cartesian grid.
        'lon' for spherical grid.
        """
        if not self._structure_initialized():
            return None
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
        if not self._structure_initialized():
            return None
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
        if not self._structure_initialized():
            return None
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

    if (x is not None) and (y is not None):
        xy = True
        native_x = 'x'
        native_y = 'y'
        xvec = x
        yvec = y

    #if (lon is not None and lon != none_tuple) and (lat is not None and lat != none_tuple):
    if (lon is not None) and (lat is not None):
        lonlat = True
        native_x = 'lon'
        native_y = 'lat'
        xvec = lon
        yvec = lat

    if xy and lonlat:
        raise ValueError("Can't set both lon/lat and x/y!")

    if not xy and not lonlat:
        raise ValueError('Have to set either lon/lat or x/y!')

    return native_x, native_y, xvec, yvec
