import numpy as np
import pandas as pd
import xarray as xr
import utm
from copy import copy

class Skeleton:
    _x_str = 'x'
    _y_str = 'y'
    _name = 'LonelySkeleton'
    _zone_number = 33
    _zone_letter = 'W'
    strict = False # If this is True, no coordinate conversions will be done (return None instead)

    def set_utm(self, zone_number: int=33, zone_letter: str='W'):
        self._zone_number = copy(zone_number)
        self._zone_letter = copy(zone_letter)

    def x(self, strict=False) -> np.ndarray:
        if self.x_str == 'lon' and (strict or self.strict):
            return None

        if self.x_str == 'lon':
            lat = np.median(self.lat())
            x, __, __, __ = utm.from_latlon(lat, self.lon(), force_zone_number=self._zone_number, force_zone_letter=self._zone_letter)
            return x

        if hasattr(self, 'data') and hasattr(self.data, 'x'):
            return self.data.x.values
        else:
            return None

    def y(self, strict=False) -> np.ndarray:
        if self.y_str == 'lat' and (strict or self.strict):
            return None

        if self.y_str == 'lat':
            lon = np.median(self.lon())
            __, y, __, __ = utm.from_latlon(self.lat(), lon, force_zone_number=self._zone_number, force_zone_letter=self._zone_letter)
            return y

        if hasattr(self, 'data') and hasattr(self.data, 'y') :
            return self.data.y.values
        else:
            return None

    def lon(self, strict=False) -> np.ndarray:
        if self.x_str == 'x' and (strict or self.strict):
            return None

        if self.x_str == 'x':
            y = np.median(self.y())
            __, lon = utm.to_latlon(self.x(), y, self._zone_number, zone_letter=self._zone_letter, strict = False)
            return lon

        if hasattr(self, 'data') and hasattr(self.data, 'lon'):
            return self.data.lon.values
        else:
            return None

    def lat(self, strict=False) -> np.ndarray:
        if self.y_str == 'y' and (strict or self.strict):
            return None

        if self.y_str == 'y':
            x = np.median(self.x())
            lat, __ = utm.to_latlon(x, self.y(), self._zone_number, zone_letter=self._zone_letter, strict = False)
            return lat

        if hasattr(self, 'data') and hasattr(self.data, 'lat') :
            return self.data.lat.values
        else:
            return None

    def native_x(self) -> np.ndarray:
        if self.x_str == 'x':
            return self.x()
        elif self.x_str == 'lon':
            return self.lon()
        else:
            return None

    def native_y(self) -> np.ndarray:
        if self.y_str == 'y':
            return self.y()
        elif self.y_str == 'lat':
            return self.lat()
        else:
            return None

    def _lonlat(self, lon: np.ndarray=None, lat: np.ndarray=None, strict=False):
        """Converts list of points to longitude and latitude if necessary.

        Input is assumed to be in the native format.
        """
        if self.x_str == 'x' and (strict or self.strict):
            return None, None

        if self.x_str == 'x':
            lat, lon = utm.to_latlon(lon, lat, self._zone_number, zone_letter=self._zone_letter, strict = False)
        return lon, lat

    def _xy(self, x: np.ndarray=None, y: np.ndarray=None, strict=False):
        """Converts list of points to x and y (UTM) if necessary.

        Input is assumed to be in the native format.
        """
        if self.x_str == 'lon' and (strict or self.strict):
            return None, None

        if self.x_str == 'lon':
            x, y, __, __ = utm.from_latlon(y, x, force_zone_number=self._zone_number, force_zone_letter=self._zone_letter)
        return x, y

    def _additional_coords(self):
        if hasattr(self, 'data'):
            return list(set(self.data.coords)-{'x','y','inds','lon','lat','time'})
        return None

    def  _additional_coord_val(self, coord):
        return eval(f"self.{coord}()")

    def lon_edges(self) -> tuple:
        return np.min(self.lon()), np.max(self.lon())

    def lat_edges(self) -> tuple:
        return np.min(self.lat()), np.max(self.lat())

    def x_edges(self) -> tuple:
        return np.min(self.x()), np.max(self.x())

    def y_edges(self) -> tuple:
        return np.min(self.y()), np.max(self.y())

    def native_x_edges(self) -> tuple:
        return np.min(self.native_x()), np.max(self.native_x())

    def native_y_edges(self) -> tuple:
        return np.min(self.native_y()), np.max(self.native_y())

    def nx(self):
        return len(self.native_x())

    def ny(self):
        return len(self.native_y())

    def nt(self):
        if self.time() is not None:
            return len(self.time())
        return None

    def dlon(self):
        if self.nx() == 1:
            return 0.
        return (self.lon_edges()[1]-self.lon_edges()[0])/(self.nx()-1)

    def dlat(self):
        if self.ny() == 1:
            return 0.
        return (self.lat_edges()[1]-self.lat_edges()[0])/(self.ny()-1)

    def dx(self):
        if self.nx() == 1:
            return 0.
        return (self.x_edges()[1]-self.x_edges()[0])/(self.nx()-1)

    def dy(self):
        if self.ny() == 1:
            return 0.
        return (self.y_edges()[1]-self.y_edges()[0])/(self.ny()-1)

    def native_dx(self):
        if self.nx() == 1:
            return 0.
        return (self.native_x_edges()[1]-self.native_x_edges()[0])/(self.nx()-1)

    def native_dy(self):
        if self.ny() == 1:
            return 0.
        return (self.native_y_edges()[1]-self.native_y_edges()[0])/(self.nx()-1)

    @property
    def x_str(self) -> str:
        return self._x_str

    @x_str.setter
    def x_str(self, new_str):
        if new_str in ['x', 'lon']:
            self._x_str = new_str
        else:
            raise ValueError("x_str need to be 'x' or 'lon'")

    @property
    def y_str(self) -> str:
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

    @x_str.setter
    def name(self, new_name):
        if isinstance(new_name, str):
            self._name = new_name
        else:
            raise ValueError("name needs to be a string")

    def time(self):
        if hasattr(self, 'data') and hasattr(self.data, 'time'):
            return self.data.time.values
        else:
            return None

    def _ds_vars_dict(self):
        """Return variable dictionary for creating xarray Dataset.
        """
        raise NotImplementedError("This method depends on the specific gridding and needs to be defined in the subclass (e.g. GriddedSkeleton or PointSkeleton)")

    def _ds_coords_dict(self):
        """Return variable dictionary for creating xarray Dataset.
        """
        raise NotImplementedError("This method depends on the specific gridding and needs to be defined in the subclass (e.g. GriddedSkeleton or PointSkeleton)")

    def _init_ds(self) -> xr.Dataset:
        """Return a Dataset with only the grid coordinates and time.
        """
        raise NotImplementedError("This method depends on the specific gridding and needs to be defined in the subclass (e.g. GriddedSkeleton or PointSkeleton)")

    def compile_to_ds(self, data: np.ndarray, data_name:str):#, additional_coords: dict=None):
        def check_consistency():
            for i, key in enumerate(coords_dict.keys()):
                if i > len(data.shape)-1:
                    raise Exception(f'{key} coordinate is {len(coords_dict[key])} long, but that dimension doesnt exist in the data!!!')
                if len(coords_dict[key]) != data.shape[i]:
                    raise Exception(f'{key} coordinate is {len(coords_dict[key])} long, but size of data in that dimension (dim {i}) is {data.shape[i]}!!!')

        # Coordinates
        coords_dict = self._ds_coords_dict()
        #if additional_coords is not None:
        for coord in self._additional_coords():
            coords_dict[coord] = self._additional_coord_val(coord)


        # Data variables
        vars_dict = self._ds_vars_dict()
        vars_dict[data_name] = (list(coords_dict.keys()),data)

        check_consistency()

        return xr.Dataset(data_vars=vars_dict, coords=coords_dict)

    def ds(self):
        if not hasattr(self, 'data'):
            return None
        return self.data

    def _set_data(self, data: np.ndarray, data_name: str) -> None:
        self.merge_in_ds(self.compile_to_ds(data, data_name))

    def get(self, data_name: str):
        """Gets data from Dataset"""
        ds = self.ds()
        if ds is None:
            return None

        data = ds.get(data_name)

        if data is None:
            return None

        return data.values.copy()

    def _create_structure(self, x=None, y=None, lon=None, lat=None, time=None, **kwargs):
        native_x, native_y, xvec, yvec = check_input_consistency(x, y, lon, lat)

        self.x_str = native_x
        self.y_str = native_y

        return self._init_ds(x=xvec,y=yvec, time=time, **kwargs)

    def merge_in_ds(self, ds_list: list[xr.Dataset]):
        if not isinstance(ds_list, list):
            ds_list = [ds_list]
        for ds in ds_list:
            self.data = ds.merge(self.data, compat='override')


def check_input_consistency(x, y, lon, lat):
    xy = False
    lonlat = False

    if x is not None and y is not None:
        xy = True
        native_x = 'x'
        native_y = 'y'
        xvec = x
        yvec = y

    if lon is not None and lat is not None:
        lonlat = True
        native_x = 'lon'
        native_y = 'lat'
        xvec = lon
        yvec = lat

    if xy and lonlat:
        raise ValueError("Can't set both lon/lat and x/y!")

    if not xy and not lonlat:
        raise ValueError('Have to set either lon/lat or x/y!')

    return native_x, native_y, np.array(xvec), np.array(yvec)
