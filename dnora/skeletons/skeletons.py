import numpy as np
import pandas as pd
import xarray as xr
import utm

class Skeleton:
    _native_x_var = 'x'
    _native_y_var = 'y'

    def x(self) -> np.ndarray:
        if self._native_x_var == 'lon':
            x, __, __, __ = utm.from_latlon(self.lat(), self.lon(), force_zone_number=33, force_zone_letter='W')
            return x

        if hasattr(self, 'data') and hasattr(self.data, 'x'):
            return self.data.x.values
        else:
            return None

    def y(self) -> np.ndarray:
        if self._native_y_var == 'lat':
            __, y, __, __ = utm.from_latlon(self.lat(), self.lon(), force_zone_number=33, force_zone_letter='W')
            return y

        if hasattr(self, 'data') and hasattr(self.data, 'y') :
            return self.data.y.values
        else:
            return None

    def lon(self) -> np.ndarray:
        if self._native_x_var == 'x':
            __, lon = utm.to_latlon(self.x(), self.y(), 33, zone_letter = 'W', strict = False)
            return lon

        if hasattr(self, 'data') and hasattr(self.data, 'lon'):
            return self.data.lon.values
        else:
            return None

    def lat(self) -> np.ndarray:
        if self._native_y_var == 'y':
            lat, __ = utm.to_latlon(self.x(), self.y(), 33, zone_letter = 'W', strict = False)
            return lat

        if hasattr(self, 'data') and hasattr(self.data, 'lat') :
            return self.data.lat.values
        else:
            return None

    def native_x(self) -> np.ndarray:
        if self._native_x_var == 'x':
            return self.x()
        elif self._native_x_var == 'lon':
            return self.lon()
        else:
            return None

    def native_y(self) -> np.ndarray:
        if self._native_y_var == 'y':
            return self.y()
        elif self._native_y_var == 'lat':
            return self.lat()
        else:
            return None

    def time(self):
        if hasattr(self, 'data') and hasattr(self.data, 'time'):
            return self.data.time.values
        else:
            return None


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
        raise Exception("Can't set both lon/lat and x/y!")

    if not xy and not lonlat:
        raise Exception('Have to set either lon/lat or x/y!')

    return native_x, native_y, xvec, yvec

class PointSkeleton(Skeleton):
    def __init__(self, x=None, y=None, lon=None, lat=None, time=None):
        native_x, native_y, xvec, yvec = check_input_consistency(x, y, lon, lat)
        if len(xvec) != len(yvec):
            raise Exception('x and y vector has to be equally long!')
        self._native_x_var = native_x
        self._native_y_var = native_y

        if time is not None and len(time) != len(xvec):
            raise Exception('Length of time vector needs to agree with spatial x/y vectors!')

        self.data = self._create_xr(station=np.arange(len(xvec)),
                                    x=xvec,
                                    y=yvec,
                                    time=time)


    def _create_xr(self, station: np.ndarray, x: np.ndarray, y: np.ndarray, time=None) -> xr.Dataset:
        coords_dict = {'station': station}
        vars_dict = {self._native_x_var: (['station'], x), self._native_y_var: (['station'], y)}
        if time is not None:
            vars_dict['time'] = (['station'], time)
        return xr.Dataset(coords=coords_dict, data_vars=vars_dict)

    def station(self) -> np.ndarray:
        if hasattr(self.data, 'station') :
            return self.data.station.values
        else:
            return None

class GriddedSkeleton(Skeleton):
    def _create_xr(self, x: np.ndarray, y: np.ndarray, time=None) -> xr.Dataset:
        coords_dict = {self._native_x_var: x, self._native_y_var: y}
        if time is not None:
            coords_dict['time'] = time
        return xr.Dataset(coords=coords_dict)
