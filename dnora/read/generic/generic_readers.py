from dnora.cacher.caching_strategies import CachingStrategy
from dnora.read.abstract_readers import PointDataReader, DataReader, SpectralDataReader
import geo_parameters as gp
import pandas as pd
import numpy as np
import xarray as xr
from dnora import utils
from pathlib import Path
from dnora.type_manager.dnora_types import DnoraDataType
from dnora.type_manager.data_sources import DataSource
from dnora.type_manager.spectral_conventions import convert_2d_to_1d, SpectralConvention

from dnora.type_manager.dnora_objects import dnora_objects
from dnora.aux_funcs import get_url
from dnora import msg, utils
from .constant_funcs import create_constant_array, print_constant_values
from dnora.read.data_var_decoding import read_data_vars, compile_data_vars
from copy import copy
from geo_skeletons import GriddedSkeleton


class PointNetcdf(SpectralDataReader):
    def default_data_source(self) -> DataSource:
        return DataSource.LOCAL

    def _caching_strategy(self) -> CachingStrategy:
        return CachingStrategy.DontCacheMe

    def __init__(self, files: list[str] = None):
        self.files = files

    def get_coordinates(
        self, grid, start_time, source, folder, filename: list[str] = None, **kwargs
    ):
        filename = filename or self.files
        if filename is None:
            raise ValueError("Provide at least one filename!")
        filepath = get_url(folder, filename, get_list=True)
        ds = xr.open_dataset(filepath[0])
        lon, lat, x, y = utils.grid.get_coordinates_from_ds(ds)
        self.set_convention(ds.attrs.get("dnora_spectral_convention", "ocean"))
        return {"lon": lon, "lat": lat, "x": x, "y": y}

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        inds: list[int],
        filename: list[str] = None,
        **kwargs,
    ):

        filename = filename or self.files

        if filename is None:
            raise ValueError("Provide at least one filename!")
        filepath = get_url(folder, filename, get_list=True)
        ds = xr.open_mfdataset(filepath)
        for fn in filepath:
            msg.from_file(fn)
        lon, lat, x, y = utils.grid.get_coordinates_from_ds(ds)

        times = slice(start_time, end_time)

        ds = ds.sel(inds=inds, time=times)

        lon, lat, x, y = utils.grid.get_coordinates_from_ds(ds)
        coord_dict = {"x": x, "y": y, "lon": lon, "lat": lat}
        for c in list(ds.coords):
            if c not in ["inds"]:
                coord_dict[c] = ds.get(c).values

        data_dict = {}
        for var in dnora_objects.get(obj_type).core.data_vars():
            meta_var = dnora_objects.get(obj_type).core.meta_parameter(var)
            ds_var = meta_var.find_me_in_ds(ds)
            ds_data = ds.get(ds_var)
            if ds_data is not None:
                data_dict[meta_var] = ds_data.values

        meta_dict = ds.attrs
        return coord_dict, data_dict, meta_dict


class Netcdf(DataReader):
    def default_data_source(self) -> DataSource:
        return DataSource.LOCAL

    def _caching_strategy(self) -> CachingStrategy:
        return CachingStrategy.DontCacheMe

    def __init__(self, files: list[str] = None):
        self.files = files

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        filename: list[str] = None,
        expansion_factor=1.0,
        **kwargs,
    ):

        filename = filename or self.files

        if filename is None:
            raise ValueError("Provide at least one filename!")
        filepath = get_url(folder, filename, get_list=True)
        ds = xr.open_mfdataset(filepath)

        msg.from_multifile(filepath)
        lon, lat, x, y = utils.grid.get_coordinates_from_ds(ds)

        times = slice(start_time, end_time)
        if x is None:
            lons, lats = utils.grid.expand_area(
                grid.edges("lon"), grid.edges("lat"), expansion_factor
            )
            try:
                ds = ds.sel(lon=slice(*lons), lat=slice(*lats), time=times)
            except:
                ds = ds.sel(longitude=slice(*lons), latitude=slice(*lats), time=times)
        else:
            xs, ys = utils.grid.expand_area(
                grid.edges("x"), grid.edges("y"), expansion_factor
            )
            ds = ds.sel(x=slice(*xs), y=slice(*ys), time=times)

        lon, lat, x, y = utils.grid.get_coordinates_from_ds(ds)
        coord_dict = {"x": x, "y": y, "lon": lon, "lat": lat}

        for c in list(ds.coords):
            if c not in ["lon", "lat", "longitudes", "latitudes"]:
                coord_dict[c] = ds.get(c).values
        data_dict = {}
        metaparameter_dict = {}
        for var, meta_var in dnora_objects.get(obj_type).meta_dict.items():
            ds_var = meta_var.find_me_in_ds(ds)
            ds_data = ds.get(ds_var)
            if ds_data is not None:
                data_dict[meta_var] = ds_data.values
                # metaparameter_dict[var] = meta_var

        meta_dict = ds.attrs

        return coord_dict, data_dict, meta_dict


class ConstantData(SpectralDataReader):
    _class_default_var_values = {"default": 1.0}
    _class_default_coord_values = {
        "freq": np.linspace(0.1, 1, 10),
        "dirs": np.linspace(0, 350, 36),
        "default": np.linspace(0.0, 10.0, 11),
    }
    _class_default_peak_values = {"freq": 0.3, "dirs": 0, "default": 0}

    def default_data_source(self) -> DataSource:
        return DataSource.CREATION

    def _caching_strategy(self) -> CachingStrategy:
        return CachingStrategy.DontCacheMe

    def get_coordinates(self, grid, start_time, source, folder, **kwargs):
        lon, lat = grid.lonlat(strict=True)
        x, y = grid.xy(strict=True)
        return {"x": x, "y": y, "lon": lon, "lat": lat}

    def __init__(
        self,
        vars: dict = None,
        coords: dict = None,
        peaks: dict = None,
        convention: SpectralConvention | str = SpectralConvention.OCEAN,
    ):
        """E.g. ConstantData(vars={'u':1, 'v':2})"""
        vars = vars or {}
        self._default_var_values = copy(self._class_default_var_values)
        self._default_var_values.update(vars)

        coords = coords or {}
        self._default_coord_values = copy(self._class_default_coord_values)
        self._default_coord_values.update(coords)

        peaks = peaks or {}
        self._default_peak_values = copy(self._class_default_peak_values)
        self._default_peak_values.update(peaks)

        self.set_convention(convention)

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        new_vars: dict = None,
        time: list[str] = None,
        dt: int = None,
        force_type: str = "",
        **kwargs,
    ):
        """Variables that are a part of the class can be given as a keyword, e.g. u=1.
        To add variables dynamically, give the values as a dict new_vars, e.g. new_vars={'hs': 6}

        time: determines start time if constant values given as an array
        dt [h]: Can be given instad of time to determine start times of array values

        If ModelRun is defined for 2020-01-01 00:00 - 2020-01-01 23:00, set 12 + 12 hour wind speed with:
        u=[1, 2, v=[1,2], time=['2020-01-01 00:00', '2020-01-01 12:00']

        OR

        u=[1, 2, v=[1,2], dt=12

        The created object will be of same type (spherical/cartesian) as the grid. To force another type, use:
        force_type = 'spherical'/'cartesian'
        """

        def non_fp_inds():
            if self._default_peak_values.get("freq") is not None and obj_type in [
                DnoraDataType.SPECTRA,
                DnoraDataType.SPECTRA1D,
            ]:
                """Set everything except dominant frequency to 0"""
                non_fp_ind = np.argmin(
                    np.abs(
                        coord_dict.get("freq") - self._default_peak_values.get("freq")
                    )
                ) != np.arange(len(coord_dict.get("freq")))
                return non_fp_ind
            return None

        def non_dirp_inds():
            if (
                self._default_peak_values.get("dirs") is not None
                and obj_type == DnoraDataType.SPECTRA
            ):
                """Set everything except dominant direction to 0"""
                non_dirp_ind = np.argmin(
                    np.abs(
                        coord_dict.get("dirs") - self._default_peak_values.get("dirs")
                    )
                ) != np.arange(len(coord_dict.get("dirs")))
                return non_dirp_ind
            return None

        coord_dict = {}

        # Determine size of the object to create
        object_coords = dnora_objects.get(obj_type).core.coords("nonspatial")
        obj_size = []

        # Time is first if exists
        if "time" in object_coords:
            time_vec = pd.date_range(start=start_time, end=end_time, freq="h")
            coord_dict["time"] = time_vec
            obj_size.append(len(time_vec))
            object_coords.remove("time")
        else:
            time_vec = None

        # Then spatial coords
        coord_dict, obj_size, zone_number, zone_letter = calculate_spatial_coords(
            coord_dict, obj_size, grid, obj_type, force_type
        )

        # Then all other coords
        for obj_coord in object_coords:
            default = self._default_coord_values.get("default")
            val = kwargs.get(
                "obj_coord", self._default_coord_values.get(obj_coord, default)
            )
            coord_dict[obj_coord] = val
            obj_size.append(len(val))

        obj_size = tuple(obj_size)

        # Get values for data variables
        dnora_obj = dnora_objects.get(obj_type)

        data_val_dict = create_datvar_val_dict(
            kwargs,
            default_values=self._default_var_values,
            object_vars=dnora_obj.core.data_vars(),
        )
        data_dict = {}
        non_fp_ind = non_fp_inds()
        non_dirp_ind = non_dirp_inds()
        for key, val in data_val_dict.items():
            data_dict[key] = create_constant_array(val, time, obj_size, time_vec)
            if non_fp_ind is not None:
                data_dict[key][:, :, non_fp_ind, ...] = 0
            elif non_dirp_ind is not None:
                data_dict[key][:, :, :, non_dirp_ind] = 0

        print_constant_values(data_dict, obj_type, time_vec)

        if zone_number is not None and zone_letter is not None:
            meta_dict = {"utm_zone": zone_number, "zone_letter": zone_letter}
        else:
            meta_dict = {}

        return coord_dict, data_dict, meta_dict


def create_datvar_val_dict(kwargs, default_values, object_vars):
    """Creates a dictionary of data variables and values"""
    new_vars_dict = kwargs.get("new_vars", {})
    new_vars = list(
        new_vars_dict.keys()
    )  # Given explicitly since not yet added to class

    datavar_list = list(object_vars) + new_vars

    data_dict = {}
    for key in datavar_list:
        default = default_values.get(key, default_values["default"])
        data_dict[key] = new_vars_dict.get(key, kwargs.get(key, default))

    return data_dict


def calculate_spatial_coords(coord_dict, obj_size, grid, obj_type, force_type):
    """Determines coords and their length while accounting for that either the object or the grid might be gridded or not."""
    obj_gridded = issubclass(dnora_objects.get(obj_type), GriddedSkeleton)
    if not obj_gridded:
        """If object is not gridded, just define it at each grid point"""
        if force_type == "spherical":
            x, y = grid.lonlat()
            x_str, y_str = "lon", "lat"
        elif force_type == "cartesian":
            x, y = grid.xy()
            x_str, y_str = "x", "y"
        else:
            x, y = grid.xy(native=True)
            x_str, y_str = grid.core.x_str, grid.core.y_str

        obj_size.append(len(y))

    elif grid.is_gridded():
        """If both not gridded, then use identical points to grid"""
        if force_type == "spherical":
            x, y = grid.lon(), grid.lat()
            x_str, y_str = "lon", "lat"
        elif force_type == "cartesian":
            x, y = grid.x, grid.y()
            x_str, y_str = "x", "y"
        else:
            x, y = grid.x(native=True), grid.y(native=True)
            x_str, y_str = grid.core.x_str, grid.core.y_str

        obj_size.append(len(y))
        obj_size.append(len(x))

    else:
        """If object is gridded, but grid is not, keep approximate number of points"""
        ny = np.ceil(np.sqrt(grid.ny())).astype(int)
        nx = np.ceil(np.sqrt(grid.nx())).astype(int)

        if force_type == "spherical":
            x_edge = grid.edges("lon")
            y_edge = grid.edges("lat")
            x = np.linspace(x_edge[0], x_edge[1], nx)
            y = np.linspace(y_edge[0], y_edge[1], nx)
            x_str, y_str = "lon", "lat"
        elif force_type == "cartesian":
            x_edge = grid.edges("x")
            y_edge = grid.edges("y")
            x = np.linspace(x_edge[0], x_edge[1], nx)
            y = np.linspace(y_edge[0], y_edge[1], nx)
            x_str, y_str = "x", "y"
        else:
            x_edge = grid.edges("x", native=True)
            y_edge = grid.edges("y", native=True)
            x = np.linspace(x_edge[0], x_edge[1], nx)
            y = np.linspace(y_edge[0], y_edge[1], nx)
            x_str, y_str = grid.core.x_str, grid.core.y_str

        obj_size.append(ny)
        obj_size.append(nx)

    coord_dict[x_str] = x
    coord_dict[y_str] = y

    if x_str == "x":
        zone_number, zone_letter = grid.utm.zone()
    else:
        zone_number, zone_letter = None, None

    return coord_dict, obj_size, zone_number, zone_letter
