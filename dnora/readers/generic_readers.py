from __future__ import annotations
from .abstract_readers import DataReader, PointDataReader
from data_sources import DataSource
import pandas as pd
import numpy as np
import xarray as xr
from dnora import aux_funcs
from pathlib import Path

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from dnora_types import DnoraDataType


class ConstantGrid(DataReader):
    def __init__(self, **kwargs):
        """E.g. ConstantGrid(u=1, v=2)"""
        self.values = kwargs

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        **kwargs,
    ):
        time = pd.date_range(start=start_time, end=end_time, freq="H").values

        coord_dict = {}
        coord_dict["lon"] = grid.lon(strict=True)
        coord_dict["lat"] = grid.lat(strict=True)
        coord_dict["x"] = grid.x(strict=True)
        coord_dict["y"] = grid.y(strict=True)
        if "time" in obj_type.value._coord_manager.added_coords():
            coord_dict["time"] = time
            obj_size = (len(time), grid.ny(), grid.nx())
        else:
            obj_size = grid.size()

        variables = obj_type.value._coord_manager.added_vars().keys()

        data_dict = {}
        for key in variables:
            val = self.values.get(key)
            if val is None:
                val = kwargs.get(key, 1)
            data_dict[key] = np.full(obj_size, val)

        meta_dict = {}

        # Create metaparameters based on standard short names
        metaparameter_dict = self.create_metaparameter_dict(data_dict.keys())

        return coord_dict, data_dict, meta_dict, metaparameter_dict


class Netcdf(DataReader):
    def __init__(self, files: str) -> None:
        self.files = files

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        **kwargs,
    ):
        ds = xr.open_mfdataset(Path(folder).joinpath(self.files))
        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds)
        coord_dict = {}
        # obj_type.value._coord_manager.added_coords()

        times = slice(start_time, end_time)
        if x is None:
            lons = slice(grid.edges("lon")[0], grid.edges("lon")[-1])
            lats = slice(grid.edges("lat")[0], grid.edges("lat")[-1])
            try:
                ds = ds.sel(lon=lons, lat=lats, time=times)
            except:
                ds = ds.sel(longitude=lons, latitude=lats, time=times)
        else:
            xs = slice(grid.edges("x")[0], grid.edges("x")[-1])
            ys = slice(grid.edges("y")[0], grid.edges("y")[-1])
            ds = ds.sel(x=xs, y=ys, time=times)

        data_dict = {}
        metaparameter_dict = {}
        for var, meta_var in obj_type.value.meta_dict.items():
            ds_var = meta_var.find_me_in_ds(ds)
            ds_data = ds.get(ds_var)
            if ds_data is not None:
                data_dict[var] = ds_data.values
                metaparameter_dict[var] = meta_var

        meta_dict = ds.attrs
        return coord_dict, data_dict, meta_dict, metaparameter_dict


class ConstantPoint(PointDataReader):
    def __init__(self, **kwargs):
        """E.g. ConstantGrid(u=1, v=2)"""
        self.values = kwargs

    def get_coordinates(self, grid, start_time, source, folder):
        lon_all, lat_all = grid.lonlat(strict=True)
        x_all, y_all = grid.xy(strict=True)

        return lon_all, lat_all, x_all, y_all

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        inds,
        **kwargs,
    ):
        time = pd.date_range(start=start_time, end=end_time, freq="H").values

        if obj_type in [DnoraDataType.SPECTRA, DnoraDataType.SPECTRA1D]:
            freq = np.linspace(0.1, 1, 10)
        if obj_type == DnoraDataType.SPECTRA:
            dirs = np.linspace(0, 350, 36).astype(int)

        coord_dict = {}

        obj_size = []
        if "time" in obj_type.value._coord_manager.added_coords():
            coord_dict["time"] = time
            obj_size.append(len(time))
        obj_size.append(len(inds))
        if "freq" in obj_type.value._coord_manager.added_coords():
            coord_dict["freq"] = freq
            obj_size.append(len(freq))
        if "dirs" in obj_type.value._coord_manager.added_coords():
            coord_dict["dirs"] = dirs
            obj_size.append(len(dirs))

        obj_size = tuple(obj_size)

        lon_all, lat_all, x_all, y_all = self.get_coordinates(
            grid, start_time, source, ""
        )
        if lon_all is not None:
            coord_dict["lon"] = lon_all[inds]
        if lat_all is not None:
            coord_dict["lat"] = lat_all[inds]
        if x_all is not None:
            coord_dict["x"] = x_all[inds]
        if y_all is not None:
            coord_dict["y"] = y_all[inds]

        variables = obj_type.value._coord_manager.added_vars().keys()

        data_dict = {}
        if variables:  # If the object has addeed variables, create those
            for key in variables:
                val = self.values.get(key)
                if val is None:
                    val = kwargs.get(key, 1)
                data_dict[key] = np.full(obj_size, val)
        else:  # If not, create the ones provided by the user and assume they will be dynamically added
            for key, val in self.values.items():
                data_dict[key] = np.full(obj_size, val)
            for key, val in kwargs.items():
                data_dict[key] = np.full(obj_size, val)

        meta_dict = {}
        metaparameter_dict = self.create_metaparameter_dict(data_dict.keys())
        return coord_dict, data_dict, meta_dict, metaparameter_dict
