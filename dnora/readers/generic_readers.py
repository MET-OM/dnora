from .abstract_readers import DataReader, SpectralDataReader
from dnora.data_sources import DataSource
import pandas as pd
import numpy as np
import xarray as xr
from dnora import aux_funcs
from pathlib import Path
from dnora.metaparameter.parameter_funcs import create_metaparameter_dict
from dnora.dnora_types import DnoraDataType

from dnora.spectral_conventions import convert_2d_to_1d, SpectralConvention


class ConstantGriddedData(DataReader):
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


class ConstantPointData(SpectralDataReader):
    def __init__(
        self,
        fp: float | None = 0.3,
        dirp: int | None = 0,
        convention: SpectralConvention | str = SpectralConvention.OCEAN,
        **kwargs,
    ):
        """**kwargs are possible extra coordinates
        E.g. ConstantGrid(z=np.linspace(0,10,11))

        fp/dirp used when importing spectra only
        """
        self.extra_coords = kwargs
        if self.extra_coords.get("freq") is None:
            self.extra_coords["freq"] = np.linspace(0.1, 1, 10)
        if self.extra_coords.get("dirs") is None:
            self.extra_coords["dirs"] = np.linspace(0, 350, 36)

        self.dirp = dirp
        self.fp = fp
        self.provided_convention = convention

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
        coord_dict = {}

        obj_size = []
        # Time is always first coord if exists
        if "time" in obj_type.value._coord_manager.added_coords():
            coord_val = self.extra_coords.get(
                "time", pd.date_range(start=start_time, end=end_time, freq="H").values
            )
            coord_dict["time"] = coord_val
            obj_size.append(len(coord_val))

        # Inds always second (or first if time doesn't exist)
        obj_size.append(len(inds))

        for added_coord in obj_type.value._coord_manager.added_coords():
            if added_coord != "time":
                coord_val = self.extra_coords.get(added_coord, np.linspace(1, 100, 100))
                coord_dict[added_coord] = coord_val
                obj_size.append(len(coord_val))

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

        if self.fp is not None and obj_type in [
            DnoraDataType.SPECTRA,
            DnoraDataType.SPECTRA1D,
        ]:
            """Set everything except dominant frequency to 0"""
            non_fp_ind = np.argmin(
                np.abs(self.extra_coords.get("freq") - self.fp)
            ) != np.arange(len(self.extra_coords.get("freq")))
        else:
            non_fp_ind = np.ones(len(self.extra_coords.get("freq"))).astype(bool)

        if self.dirp is not None and obj_type == DnoraDataType.SPECTRA:
            """Set everything except dominant direction to 0"""
            non_dirp_ind = np.argmin(
                np.abs(self.extra_coords.get("dirs") - self.dirp)
            ) != np.arange(len(self.extra_coords.get("dirs")))
        else:
            non_dirp_ind = np.ones(len(self.extra_coords.get("dirs"))).astype(bool)

        variables = obj_type.value._coord_manager.added_vars().keys()

        data_dict = {}
        if variables:  # If the object has addeed variables, create those
            for key in variables:
                val = kwargs.get(key, 1)
                data_dict[key] = np.full(obj_size, val)

        else:  # If not, create the ones provided by the user and assume they will be dynamically added (if not, they are just dumped)
            for key, val in kwargs.items():
                data_dict[key] = np.full(obj_size, val)

        for key in data_dict:
            if obj_type == DnoraDataType.SPECTRA:
                data_dict[key][:, :, non_fp_ind, :] = 0
                data_dict[key][:, :, :, non_dirp_ind] = 0
            elif obj_type == DnoraDataType.SPECTRA1D:
                data_dict[key][:, :, non_fp_ind] = 0

        if obj_type == DnoraDataType.SPECTRA:
            self.set_convention(self.provided_convention)
        if obj_type == DnoraDataType.SPECTRA1D:
            self.set_convention(convert_2d_to_1d(self.provided_convention))

        meta_dict = {}
        metaparameter_dict = create_metaparameter_dict(data_dict.keys())
        return coord_dict, data_dict, meta_dict, metaparameter_dict
