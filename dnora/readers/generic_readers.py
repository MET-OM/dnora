from .abstract_readers import DataReader, PointDataReader
from ..dnora_object_type import DnoraDataType
from ..data_sources import DataSource
import pandas as pd
import numpy as np


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
        **kwargs
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

        # # If not found, use the metaparameter defined in the object
        # for key in variables:
        #     if metaparameter_dict.get(key) is None:
        #         metaparameter_dict[key] = obj_type.value.meta_dict.get(key)

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
        **kwargs
    ):
        time = pd.date_range(start=start_time, end=end_time, freq="H").values

        coord_dict = {}
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

        if "time" in obj_type.value._coord_manager.added_coords():
            coord_dict["time"] = time
            obj_size = (len(time), len(inds))
        else:
            obj_size = (len(inds),)

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
