from .abstract_readers import DataReader
from ..dnora_object_type import DnoraDataType
from ..data_sources import DataSource
import pandas as pd


class ConstantGrid(DataReader):
    def __init__(self, u: float = 1, v: float = 2, metadata: dict = None):
        self.u = u
        self.v = v
        self.metadata = metadata

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid,
        start_time,
        end_time,
        source: DataSource,
        **kwargs
    ):
        time = pd.date_range(start=start_time, end=end_time, freq="H").values

        coord_dict = {}
        coord_dict["lon"] = grid.lon(strict=True)
        coord_dict["lat"] = grid.lat(strict=True)
        coord_dict["x"] = grid.x(strict=True)
        coord_dict["y"] = grid.y(strict=True)

        u = np.full((len(time), grid.ny(), grid.nx()), self.u)
        v = np.full((len(time), grid.ny(), grid.nx()), self.v)
        metadata = {"metadata": "this is a constant forcing"}

        return time, u, v, lon, lat, x, y, metadata
