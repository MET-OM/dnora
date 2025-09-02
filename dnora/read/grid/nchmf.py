from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from dnora.grid import Grid, TriGrid
from dnora.read.abstract_readers import DataReader
from dnora.utils.io import get_url
from dnora import utils
from dnora import msg
from typing import Union
import os
from dnora.type_manager.data_sources import DataSource
from dnora.type_manager.dnora_types import DnoraDataType


class DepFile(DataReader):
    """Reads bathymetry from ascii .dep-files"""

    def default_data_source(self) -> DataSource:
        return DataSource.LOCAL

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid: Union[Grid, TriGrid],
        start_time,
        end_time,
        source: DataSource,
        expansion_factor: float = 1.2,
        folder: str = None,
        filename: str = None,
        **kwargs,
    ) -> tuple:

        folder = get_url(folder, f"DepFile")
        self.filename = get_url(folder, filename)
        # Area is expanded a bit to not get in trouble in the meshing stage
        # when we interpoolate or filter

        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        lon, lat = utils.grid.expand_area(
            grid.edges("lon"), grid.edges("lat"), expansion_factor
        )

        with open(self.filename, "r") as file:
            line = True
            ncols = int(file.readline()[5:])
            nrows = int(file.readline()[5:])
            lon0 = float(file.readline()[9:])
            lat0 = float(file.readline()[9:])
            dcell = float(file.readline()[8:])
            nodata = int(file.readline()[12:])
            topo = np.zeros((nrows, ncols))
            for n in range(nrows):
                topo[n, :] = np.array(
                    [int(val) for val in file.readline().split(" ") if val]
                )
        topo[topo == nodata] == 0
        lons = np.arange(lon0, lon0 + dcell * (ncols - 1), dcell)
        lats = np.arange(lat0, lat0 + dcell * (nrows - 1), dcell)

        msg.plain(
            f"Grid in file start at lon={lon0} ({ncols} points), lat={lat0} ({nrows} points) and has {dcell} resolution."
        )

        msg.plain(
            f"Reading bathymetry for: {lon[0]:10.7f}-{lon[1]:10.7f}, {lat[0]:10.7f}-{lat[1]:10.7f}."
        )

        lon_inds = np.where(np.logical_and(lons > lon[0], lons < lon[1]))[0]
        lat_inds = np.where(np.logical_and(lats > lat[0], lats < lat[1]))[0]
        topo = np.flipud(topo)
        topo = topo[lat_inds, :]
        topo = topo[:, lon_inds]

        coord_dict = {"lon": lons[lon_inds], "lat": lats[lat_inds]}

        data_dict = {"topo": -topo}
        meta_dict = {}

        return coord_dict, data_dict, meta_dict

    def __str__(self):
        return f"Reading bathymetry fron {self.filename}."
