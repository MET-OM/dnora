from __future__ import annotations
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from dnora.read.abstract_readers import DataReader
from geo_skeletons import GriddedSkeleton, PointSkeleton
from geo_skeletons.decorators import add_datavar, add_mask


from dnora.type_manager.data_sources import (
    DataSource,
    data_source_from_string,
)
from dnora import msg, utils
from dnora.defaults import read_environment_variable
import os

import geo_parameters as gp
from dnora.type_manager.dnora_types import DnoraDataType

from dnora.utils.io import get_url


@add_mask(name="sea", coord_group="grid", default_value=1, opposite_name="land")
@add_datavar(name=gp.ocean.WaterDepth("topo"), default_value=999.0, coord_group="grid")
class GriddedTopo(GriddedSkeleton):
    _default_reader = None
    pass


@add_mask(name="sea", coord_group="grid", default_value=1, opposite_name="land")
@add_datavar(name=gp.ocean.WaterDepth("topo"), default_value=999.0, coord_group="grid")
class PointTopo(PointSkeleton):
    _default_reader = None
    pass


def import_topo(
    grid,
    topo_reader: DataReader = None,
    source: Union[str, DataSource] = None,
    folder: str = None,
    **kwargs,
) -> GriddedTopo:
    """Reads the raw bathymetrical data."""
    if topo_reader is None:
        raise ValueError("Define a DataReader!")
    msg.header(topo_reader, "Importing topography...")

    source = source or topo_reader.default_data_source()
    source = data_source_from_string(source)

    if folder is None:
        folder = read_environment_variable(DnoraDataType.GRID, source)

    if folder is None:
        folder = ""

    if folder and source == DataSource.LOCAL:
        if not os.path.exists(os.path.expanduser(folder)):
            os.mkdir(folder)
        if not os.path.exists(get_url(folder, topo_reader.name())):
            os.mkdir(get_url(folder, topo_reader.name()))

    # start_time, end_time and inds given to satisfy the generic Netcdf-readers
    # Will need to fix this at some point
    ds = topo_reader(
        obj_type=DnoraDataType.GRID,
        grid=grid,
        start_time=None,
        end_time=None,
        source=source,
        folder=folder,
        inds=None,
        **kwargs,
    )

    if not isinstance(ds, tuple):
        coord_dict = {key: ds.get(key) for key in ["lon", "lat", "x", "y"]}
        data_dict = {"topo": ds.get("topo").values}
        meta_dict = ds.attrs
    else:
        coord_dict, data_dict, meta_dict = ds

    lon, lat, x, y = (
        coord_dict.get("lon"),
        coord_dict.get("lat"),
        coord_dict.get("x"),
        coord_dict.get("y"),
    )

    topo = data_dict.get("topo")
    zone_number, zone_letter = data_dict.get("zone_number"), data_dict.get(
        "zone_letter"
    )

    if 0 in topo.shape:
        msg.warning("Imported topography seems to be empty. Maybe using wrong tile?")
        return

    if utils.grid.is_gridded(topo, lon, lat) or utils.grid.is_gridded(topo, x, y):
        topo_grid = GriddedTopo(lon=lon, lat=lat, x=x, y=y, name=topo_reader.name())
        topo_grid.set_spacing(nx=len(x or lon), ny=len(y or lat))
    else:
        topo_grid = PointTopo(lon=lon, lat=lat, x=x, y=y, name=topo_reader.name())
    if zone_number is not None:
        topo_grid.utm.set((zone_number, zone_letter))

    if (
        grid.edges("lon", native=True)[0] < topo_grid.edges(grid.core.x_str)[0]
        or grid.edges("lon", native=True)[1] > topo_grid.edges(grid.core.x_str)[1]
    ):
        msg.warning(
            f"The data gotten from the DataReader doesn't cover the grid in the {grid.core.x_str} direction. Grid: {grid.edges('lon', native=True)}, imported topo: {topo_grid.edges(grid.core.x_str)}"
        )

    if (
        grid.edges("lat", native=True)[0] < topo_grid.edges(grid.core.y_str)[0]
        or grid.edges("lat", native=True)[1] > topo_grid.edges(grid.core.y_str)[1]
    ):
        msg.warning(
            f"The data gotten from the DataReader doesn't cover the grid in the {grid.core.y_str} direction. Grid: {grid.edges('lat', native=True)}, imported topo: {topo_grid.edges(grid.core.y_str)}"
        )

    topo_grid.set_topo(topo)
    topo_grid.meta.set(meta_dict)

    return topo_grid
