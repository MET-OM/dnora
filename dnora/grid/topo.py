from geo_skeletons import GriddedSkeleton, PointSkeleton
from geo_skeletons.decorators import add_datavar, add_mask
from dnora.readers.abstract_readers import DataReader
from dnora.dnora_type_manager.data_sources import (
    DataSource,
    data_source_from_string,
)
from dnora import msg, aux_funcs
from dnora.defaults import read_environment_variable
import os

from dnora.metaparameter.parameter_funcs import set_metaparameters_in_object
from dnora.metaparameter import WaterDepth
from dnora.dnora_type_manager.dnora_types import DnoraDataType


@add_mask(name="sea", coords="grid", default_value=1, opposite_name="land")
@add_datavar(name="topo", default_value=999.0, coords="grid")
class GriddedTopo(GriddedSkeleton):
    meta_dict = {"topo": WaterDepth}
    _default_reader = None
    pass


@add_mask(name="sea", coords="grid", default_value=1, opposite_name="land")
@add_datavar(name="topo", default_value=999.0, coords="grid")
class PointTopo(PointSkeleton):
    meta_dict = {"topo": WaterDepth}
    _default_reader = None
    pass


def import_topo(
    grid,
    topo_reader: DataReader = None,
    source: str | DataSource = None,
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
        if not os.path.exists(aux_funcs.get_url(folder, topo_reader.name())):
            os.mkdir(aux_funcs.get_url(folder, topo_reader.name()))

    coord_dict, data_dict, meta_dict, metaparameter_dict = topo_reader(
        obj_type=DnoraDataType.GRID,
        grid=grid,
        start_time=None,
        end_time=None,
        source=source,
        folder=folder,
        **kwargs,
    )

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

    if aux_funcs.is_gridded(topo, lon, lat) or aux_funcs.is_gridded(topo, x, y):
        topo_grid = GriddedTopo(lon=lon, lat=lat, x=x, y=y)
        topo_grid.set_spacing(nx=len(x or lon), ny=len(y or lat))
    else:
        topo_grid = PointTopo(lon=lon, lat=lat, x=x, y=y)

    if (
        grid.edges("lon", native=True)[0] < topo_grid.edges(grid.x_str)[0]
        or grid.edges("lon", native=True)[1] > topo_grid.edges(grid.x_str)[1]
    ):
        msg.warning(
            f"The data gotten from the DataReader doesn't cover the grid in the {grid.x_str} direction. Grid: {grid.edges('lon', native=True)}, imported topo: {topo_grid.edges(grid.x_str)}"
        )

    if (
        grid.edges("lat", native=True)[0] < topo_grid.edges(grid.y_str)[0]
        or grid.edges("lat", native=True)[1] > topo_grid.edges(grid.y_str)[1]
    ):
        msg.warning(
            f"The data gotten from the DataReader doesn't cover the grid in the {grid.y_str} direction. Grid: {grid.edges('lat', native=True)}, imported topo: {topo_grid.edges(grid.y_str)}"
        )

    if zone_number is not None:
        topo_grid.set_utm((zone_number, zone_letter))

    topo_grid.set_topo(topo)
    topo_grid = set_metaparameters_in_object(topo_grid, metaparameter_dict, data_dict)
    topo_grid.set_metadata(meta_dict)

    return topo_grid
