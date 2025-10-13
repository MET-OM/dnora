from dnora.type_manager.dnora_types import DnoraDataType
from typing import Union
from dnora.read.abstract_readers import (
    DataReader,
    SpectralDataReader,
    ReaderFunction,
)
from dnora.type_manager.dnora_objects import DnoraObject, Grid, dnora_objects
import geo_parameters as gp
import pandas as pd
import numpy as np
from dnora.pick import PointPicker, Area, NearestGridPoint
from geo_skeletons import PointSkeleton
from dnora import msg


def import_data(
    grid: Grid,
    start_time: str,
    end_time: str,
    obj_type: DnoraDataType,
    name,
    dry_run,
    reader,
    expansion_factor: float,
    source,
    folder: str,
    filename: str,
    point_mask=None,
    point_picker=None,
    **kwargs,
) -> DnoraObject:
    """Imports data using DataReader and creates and returns a DNORA object"""
    msg.plain("")
    msg.print_line(marker="+")
    msg.plain(f"Starting import of data: {obj_type.name}")
    msg.print_line(marker="+")
    msg.plain(
        f"Area: {grid.core.x_str}: {grid.edges('lon',native=True)}, {grid.core.y_str}: {grid.edges('lat',native=True)}"
    )
    msg.plain(f"{start_time} - {end_time}")

    if dry_run:
        msg.info("Dry run! No data will be imported.")
        return

    if not dnora_objects.get(obj_type).is_gridded():
        msg.header(point_picker, "Choosing points to import...")
        inds = pick_points(
            grid,
            reader,
            start_time,
            point_picker,
            expansion_factor,
            point_mask,
            source,
            folder,
            filename,
            **kwargs,
        )
        if len(inds) < 1:
            msg.warning("PointPicker didn't find any points. Aborting import of data.")
            return
    else:
        inds = np.array([])

    msg.header(reader, f"Importing {obj_type.name}...")

    obj = read_data_and_create_object(
        obj_type,
        reader,
        expansion_factor,
        grid,
        start_time,
        end_time,
        name,
        source,
        folder,
        filename,
        inds,
        **kwargs,
    )
    return obj


def read_data_and_create_object(
    obj_type: DnoraDataType,
    reader: Union[DataReader, SpectralDataReader],
    expansion_factor: float,
    grid: Grid,
    start_time: str,
    end_time: str,
    name: str,
    source: DnoraDataType,
    folder: str,
    filename: str,
    inds: list[int] = None,
    **kwargs,
) -> DnoraObject:
    """Reads data using the reader, creates the objects and sets data and metadata in object"""
    dnora_class = dnora_objects.get(obj_type)

    ds = reader(
        obj_type=obj_type,
        grid=grid,
        start_time=pd.to_datetime(start_time),
        end_time=pd.to_datetime(end_time),
        source=source,
        folder=folder,
        filename=filename,
        inds=inds,
        dnora_class=dnora_class,
        expansion_factor=expansion_factor,
        **kwargs,
    )
    if ds is None:
        raise ValueError("No data found!")

    if not isinstance(ds, tuple):
        obj = dnora_class.from_ds(ds, name=name, dynamic=False)
    else:
        coord_dict, data_dict, meta_dict = ds
        obj = dnora_class(name=name, **coord_dict)
        existing_vars = obj.core.non_coord_objects()

        for key, value in data_dict.items():
            # Give (name[str], parmater[gp]) or (name[str], None) is values is a str
            name, param = gp.decode(key)

            # Only a string identifier was given, e.g. 'hs'
            if name not in existing_vars and param is not None:
                names = obj.core.find_cf(param.standard_name()) + obj.core.find_cf(
                    param.standard_name(alias=True)
                )
                names = list(set(names))  # Remove duplicates

                if len(names) > 1:
                    raise KeyError(
                        f"The standard_name '{param.standard_name()} of class {param.__name__} matches several variables: {names}'. Try specifying with e.g. {param.__name__}('{names[0]}') "
                    )
                elif len(names) < 1:
                    name = 0
                else:
                    name = names[0]

            if name in existing_vars:
                if isinstance(value, tuple):
                    obj.set(name, value[0], coords=value[1])
                else:
                    obj.set(name, value)

        obj.meta.append(meta_dict)
        if (
            meta_dict.get("zone_number") is not None
            and meta_dict.get("zone_letter") is not None
        ):
            obj.utm.set((meta_dict.get("zone_number"), meta_dict.get("zone_letter")))

    try:
        obj._mark_convention(reader.convention())
    except AttributeError:
        pass

    return obj


def pick_points(
    grid: Grid,
    reader: ReaderFunction,
    start_time: str,
    point_picker: PointPicker,
    expansion_factor: float,
    point_mask: np.ndarray[bool],
    source: str,
    folder: str,
    filename: str,
    **kwargs,
):
    """Gets the indeces of the point defined by the logical mask with respect to all points available from the reader function."""
    available_points = reader.get_coordinates(
        grid=grid,
        start_time=start_time,
        source=source,
        folder=folder,
        filename=filename,
        **kwargs,
    )

    all_points = PointSkeleton(
        lon=available_points.get("lon"),
        lat=available_points.get("lat"),
        x=available_points.get("x"),
        y=available_points.get("y"),
    )

    ## Only take points that are reasonable close to the wanted grid
    ## This speeds up the searcg considerably, especially if we have points over 84 lat
    ## since then the fast cartesian searhc is not possible
    slon, slat = grid.edges("lon"), grid.edges("lat")
    search_grid = Grid(lat=(slat[0] - 3, slat[1] + 3), lon=(slon[0] - 6, slon[1] + 6))
    if isinstance(point_picker, NearestGridPoint):
        search_inds = Area()(search_grid, all_points, expansion_factor=1)
    else:
        search_inds = all_points.inds()

    # if np.all(np.logical_not(point_mask)):
    #     msg.warning(
    #         "None of the points set to interest points! Aborting import of data."
    #     )
    #     return
    # else:
    #     interest_points = PointSkeleton.from_skeleton(grid, mask=point_mask)

    if not np.all(np.logical_not(point_mask)):
        interest_points = PointSkeleton.from_skeleton(grid, mask=point_mask)
    else:
        interest_points = None

    inds = point_picker(
        grid=grid,
        all_points=all_points.sel(inds=search_inds),
        expansion_factor=expansion_factor,
        selected_points=interest_points,
        fast=True,
        **kwargs,
    )

    if len(inds) < 1:
        return np.array([])
    return search_inds[inds]
