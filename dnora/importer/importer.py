from dnora.readers.abstract_readers import (
    ReaderFunction,
    DataReader,
    SpectralDataReader,
)
from dnora.dnora_type_manager.dnora_types import DnoraDataType
from dnora import msg
from dnora.pick import PointPicker
from geo_skeletons import PointSkeleton
from geo_skeletons.decorators import add_datavar
from dnora.dnora_type_manager.dnora_objects import DnoraObject, Grid, dnora_objects
import numpy as np
from dnora.aux_funcs import set_metaparameters_in_object
import pandas as pd
from typing import Union


class DataImporter:

    def _pick_points(
        self,
        grid: Grid,
        reader: ReaderFunction,
        start_time: str,
        point_picker: PointPicker,
        point_mask: np.ndarray[bool],
        source: str,
        folder: str,
        **kwargs,
    ):
        """Gets the indeces of the point defined by the logical mask with respect to all points available from the reader function."""
        available_points = reader.get_coordinates(
            grid=grid,
            start_time=start_time,
            source=source,
            folder=folder,
            **kwargs,
        )

        all_points = PointSkeleton(
            lon=available_points.get("lon"),
            lat=available_points.get("lat"),
            x=available_points.get("x"),
            y=available_points.get("y"),
        )

        if np.all(np.logical_not(point_mask)):
            msg.warning(
                "None of the points set to interest points! Aborting import of data."
            )
            return
        else:
            interest_points = PointSkeleton.from_skeleton(grid, mask=point_mask)

        inds = point_picker(
            grid=grid,
            all_points=all_points,
            selected_points=interest_points,
            fast=True,
            **kwargs,
        )

        if len(inds) < 1:
            msg.warning("PointPicker didn't find any points. Aborting import of data.")
            return
        return inds

    @staticmethod
    def _read_data_and_create_object(
        obj_type: DnoraDataType,
        reader: Union[DataReader, SpectralDataReader],
        grid: Grid,
        start_time: str,
        end_time: str,
        name: str,
        source: DnoraDataType,
        folder: str,
        inds: list[int] = None,
        **kwargs,
    ) -> DnoraObject:
        """Reads data using the reader, creates the objects and sets data and metadata in object"""

        coord_dict, data_dict, meta_dict, metaparameter_dict = reader(
            obj_type=obj_type,
            grid=grid,
            start_time=pd.to_datetime(start_time),
            end_time=pd.to_datetime(end_time),
            source=source,
            folder=folder,
            inds=inds,
            **kwargs,
        )

        obj = dnora_objects.get(obj_type)(name=name, **coord_dict)

        for key, value in data_dict.items():
            if obj.get(key) is None:
                obj = add_datavar(key, append=True)(obj)  # Creates .hs() etc. methods
            obj.set(key, value, allow_reshape=True)

        obj = set_metaparameters_in_object(obj, metaparameter_dict, data_dict)

        obj.set_metadata(meta_dict)
        if (
            meta_dict.get("zone_number") is not None
            and meta_dict.get("zone_letter") is not None
        ):
            obj.set_utm((meta_dict.get("zone_number"), meta_dict.get("zone_letter")))

        if obj_type in [DnoraDataType.SPECTRA, DnoraDataType.SPECTRA1D]:
            obj._mark_convention(reader.convention())

        return obj

    def import_data(
        self,
        grid: Grid,
        start_time: str,
        end_time: str,
        obj_type: DnoraDataType,
        name,
        dry_run,
        reader,
        source,
        folder,
        point_mask=None,
        point_picker=None,
        **kwargs,
    ) -> DnoraObject:
        """Imports data using DataReader and creates and returns a DNORA object"""

        msg.header(point_picker, "Choosing points to import...")
        if dry_run:
            msg.info("Dry run! No data will be imported.")
        else:
            if point_mask is not None:

                inds = self._pick_points(
                    grid,
                    reader,
                    start_time,
                    point_picker,
                    point_mask,
                    source,
                    folder,
                    **kwargs,
                )
            else:
                inds = None

        msg.header(reader, f"Importing {obj_type.name}...")
        msg.plain(
            f"Area: {grid.x_str}: {grid.edges('lon',native=True)}, {grid.y_str}: {grid.edges('lat',native=True)}"
        )
        msg.plain(f"{start_time} - {end_time}")

        if dry_run:
            msg.info("Dry run! No data will be imported.")
            return

        obj = self._read_data_and_create_object(
            obj_type,
            reader,
            grid,
            start_time,
            end_time,
            name,
            source,
            folder,
            inds,
            **kwargs,
        )
        return obj
