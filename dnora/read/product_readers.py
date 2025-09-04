from .abstract_readers import SpectralDataReader, DataReader
from dnora.read.file_structure import FileStructure
from dnora.read.product_configuration import ProductConfiguration
from functools import partial
from dnora.type_manager.data_sources import DataSource
from dnora.read.ds_read_functions import (
    read_ds_list,
    read_first_ds,
    find_time_var_in_ds,
    setup_temp_dir,
    read_list_of_spatial_ds,
)

from geo_skeletons import PointSkeleton
from dnora import msg
import xarray as xr
from dnora import utils

from typing import Union, Optional
class ProductReader(DataReader):
    """This should not be used directly, but is a template for model specific implementations"""

    # This defines filenames, data sources etc.
    product_configuration = ProductConfiguration(default_data_source=DataSource.LOCAL)

    # This defines how the file structure of the model output is set up
    file_structure = FileStructure(
        stride="month",
        hours_per_file=None,
        last_file="",
        lead_time=0,
        offset=0,
    )

    def __init__(
        self,
        stride: Union[int, str, None] = None,  # Integer is number of hours, 'month' for monthly files
        hours_per_file: Optional[int] = None,  # None for stride = 'month'
        last_file: Optional[str] = None,
        lead_time: Optional[int] = None,
        offset: Optional[int] = None,
        tile: Optional[str] = None,
    ) -> None:
        if stride is not None:
            self.file_structure.stride = stride

        if hours_per_file is not None:
            self.file_structure.hours_per_file = hours_per_file

        if last_file is not None:
            self.file_structure.last_file = last_file

        if lead_time is not None:
            self.file_structure.lead_time = lead_time

        if offset is not None:
            self.file_structure.offset = offset

        self._default_filename = self.product_configuration.filename
        self._default_filenames = self.product_configuration.default_filenames
        self._default_folders = self.product_configuration.default_folders
        self.set_default_data_source(self.product_configuration.default_data_source)

        self._tile = tile or self.product_configuration.tile

    def __call__(
        self,
        obj_type,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        filename: str,
        expansion_factor: float = 1.2,
        program: str = "pyfimex",
        dnora_class=None,
        tile: Optional[str] = None,
        **kwargs,
    ) -> tuple[dict]:
        """Reads in all gridded data for given area and time"""
        tile = tile or self._tile
        tile_name = self.product_configuration.tile_names.get(tile)
        folder = self._folder(folder, source, tile=tile, tile_name=tile_name, strict=False)
        filename = self._filename(filename, source, tile=tile, tile_name=tile_name, strict=False)
        start_times, end_times, file_times = self.file_structure.create_time_stamps(
            start_time, end_time, kwargs.get('last_file', '')
        )
        msg.info(
            f"Getting {obj_type} from {self.name()} from {start_time} to {end_time}"
        )
        setup_temp_dir(obj_type, self.name())
        # Define area to search in
        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        lon, lat = utils.grid.expand_area(
            grid.edges("lon"), grid.edges("lat"), expansion_factor
        )

        ds_creator_function = partial(
            self.product_configuration.ds_creator_function,
            lon=lon,
            lat=lat,
            data_type=obj_type,
            name=self.name(),
            program=program,
            data_vars=self.product_configuration.data_vars,
        )
        ds_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
            url_function=self.product_configuration.url_function,
            hours_per_file=self.file_structure.hours_per_file,
            lead_time=self.file_structure.lead_time,
        )
        if not ds_list:
            msg.warning("No data found!")
            return None
        msg.info("Merging dataset together (this might take a while)...")
        time_var = self.product_configuration.time_var or find_time_var_in_ds(
            ds_list[0]
        )
        ds = xr.concat(ds_list, dim=time_var, coords="minimal")

        ds, kwargs = self.product_configuration.ds_pre_processor(ds)
        ds_aliases = self.product_configuration.ds_aliases
        core_aliases = self.product_configuration.core_aliases

        points = dnora_class.from_ds(
            ds,
            dynamic=False,
            ignore_vars=[],
            only_vars=[],
            ds_aliases=ds_aliases,
            core_aliases=core_aliases,
            time=ds.get(time_var).values,
            **kwargs,
        )

        # in case the ds_creator_function didn't crop it
        ds = points.ds().sel(**{"lat": slice(*lat), "lon": slice(*lon)})
        return ds


class SpectralProductReader(SpectralDataReader):
    """This should not be used directly, but is a template for model specific implementations"""

    # This defines filenames, data sources etc.
    product_configuration = ProductConfiguration(default_data_source=DataSource.LOCAL)

    # This defines how the file structure of the model output is set up
    file_structure = FileStructure(
        stride="month",
        hours_per_file=None,
        last_file="",
        lead_time=0,
        offset=0,
    )

    def __init__(
        self,
        stride: Union[int, str, None] = None,  # Integer is number of hours, 'month' for monthly files
        hours_per_file: Optional[int] = None,  # None for stride = 'month'
        last_file: Optional[str] = None,
        lead_time: Optional[int] = None,
        offset: Optional[int] = None,
        tile: Optional[str] = None,
    ) -> None:
        if stride is not None:
            self.file_structure.stride = stride

        if hours_per_file is not None:
            self.file_structure.hours_per_file = hours_per_file

        if last_file is not None:
            self.file_structure.last_file = last_file

        if lead_time is not None:
            self.file_structure.lead_time = lead_time

        if offset is not None:
            self.file_structure.offset = offset

        self.set_convention(self.product_configuration.convention)

        self._default_filename = self.product_configuration.filename
        self._default_filenames = self.product_configuration.default_filenames
        self._default_folders = self.product_configuration.default_folders
        self.set_default_data_source(self.product_configuration.default_data_source)

        self._tile = tile or self.product_configuration.tile

    def get_coordinates(
        self,
        grid,
        start_time,
        source: DataSource,
        folder: str,
        filename: str,
        tile: Optional[str] = None,
        **kwargs,
    ) -> dict:
        """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""
        tile = tile or self._tile
        tile_name = self.product_configuration.tile_names.get(tile)

        folder = self._folder(folder, source, tile=tile, tile_name=tile_name)
        filename = self._filename(filename, source, tile=tile, tile_name=tile_name)

        if self.file_structure.stride is None:
            # Points scattered between files with all times in each file
            ds_list = read_list_of_spatial_ds(folder, filename)
            points = None
            for ds in ds_list:
                if points is None:
                    points = PointSkeleton.from_ds(ds)
                else:
                    points = points.absorb(PointSkeleton.from_ds(ds), dim="inds")
        else:
            ds = read_first_ds(folder, filename, start_time, self.file_structure)
            points = PointSkeleton.from_ds(ds)
        return points.coord_dict()

    def __call__(
        self,
        obj_type,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        filename: str,
        inds,
        dnora_class=None,
        tile: Optional[str] = None,
        verbose: bool = False,
        add_datavars: list = None,
        **kwargs,
    ) -> tuple[dict]:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        tile = tile or self._tile
        tile_name = self.product_configuration.tile_names.get(tile)
        folder = self._folder(folder, source, tile=tile, tile_name=tile_name)
        filename = self._filename(filename, source, tile=tile, tile_name=tile_name)
        msg.info(
            f"Getting boundary spectra from {self.name()} from {start_time} to {end_time}"
        )

        if add_datavars is not None:
            for data_var in add_datavars:
                dnora_class = dnora_class.add_datavar(data_var)

        ds_aliases = self.product_configuration.ds_aliases
        core_aliases = self.product_configuration.core_aliases

        if self.file_structure.stride is None:
            # Points scattered between files with all times in each file
            ds_list = read_list_of_spatial_ds(folder, filename)
            points = None
            for ds in ds_list:
                data = dnora_class.from_ds(
                    ds,
                    ds_aliases=ds_aliases,
                    core_aliases=core_aliases,
                ).sel(time=slice(start_time, end_time))
                if points is None:
                    points = data
                else:
                    points = points.absorb(data, dim="inds")
        else:
            start_times, end_times, file_times = self.file_structure.create_time_stamps(
                start_time, end_time
            )

            ds_creator_function = partial(
                self.product_configuration.ds_creator_function, inds=inds
            )
            ds_list = read_ds_list(
                start_times,
                end_times,
                file_times,
                folder,
                filename,
                ds_creator_function,
                url_function=self.product_configuration.url_function,
                hours_per_file=self.file_structure.hours_per_file,
                lead_time=self.file_structure.lead_time,
            )

            msg.info("Merging dataset together (this might take a while)...")
            time_var = self.product_configuration.time_var or find_time_var_in_ds(
                ds_list[0]
            )
            ds = xr.concat(ds_list, dim=time_var)
            ds, kwargs = self.product_configuration.ds_pre_processor(ds)

            points = dnora_class.from_ds(
                ds,
                dynamic=False,
                ignore_vars=[],
                only_vars=[],
                ds_aliases=ds_aliases,
                core_aliases=core_aliases,
                verbose=verbose,
                **kwargs,
            )

        return points.ds()
