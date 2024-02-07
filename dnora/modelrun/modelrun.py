from __future__ import annotations  # For TYPE_CHECKING
from copy import copy
import pandas as pd
import numpy as np
from geo_skeletons import PointSkeleton
from typing import Union

# Import objects
from dnora.grid import Grid, TriGrid

from dnora.file_module import FileNames
from dnora.data_sources import DataSource

# Import abstract classes and needed instances of them
from dnora.pick.point_pickers import PointPicker, NearestGridPoint


from dnora.spectral_grid import SpectralGrid

from geo_skeletons.decorators import add_datavar

from dnora import msg
from dnora.cacher.cache_decorator import cached_reader

from pathlib import Path
from dnora.defaults.default_reader import data_sources

from dnora.spectra1d.read import SpectraTo1D
from dnora.waveseries.read import SpectraToWaveSeries

from dnora.pick import TrivialPicker

from dnora.export.templates import Cacher

from dnora.readers import generic_readers
from dnora.readers.abstract_readers import (
    DataReader,
    PointDataReader,
    SpectralDataReader,
)

from dnora.dnora_types import (
    ReaderFunction,
    DnoraDataType,
    DnoraObject,
    data_type_from_string,
)
from dnora.dnora_types import (
    Grid,
    Wind,
    Spectra,
    Spectra1D,
    WaterLevel,
    WaveSeries,
    Current,
    Ice,
)


class ModelRun:
    _reader_dict: dict[DnoraDataType:ReaderFunction] = {}
    _point_picker: PointPicker = NearestGridPoint()

    def __init__(
        self,
        grid: Grid | None = None,
        start_time: str = "1970-01-01T00:00",
        end_time: str = "2030-12-31T23:59",
        name: str = "AnonymousModelRun",
        dry_run: bool = False,
    ):
        self.name = copy(name)
        self._grid = copy(grid)
        self._time = pd.date_range(start_time, end_time, freq="H")
        self._exported_to: dict[DnoraDataType : list[str]] = {}
        self._global_dry_run = dry_run
        self._dry_run = False  # Set by methods

        self._dnora_objects: dict[DnoraDataType:DnoraObject] = {
            DnoraDataType.GRID: grid,
        }
        self._consistency_check(
            objects_to_ignore_get=[
                DnoraDataType.TRIGRID,
            ],
            objects_to_ignore_import=[
                DnoraDataType.GRID,
                DnoraDataType.TRIGRID,
            ],
        )

    def _consistency_check(
        self, objects_to_ignore_get: list[str], objects_to_ignore_import: list[str]
    ):
        """Checks that the class contains the proper import and getter methods. This is a safety feature in case more object types are added."""
        for obj_type in DnoraDataType:
            if obj_type not in objects_to_ignore_get:
                if not hasattr(self, obj_type.name.lower()):
                    raise SyntaxError(
                        f"No getter method self.{obj_type.name.lower()}() defined for object {obj_type.name}!"
                    )

            if obj_type not in objects_to_ignore_import:
                if not hasattr(self, f"import_{obj_type.name.lower()}"):
                    raise SyntaxError(
                        f"No import method self.import_{obj_type.name.lower()}() defined for object {obj_type.name}!"
                    )

    def _setup_import(
        self,
        obj_type: DnoraDataType,
        name: str,
        dry_run: bool,
        reader: ReaderFunction,
        source: DataSource | str,
        folder: str,
    ) -> tuple[ReaderFunction, str]:
        """Sets up readers, names and dry runs porperties for import of object."""

        self._dry_run = dry_run

        reader = reader or self._get_reader(obj_type)
        if reader is None:
            raise Exception(f"Define a {obj_type.name}Reader!")

        name = name or reader.name()

        if name is None:
            raise ValueError(
                f"Provide either a name or a {obj_type.name}Reader that will then define the name!"
            )

        msg.header(reader, f"Importing {obj_type.name}...")

        if isinstance(source, str):
            try:
                source = DataSource[source.upper()]
            except KeyError as ke:
                raise ke(
                    f"source should be 'local' (DataSource.LOCAL), 'internal' (DataSource.INTERNAL), 'remote' (DataSource.REMOTE),  or 'undefined' (DataSource.UNDEFINED), not {source}!"
                )

        if folder and source is DataSource.UNDEFINED:
            source = DataSource.LOCAL

        if source in [DataSource.INTERNAL, DataSource.LOCAL]:
            if folder is None:
                folder = data_sources(source)
        elif source == DataSource.REMOTE:
            folder = ""

        if source == DataSource.UNDEFINED:
            source = reader.default_data_source()

        return reader, name, source, folder

    def _setup_point_picker(self, point_picker: PointPicker):
        """Sets up point picker using possible default values."""
        point_picker = point_picker or self._get_point_picker()
        if point_picker is None:
            raise Exception("Define a PointPicker!")
        return point_picker

    def _pick_points(
        self,
        reader: ReaderFunction,
        point_picker: PointPicker,
        mask_of_points: np.ndarray[bool],
        source: str,
        folder: str,
        **kwargs,
    ):
        """Gets the indeces of the point defined by the logical mask with respect to all points available from the reader function."""
        if not self.dry_run():
            available_points = reader.get_coordinates(
                grid=self.grid(),
                start_time=self.start_time(),
                source=source,
                folder=folder,
            )

            all_points = PointSkeleton(
                lon=available_points.get("lon"),
                lat=available_points.get("lat"),
                x=available_points.get("x"),
                y=available_points.get("y"),
            )

            if np.all(np.logical_not(mask_of_points)):
                interest_points = None
            else:
                interest_points = PointSkeleton.from_skeleton(
                    self.grid(), mask=mask_of_points
                )

            msg.header(point_picker, "Choosing points to import...")
            inds = point_picker(
                grid=self.grid(),
                all_points=all_points,
                selected_points=interest_points,
                fast=True,
                **kwargs,
            )

            if len(inds) < 1:
                msg.warning(
                    "PointPicker didn't find any points. Aborting import of data."
                )
                return
            return inds

    def _read_data(
        self, obj_type: DnoraDataType, reader, name, source, folder, inds, **kwargs
    ):
        """Reads data using the reader, creates the objects and sets data and metadata in object"""
        if self.dry_run():
            msg.info("Dry run! No forcing will be imported.")
            return

        coord_dict, data_dict, meta_dict, metaparameter_dict = reader(
            obj_type=obj_type,
            grid=self.grid(),
            start_time=self.start_time(),
            end_time=self.end_time(),
            source=source,
            folder=folder,
            inds=inds,
            **kwargs,
        )

        obj = obj_type.value(name=name, **coord_dict)
        for key, value in data_dict.items():
            if obj.get(key) is None:
                obj = add_datavar(key, append=True)(obj)  # Creates .hs() etc. methods
            obj.set(key, value)

            metaparameter = metaparameter_dict.get(
                key
            )  # Check if metaparameter provided by reader
            if metaparameter is None:
                # DNORA object usually has specified the metaparameters
                if hasattr(obj_type.value, "meta_dict"):
                    metaparameter = obj_type.value.meta_dict.get(key)

            if metaparameter is not None:
                obj.set_metadata(metaparameter.meta_dict(), data_array_name=key)

        obj.set_metadata(meta_dict)

        if obj_type in [DnoraDataType.SPECTRA, DnoraDataType.SPECTRA1D]:
            obj._mark_convention(reader.convention())

        self[obj_type] = obj

    def _import_data(
        self,
        obj_type: DnoraDataType,
        name,
        dry_run,
        reader,
        source,
        folder,
        mask=None,
        point_picker=None,
        **kwargs,
    ):
        """Performs import of DNORA object"""
        reader, name, source, folder = self._setup_import(
            obj_type, name, dry_run, reader, source, folder
        )

        if mask is not None:
            point_picker = self._setup_point_picker(point_picker)

            inds = self._pick_points(
                reader,
                point_picker,
                mask,
                source,
                folder,
                **kwargs,
            )
        else:
            inds = None

        self._read_data(obj_type, reader, name, source, folder, inds, **kwargs)

    def cache(self, obj_type: DnoraDataType | str) -> None:
        Cacher(self).export(obj_type)

    @cached_reader(DnoraDataType.WIND, generic_readers.Netcdf)
    def import_wind(
        self,
        reader: DataReader | None = None,
        name: str | None = None,
        dry_run: bool = False,
        source: str | DataSource = DataSource.UNDEFINED,
        folder: str = None,
        **kwargs,
    ) -> None:
        """Import wind data from a source using the given reader"""
        self._import_data(
            DnoraDataType.WIND, name, dry_run, reader, source, folder, **kwargs
        )

    @cached_reader(DnoraDataType.WATERLEVEL, generic_readers.Netcdf)
    def import_waterlevel(
        self,
        reader: DataReader | None = None,
        name: str | None = None,
        dry_run: bool = False,
        source: str | DataSource = DataSource.UNDEFINED,
        folder: str = None,
        **kwargs,
    ) -> None:
        """Import waterlevel data from a source using the given reader"""
        self._import_data(
            DnoraDataType.WATERLEVEL, name, dry_run, reader, source, folder, **kwargs
        )

    @cached_reader(DnoraDataType.SPECTRA, generic_readers.Netcdf)
    def import_spectra(
        self,
        reader: SpectralDataReader | None = None,
        point_picker: PointPicker | None = None,
        name: str | None = None,
        dry_run: bool = False,
        source: str | DataSource = DataSource.UNDEFINED,
        folder: str = None,
        **kwargs,
    ) -> None:
        self._import_data(
            DnoraDataType.SPECTRA,
            name,
            dry_run,
            reader,
            source,
            folder,
            mask=self.grid().sea_mask(),
            point_picker=point_picker,
            **kwargs,
        )

    @cached_reader(DnoraDataType.SPECTRA1D, generic_readers.Netcdf)
    def import_spectra1d(
        self,
        reader: SpectralDataReader | None = None,
        point_picker: PointPicker | None = None,
        name: str | None = None,
        dry_run: bool = False,
        source: str | DataSource = DataSource.UNDEFINED,
        folder: str = None,
        **kwargs,
    ) -> None:
        self._import_data(
            DnoraDataType.SPECTRA1D,
            name,
            dry_run,
            reader,
            source,
            folder,
            mask=self.grid().sea_mask(),
            point_picker=point_picker,
            **kwargs,
        )

    @cached_reader(DnoraDataType.WAVESERIES, generic_readers.Netcdf)
    def import_waveseries(
        self,
        reader: PointDataReader | None = None,
        point_picker: PointPicker | None = None,
        name: str | None = None,
        dry_run: bool = False,
        source: str | DataSource = DataSource.UNDEFINED,
        folder: str = None,
        **kwargs,
    ) -> None:
        self._import_data(
            DnoraDataType.WAVESERIES,
            name,
            dry_run,
            reader,
            source,
            folder,
            mask=self.grid().sea_mask(),
            point_picker=point_picker,
            **kwargs,
        )

    @cached_reader(DnoraDataType.CURRENT, generic_readers.Netcdf)
    def import_current(
        self,
        reader: DataReader | None = None,
        name: str | None = None,
        dry_run: bool = False,
        source: str | DataSource = DataSource.UNDEFINED,
        folder: str = None,
        **kwargs,
    ) -> None:
        self._import_data(
            DnoraDataType.CURRENT, name, dry_run, reader, source, folder, **kwargs
        )

    @cached_reader(DnoraDataType.ICE, generic_readers.Netcdf)
    def import_ice(
        self,
        reader: DataReader | None = None,
        name: str | None = None,
        dry_run: bool = False,
        source: str | DataSource = DataSource.UNDEFINED,
        folder: str = None,
        **kwargs,
    ) -> None:
        self._import_data(
            DnoraDataType.ICE, name, dry_run, reader, source, folder, **kwargs
        )

    def spectra_to_1d(
        self,
        dry_run: bool = False,
        name: str | None = None,
        write_cache=False,
        **kwargs,
    ):
        self._dry_run = dry_run
        if self.spectra() is None:
            msg.warning("No Spectra to convert to Spectra!")
            return

        spectral_reader = SpectraTo1D(self.spectra())
        msg.header(
            spectral_reader,
            "Converting the boundary spectra to omnidirectional spectra...",
        )

        name = self.spectra().name

        if self.dry_run():
            msg.info("Dry run! No boundary will not be converted to spectra.")
            return

        self.import_spectra1d(
            reader=spectral_reader,
            point_picker=TrivialPicker(),
            name=name,
            write_cache=write_cache,
            **kwargs,
        )

    def spectra_to_waveseries(
        self,
        dry_run: bool = False,
        write_cache=False,
        freq: tuple = (0, 10_000),
        **kwargs,
    ):
        self._dry_run = dry_run
        if self.spectra1d() is None:
            if self.spectra() is not None:
                self.spectra_to_1d(dry_run=dry_run, write_cache=write_cache, **kwargs)
            else:
                msg.warning("No Spectra to convert to WaveSeries!")
                return

        name = self.spectra1d().name

        if self.dry_run():
            msg.info("Dry run! No boundary will not be converted to spectra.")
            return

        waveseries_reader = SpectraToWaveSeries(self.spectra1d(), freq)
        msg.header(waveseries_reader, "Converting the spectra to wave series data...")

        self.import_waveseries(
            reader=waveseries_reader,
            point_picker=TrivialPicker(),
            name=name,
            write_cache=write_cache,
            **kwargs,
        )

    def set_spectral_grid_from_spectra(self, **kwargs):
        if self.spectra() is None:
            msg.warning("No Spectra exists. Can't set spectral grid.")
            return
        self.set_spectral_grid(freq=self.spectra().freq(), **kwargs)

    def set_spectral_grid(
        self,
        freq: np.ndarray | None = None,
        dirs: np.ndarray | None = None,
        freq0: float = 0.04118,
        nfreq: int = 32,
        ndir: int = 36,
        finc: float = 1.1,
        dirshift: float | None = None,
    ):
        """Sets spectral grid for model run. Will be used to write input files."""
        if freq is None:
            freq = np.array([freq0 * finc**n for n in np.linspace(0, nfreq - 1, nfreq)])
        if len(freq) < nfreq:
            start_freq = freq[-1]
            add_freq = np.array(
                [
                    start_freq * finc**n
                    for n in np.linspace(1, nfreq - len(freq), nfreq - len(freq))
                ]
            )
            freq = np.concatenate((freq, add_freq))
        if dirs is None:
            if dirshift is None:
                dirshift = 0.0
            dirs = np.linspace(0, 360, ndir + 1)[0:-1] + dirshift
        if len(dirs) < ndir:
            if dirshift is None:
                dirshift = np.min(np.abs(dirs))
            dirs = np.linspace(0, 360, ndir + 1)[0:-1] + dirshift

        self[DnoraDataType.SpectralGrid] = SpectralGrid(
            name=DnoraDataType.SpectralGrid.value, freq=freq, dirs=dirs
        )

    def dry_run(self):
        """Checks if method or global ModelRun dryrun is True."""
        return self._dry_run or self._global_dry_run

    def grid(self) -> Union[Grid, TriGrid]:
        """Returns the grid object."""
        return self._dnora_objects.get(DnoraDataType.GRID)

    def wind(self) -> Wind:
        """Returns the forcing object if exists."""
        return self._dnora_objects.get(DnoraDataType.WIND)

    def spectra(self) -> Spectra:
        """Returns the boundary object if exists."""
        return self._dnora_objects.get(DnoraDataType.SPECTRA)

    def spectra1d(self) -> Spectra1D:
        """Returns the spectral object if exists."""
        return self._dnora_objects.get(DnoraDataType.SPECTRA1D)

    def waveseries(self) -> WaveSeries:
        """Returns the wave series object if exists."""
        return self._dnora_objects.get(DnoraDataType.WAVESERIES)

    def waterlevel(self) -> WaterLevel:
        """Returns the water level object if exists."""
        return self._dnora_objects.get(DnoraDataType.WATERLEVEL)

    def current(self) -> Current:
        """Returns the ocean current object if exists."""
        return self._dnora_objects.get(DnoraDataType.CURRENT)

    def ice(self) -> Ice:
        """Returns the ocean current object if exists."""
        return self._dnora_objects.get(DnoraDataType.ICE)

    # def spectral_grid(self) -> SpectralGrid:
    #     """Returns the spectral grid object if exists."""
    #     return self._dnora_objects.get(DnoraDataType.SpectralGrid)

    # def input_file(self) -> None:
    #     """Only defined to have method for all objects"""
    #     return None

    def list_of_objects(
        self,
    ) -> list[DnoraObject]:
        """[ModelRun, Boundary] etc."""
        return [self[x] for x in DnoraDataType if self[x] is not None]

    def dict_of_object_names(self) -> dict[str:str]:
        """{'Boundary': 'NORA3'} etc."""
        dict_of_object_names = {}
        for obj_type in DnoraDataType:
            if self[obj_type] is not None:
                dict_of_object_names[obj_type] = self[obj_type].name
        return dict_of_object_names

    def exported_to(self, obj_type: DnoraDataType) -> str:
        """Returns the path the object (e.g. grid) was exported to.

        If object has not been exported, the default filename is returned as
        a best guess
        """

        if self[obj_type] is None:
            return [""]

        default_name = FileNames(
            model=self, format=self._get_default_format(), obj_type=obj_type
        ).get_filepath()
        return self._exported_to.get(obj_type, default_name)

    def time(self, crop_with: list[DnoraDataType | str] = None):
        """Returns times of ModelRun
        crop_with = ['Forcing', 'Boundary'] gives time period covered by those objects
        crop_with = 'all' crops with all objects"""
        t0 = self._time[0]
        t1 = self._time[-1]

        if crop_with is not None:
            if type(crop_with) is not list:
                crop_with = [crop_with]

            if "all" in crop_with:
                crop_with = self.list_of_object_strings()

            for obj_str in crop_with:
                dnora_obj = self[data_type_from_string(obj_str)]
                if dnora_obj is not None:
                    time = dnora_obj.time()
                    if time[0] is not None:
                        t0 = pd.to_datetime([t0, time[0]]).max()
                    if time[-1] is not None:
                        t1 = pd.to_datetime([t1, time[-1]]).min()
        time = pd.date_range(t0, t1, freq="H")
        return time[:: len(time) - 1]

    def start_time(self, crop_with: list[str] = None):
        """Returns start time of ModelRun
        crop = True: Give the period that is covered by all objects (Forcing etc.)"""
        return self.time(crop_with=crop_with)[0]

    def end_time(self, crop_with: list[str] = None):
        """Returns start time of ModelRun
        crop = True: Give the period that is covered by all objects (Forcing etc.)"""
        return self.time(crop_with=crop_with)[-1]

    def __getitem__(self, obj_type: DnoraDataType) -> DnoraObject:
        try:
            obj_type = data_type_from_string(obj_type)
        except:
            pass

        if obj_type is None:
            return self
        if obj_type == DnoraDataType.TRIGRID:
            if type(self.grid()) == TriGrid:
                return self.grid()
            else:
                return None

        return self._dnora_objects.get(obj_type)

    def __setitem__(self, key: DnoraDataType, value: DnoraObject) -> None:
        self._dnora_objects[key] = value

    def _get_reader(self, obj_type: DnoraDataType) -> ReaderFunction:
        return self._reader_dict.get(obj_type)

    def _get_point_picker(self) -> PointPicker:
        return self._point_picker


# def camel_to_snake(string: str) -> str:
#     return re.sub(r"(?<!^)(?=[A-Z])", "_", string).lower()
