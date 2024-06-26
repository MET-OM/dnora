from __future__ import annotations  # For TYPE_CHECKING
from copy import copy
import pandas as pd
import numpy as np

from typing import Union, TYPE_CHECKING
import os


# Import objects
from dnora.grid import Grid, TriGrid

from dnora.file_module import FileNames

# Import abstract classes and needed instances of them
from dnora.pick.point_pickers import PointPicker, NearestGridPoint
from dnora.importer import DataImporter

from dnora.spectral_grid import SpectralGrid


from dnora import msg
from dnora.cacher.cache_decorator import cached_reader

from dnora.defaults import read_environment_variable
from dnora.spectra1d.read import SpectraTo1D
from dnora.waveseries.read import Spectra1DToWaveSeries
from dnora.spectral_conventions import SpectralConvention
from dnora.pick import TrivialPicker

from dnora.export.templates import Cacher
from dnora.aux_funcs import get_url, get_first_file
from dnora.readers import generic_readers
from dnora.readers.abstract_readers import (
    DataReader,
    PointDataReader,
    SpectralDataReader,
)

from dnora.dnora_type_manager.dnora_types import (
    DnoraDataType,
    DnoraFileType,
    data_type_from_string,
    file_type_from_string,
)
from dnora.dnora_type_manager.data_sources import DataSource
from calendar import monthrange

if TYPE_CHECKING:
    from dnora.dnora_type_manager.dnora_objects import (
        Grid,
        Wind,
        Spectra,
        Spectra1D,
        WaterLevel,
        WaveSeries,
        Current,
        Ice,
        DnoraObject,
    )

    from dnora.readers.abstract_readers import ReaderFunction

from typing import Optional


def start_and_end_time_of_run(
    start_time: Optional[str],
    end_time: Optional[str],
    year: Optional[int],
    month: Optional[int],
    day: Optional[int],
    hotstart_hour: bool,
) -> tuple[str, str]:
    """Determined the start and end time of the model run:

    If year is given, start is 01 January 00:00 and end is 31 December 23:00
    If month is given, start if 01 of month to last of month

    If hotstart_hour = True, then
    year is 01 January 00:00 to 01 January next year 00:00
    month is 01 month 00:00 to 01 next month 00:00"""
    if year is not None:
        year = int(year)
        if month is None:
            start_time = f"{year}-01-01 00:00"
            end_time = f"{year}-12-31 23:00"
        else:
            if day is None:
                start_time = f"{year}-{month}-01 00:00"
                nofdays = monthrange(year, month)[1]
                end_time = f"{year}-{month}-{nofdays} 23:00"
            else:
                start_time = f"{year}-{month}-{day} 00:00"
                end_time = f"{year}-{month}-{day} 23:00"

        if hotstart_hour:
            end_time = pd.to_datetime(end_time) + pd.Timedelta(hours=1)

    if start_time is None:
        start_time = pd.Timestamp.now().strftime("%Y-%m-%d %H:00")
    if end_time is None:
        end_time = (pd.to_datetime(start_time) + pd.Timedelta(hours=240)).strftime(
            "%Y-%m-%d %H:00"
        )
    return start_time, end_time


class ModelRun:
    _reader_dict: dict[DnoraDataType:ReaderFunction] = {}
    _point_picker: PointPicker = NearestGridPoint()

    def __init__(
        self,
        grid: Grid | None = None,
        start_time: str = None,
        end_time: str = None,
        year: int = None,
        month: int = None,
        day: int = None,
        hotstart_hour: bool = False,
        dry_run: bool = False,
        name: str = "DnoraModelRun",
    ):
        self._grid = grid
        start_time, end_time = start_and_end_time_of_run(
            start_time, end_time, year, month, day, hotstart_hour
        )
        self._time = pd.date_range(start_time, end_time, freq="h")
        self._data_exported_to: dict[DnoraDataType, list[str]] = {}
        self._input_file_exported_to: dict[DnoraFileType, list[str]] = {}
        self._global_dry_run = dry_run
        self._dry_run = False  # Set by methods
        self._source = DataSource.UNDEFINED
        self._reference_time = None
        self.name = name
        self._dnora_objects: dict[DnoraDataType, DnoraObject] = {
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

        if isinstance(source, str):
            try:
                source = DataSource[source.upper()]
            except KeyError as ke:
                raise ke(
                    f"source should be 'local' (DataSource.LOCAL), 'internal' (DataSource.INTERNAL), 'remote' (DataSource.REMOTE),  or 'undefined' (DataSource.UNDEFINED), not {source}!"
                )

        if folder and source is DataSource.UNDEFINED:
            source = DataSource.LOCAL

        if source == DataSource.UNDEFINED:
            source = self._source  # Internal mode might have been activated

        if source == DataSource.UNDEFINED:
            source = reader.default_data_source()

        if source in [DataSource.INTERNAL, DataSource.LOCAL]:
            if folder is None:
                folder = read_environment_variable(
                    obj_type=obj_type, data_source=source
                )
        elif source == DataSource.REMOTE:
            folder = None

        folder = folder or ""

        if folder and source == DataSource.LOCAL:
            if not os.path.exists(os.path.expanduser(folder)):
                os.mkdir(folder)
            if not os.path.exists(get_url(folder, reader.name())):
                os.mkdir(get_url(folder, reader.name()))
        return reader, name, source, folder

    def _setup_point_picker(self, point_picker: PointPicker):
        """Sets up point picker using possible default values."""
        point_picker = point_picker or self._get_point_picker()
        if point_picker is None:
            raise Exception("Define a PointPicker!")
        return point_picker

    def _import_data(
        self,
        obj_type: DnoraDataType,
        name,
        dry_run,
        reader,
        source,
        folder,
        point_mask=None,
        point_picker=None,
        **kwargs,
    ):
        """Performs import and returns DNORA object"""
        reader, name, source, folder = self._setup_import(
            obj_type, name, dry_run, reader, source, folder
        )

        point_picker = self._setup_point_picker(point_picker)
        point_mask = point_mask or self.grid().sea_mask()
        if self.forecast_mode():
            if hasattr(reader, "hours_per_file"):
                start_time = self._reference_time
                kwargs["lead_time"] = kwargs.get("lead_time", self._lead_time)
                kwargs["last_file"] = kwargs.get(
                    "last_file",
                    get_first_file(
                        start_time, reader.stride, lead_time=kwargs.get("lead_time")
                    ),
                )
                kwargs["lead_time"] = kwargs.get("lead_time", self._lead_time)
                end_time = min(
                    kwargs.get("last_file")
                    + pd.Timedelta(hours=reader.hours_per_file - 1),
                    self.end_time(),
                )
        else:
            start_time, end_time = self.start_time(), self.end_time()

        data_importer = DataImporter()
        obj = data_importer.import_data(
            grid=self.grid(),
            start_time=start_time,
            end_time=end_time,
            obj_type=obj_type,
            name=name,
            dry_run=self.dry_run(),
            reader=reader,
            source=source,
            folder=folder,
            point_picker=point_picker,
            point_mask=point_mask,
            **kwargs,
        )
        self[obj_type] = obj

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

    @cached_reader(DnoraDataType.SPECTRA, generic_readers.PointNetcdf)
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
            point_mask=self.grid().sea_mask(),
            point_picker=point_picker,
            **kwargs,
        )

    @cached_reader(DnoraDataType.SPECTRA1D, generic_readers.PointNetcdf)
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
            point_mask=self.grid().sea_mask(),
            point_picker=point_picker,
            **kwargs,
        )

    @cached_reader(DnoraDataType.WAVESERIES, generic_readers.PointNetcdf)
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
            point_mask=self.grid().sea_mask(squeeze=False),
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
        **kwargs,
    ):
        if self.spectra() is None:
            msg.warning("No Spectra to convert to Spectra!")
            return

        spectral_reader = SpectraTo1D(self.spectra())

        name = self.spectra().name

        self.import_spectra1d(
            reader=spectral_reader,
            point_picker=TrivialPicker(),
            name=name,
            dry_run=dry_run,
            **kwargs,
        )

    def spectra_to_waveseries(
        self,
        dry_run: bool = False,
        freq: tuple = (0, 10_000),
        **kwargs,
    ):
        if self.spectra1d() is None:
            if self.spectra() is not None:
                self.spectra_to_1d(dry_run=dry_run, **kwargs)
            else:
                msg.warning("No Spectra to convert to WaveSeries!")
                return

        name = self.spectra1d().name
        self.spectra1d().set_convention(SpectralConvention.MET)
        waveseries_reader = Spectra1DToWaveSeries(self.spectra1d(), freq)

        self.import_waveseries(
            reader=waveseries_reader,
            point_picker=TrivialPicker(),
            name=name,
            dry_run=dry_run,
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

    def data_exported_to(self, obj_type: DnoraDataType | str) -> str:
        """Returns the path the object (e.g. grid) was exported to.

        If object has not been exported, the default filename is returned as
        a best guess
        """
        obj_type = data_type_from_string(obj_type)
        if self[obj_type] is None:
            return [""]

        default_name = FileNames(
            model=self, format=self._get_default_format(), obj_type=obj_type
        ).get_filepath()
        return self._data_exported_to.get(obj_type, default_name)

    def input_file_exported_to(self, file_type: DnoraFileType | str) -> str:
        """Returns the path the object (e.g. grid) was exported to.

        If object has not been exported, the default filename is returned as
        a best guess
        """
        file_type = file_type_from_string(file_type)
        if self[file_type] is None:
            return [""]

        default_name = FileNames(
            model=self, format=self._get_default_format(), obj_type=file_type
        ).get_filepath()
        return self._input_file_exported_to.get(file_type, default_name)

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
                crop_with = DnoraDataType

            for obj_str in crop_with:
                dnora_obj = self[data_type_from_string(obj_str)]
                if dnora_obj is not None:
                    time = dnora_obj.time()
                    if time[0] is not None:
                        t0 = pd.to_datetime([t0, time[0]]).max()
                    if time[-1] is not None:
                        t1 = pd.to_datetime([t1, time[-1]]).min()
        time = pd.date_range(t0, t1, freq="h")
        return time
        # return time[:: len(time) - 1]

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

    def _get_reader(self, obj_type: DnoraDataType | str) -> ReaderFunction:
        obj_type = data_type_from_string(obj_type)
        return self._reader_dict.get(obj_type)

    def _get_point_picker(self) -> PointPicker:
        return self._point_picker

    def activate_internal_mode(self, folder: str = None) -> None:
        self._activate_source_mode(DataSource.INTERNAL, folder)

    def activate_remote_mode(self) -> None:
        self._activate_source_mode(DataSource.REMOTE, folder=None)

    def activate_local_mode(self, folder: str = None) -> None:
        self._activate_source_mode(DataSource.LOCAL, folder)

    def _activate_source_mode(self, source: DataSource, folder: str):
        self._source = source
        if folder is not None:
            os.environ[f"DNORA_{source.name}_PATH"] = folder

    def deactivate_source_mode(self) -> None:
        self._source = DataSource.UNDEFINED

    def activate_forecast_mode(
        self, reference_time: str = None, lead_time: int = 0
    ) -> None:
        reference_time = reference_time or self.start_time()
        self._reference_time = pd.to_datetime(reference_time)
        self._lead_time = lead_time
        if lead_time > 0:
            self._time = pd.date_range(
                reference_time,
                pd.to_datetime(reference_time) + pd.Timedelta(hours=lead_time),
                freq="h",
            )
        msg.info(f"Activating forecast mode with reference time {reference_time}")

    def deactivate_forecast_mode(self) -> None:
        self._reference_time = None
        msg.info(f"Deactivating forecast mode")

    def forecast_mode(self) -> None:
        return self._reference_time is not None

    def empty_copy(
        self,
        grid: Grid = None,
        start_time: str = None,
        end_time: str = None,
        source: DataSource = None,
    ):
        new_model = type(self)(
            grid=grid or self.grid(),
            start_time=start_time or self.start_time(),
            end_time=end_time or self.end_time(),
        )
        new_model._source = source or self._source
        new_model._dry_run = self._dry_run
        new_model._global_dry_run = self._global_dry_run
        return new_model
