from __future__ import annotations  # For TYPE_CHECKING
from copy import copy
import pandas as pd
import numpy as np
from geo_skeletons import PointSkeleton
import re
from typing import Union

# Import objects
from ..grid.grid import Grid, TriGrid

from ..file_module import FileNames
from ..dnora_object_type import DnoraDataType
from ..data_sources import DataSource

# Import abstract classes and needed instances of them
from ..wind import Wind
from ..wind.read import ForcingReader
from .. import wind

from ..spectra import Spectra
from ..spectra.read import SpectraReader
from ..pick.point_pickers import PointPicker, NearestGridPoint
from .. import spectra
from .. import pick


from ..spectra1d import Spectra1D
from ..spectra1d.read import Spectra1DReader
from .. import spectra1d

from ..waveseries import WaveSeries
from ..waveseries.read import WaveSeriesReader, SpectraToWaveSeries
from .. import waveseries

from ..waterlevel import WaterLevel
from ..waterlevel.read import WaterLevelReader
from .. import waterlevel

from ..current import Current
from ..current.read import CurrentReader
from .. import current

from ..ice import Ice
from ..ice.read import IceReader
from .. import ice

from ..run.model_executers import ModelExecuter
from ..spectral_grid import SpectralGrid

from geo_skeletons.decorators import add_datavar

from .. import msg
from ..cacher.cache_decorator import cached_reader
from ..converters import convert_swash_mat_to_netcdf
from pathlib import Path
from ..defaults.default_reader import data_sources

from ..export.exporters import Cacher


ReaderFunction = Union[
    ForcingReader,
    SpectraReader,
    Spectra1DReader,
    WaveSeriesReader,
    WaterLevelReader,
    CurrentReader,
    IceReader,
]


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
        # self._consistency_check(
        #     objects_to_ignore_get=[
        #         DnoraDataType.ModelRun,
        #         DnoraDataType.TriGrid,
        #         DnoraDataType.InputFile,
        #     ],
        #     objects_to_ignore_import=[
        #         DnoraDataType.ModelRun,
        #         DnoraDataType.Grid,
        #         DnoraDataType.SpectralGrid,
        #         DnoraDataType.TriGrid,
        #         DnoraDataType.InputFile,
        #     ],
        # )

    def _consistency_check(
        self, objects_to_ignore_get: list[str], objects_to_ignore_import: list[str]
    ):
        """Checks that the class contains the proper import and getter methods. This is a safety feature in case more object types are added."""
        for obj_type in DnoraDataType:
            if obj_type not in objects_to_ignore_get:
                if not hasattr(self, obj_type.value):
                    raise SyntaxError(
                        f"No getter method self.{obj_type.value}() defined for object {obj_type.name}!"
                    )

            if obj_type not in objects_to_ignore_import:
                if not hasattr(self, f"import_{obj_type.value}"):
                    raise SyntaxError(
                        f"No import method self.import_{obj_type.value}() defined for object {obj_type.name}!"
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
                    f"source should be 'local' (DataSource.LOCAL), 'internal' (DataSource.INTERNAL), or 'remote' (DataSource.REMOTE), not {source}!"
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
            lon_all, lat_all, x_all, y_all = reader.get_coordinates(
                grid=self.grid(),
                start_time=self.start_time(),
                source=source,
                folder=folder,
            )

            all_points = PointSkeleton(lon=lon_all, lat=lat_all, x=x_all, y=y_all)

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
                    "PointPicker didn't find any points. Aborting import of boundary."
                )
                return
            return inds

    def _import_data(
        self, obj_type: DnoraDataType, reader, name, source, folder, **kwargs
    ):
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
            **kwargs,
        )

        obj = obj_type.value(name=name, **coord_dict)
        for key, value in data_dict.items():
            if obj.get(key) is None:
                obj = add_datavar(key, append=True)(obj)  # Creates .hs() etc. methods
            obj.set(key, value)

            metaparameter = metaparameter_dict.get(key)
            if metaparameter is None:
                if hasattr(obj_type.value, "meta_dict"):
                    metaparameter = obj_type.value.meta_dict.get(key)

            if metaparameter is not None:
                obj.set_metadata(metaparameter.meta_dict(), data_array_name=key)

        obj.set_metadata(meta_dict)

        self[obj_type] = obj

    def cache(self, obj_type: DnoraDataType | str) -> None:
        Cacher(self).export(obj_type)

    @cached_reader(DnoraDataType.WIND, wind.read.DnoraNc)
    def import_wind(
        self,
        reader: ForcingReader | None = None,
        name: str | None = None,
        dry_run: bool = False,
        source: str | DataSource = DataSource.UNDEFINED,
        folder: str = None,
        **kwargs,
    ) -> None:
        reader, name, source, folder = self._setup_import(
            DnoraDataType.WIND, name, dry_run, reader, source, folder
        )
        self._import_data(DnoraDataType.WIND, reader, name, source, folder, **kwargs)

    @cached_reader(DnoraDataType.WATERLEVEL, waterlevel.read.DnoraNc)
    def import_waterlevel(
        self,
        reader: WaterLevelReader | None = None,
        name: str | None = None,
        dry_run: bool = False,
        source: str | DataSource = DataSource.UNDEFINED,
        folder: str = None,
        **kwargs,
    ) -> None:
        reader, name, source, folder = self._setup_import(
            DnoraDataType.WATERLEVEL, name, dry_run, reader, source, folder
        )
        self._import_data(
            DnoraDataType.WATERLEVEL, reader, name, source, folder, **kwargs
        )

    def import_waveseries(
        self,
        reader: WaveSeriesReader | None = None,
        point_picker: PointPicker | None = None,
        name: str | None = None,
        dry_run: bool = False,
        source: str | DataSource = DataSource.UNDEFINED,
        folder: str = None,
        **kwargs,
    ) -> None:
        reader, name, source, folder = self._setup_import(
            DnoraDataType.WATERLEVEL, name, dry_run, reader, source, folder
        )

        point_picker = self._setup_point_picker(point_picker)

        inds = self._pick_points(
            reader,
            point_picker,
            self.grid().sea_mask(),
            source,
            folder,
            **kwargs,
        )

        self._import_data(
            DnoraDataType.WAVESERIES, reader, name, source, folder, inds=inds, **kwargs
        )

    @cached_reader(DnoraDataType.CURRENT, current.read.DnoraNc)
    def import_current(
        self,
        reader: CurrentReader | None = None,
        name: str | None = None,
        dry_run: bool = False,
        source: str | DataSource = DataSource.UNDEFINED,
        folder: str = None,
        **kwargs,
    ) -> None:
        reader, name, source, folder = self._setup_import(
            DnoraDataType.CURRENT, name, dry_run, reader, source, folder
        )
        self._import_data(DnoraDataType.CURRENT, reader, name, source, folder, **kwargs)

    @cached_reader(DnoraDataType.ICE, ice.read.DnoraNc)
    def import_ice(
        self,
        reader: IceReader | None = None,
        name: str | None = None,
        dry_run: bool = False,
        source: str | DataSource = DataSource.UNDEFINED,
        folder: str = None,
        **kwargs,
    ) -> None:
        reader, name, source, folder = self._setup_import(
            DnoraDataType.ICE, name, dry_run, reader, source, folder
        )
        self._import_data(DnoraDataType.ICE, reader, name, source, folder, **kwargs)

    # @cached_reader(DnoraDataType.Forcing, wind.read.DnoraNc)
    # def import_forcing(
    #     self,
    #     forcing_reader: ForcingReader | None = None,
    #     name: str | None = None,
    #     dry_run: bool = False,
    #     source: str | DataSource = DataSource.UNDEFINED,
    #     folder: str = None,
    #     **kwargs,
    # ) -> None:
    #     """Imports wind forcing.

    #     source = 'remote' (default) / '<folder>' / 'met'

    #     The implementation of this is up to the ForcingReader, and all options might not be functional.
    #     'met' options will only work in MET Norway internal networks.

    #     To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
    #     """

    #     forcing_reader, name, source, folder = self._setup_import(
    #         DnoraDataType.Forcing, name, dry_run, forcing_reader, source, folder
    #     )

    #     if self.dry_run():
    #         msg.info("Dry run! No forcing will be imported.")
    #         return

    #     time, u, v, lon, lat, x, y, attributes = forcing_reader(
    #         grid=self.grid(),
    #         start_time=self.start_time(),
    #         end_time=self.end_time(),
    #         source=source,
    #         folder=folder,
    #         **kwargs,
    #     )

    #     self[DnoraDataType.Forcing] = Forcing(
    #         lon=lon, lat=lat, x=x, y=y, time=time, name=name
    #     )
    #     x = x or lon
    #     y = y or lat
    #     self.forcing().set_spacing(nx=len(x), ny=len(y))

    #     self.forcing().name = name
    #     self.forcing().set_u(u)
    #     self.forcing().set_v(v)
    #     self.forcing().set_metadata(attributes)

    # @cached_reader(DnoraDataType.Boundary, spectra.read.DnoraNc)
    # def import_boundary(
    #     self,
    #     boundary_reader: BoundaryReader | None = None,
    #     point_picker: PointPicker | None = None,
    #     name: str | None = None,
    #     dry_run: bool = False,
    #     source: str | DataSource = DataSource.UNDEFINED,
    #     folder: str = None,
    #     **kwargs,
    # ):
    #     """Imports boundary spectra. Which spectra to choose spatically
    #     are determined by the point_picker.

    #     source = 'remote' / 'internal' / 'local'

    #     If 'folder' is set, then source is assumed to be 'local'

    #     The implementation of this is up to the BoundaryReader, and all options might not be functional.

    #     To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
    #     """

    #     boundary_reader, name, source, folder = self._setup_import(
    #         DnoraDataType.Boundary, name, dry_run, boundary_reader, source, folder
    #     )
    #     point_picker = self._setup_point_picker(point_picker)

    #     inds = self._pick_points(
    #         boundary_reader,
    #         point_picker,
    #         self.grid().boundary_mask(),
    #         source,
    #         folder,
    #         **kwargs,
    #     )

    #     if self.dry_run():
    #         msg.info("Dry run! No boundary spectra will be imported.")
    #         return

    #     msg.header(boundary_reader, "Loading boundary spectra...")

    #     time, freq, dirs, spec, lon, lat, x, y, metadata = boundary_reader(
    #         grid=self.grid(),
    #         start_time=self.start_time(),
    #         end_time=self.end_time(),
    #         inds=inds,
    #         source=source,
    #         folder=folder,
    #         **kwargs,
    #     )
    #     self[DnoraDataType.Boundary] = Boundary(
    #         x=x, y=y, lon=lon, lat=lat, time=time, freq=freq, dirs=dirs, name=name
    #     )

    #     self.boundary().set_spec(spec)
    #     self.boundary().set_metadata(metadata)
    #     # E.g. are the spectra oceanic convention etc.
    #     self.boundary()._mark_convention(boundary_reader.convention())

    #     if boundary_reader.post_processing() is not None:
    #         self.boundary().process_boundary(boundary_reader.post_processing())

    # @cached_reader(DnoraDataType.Spectra, spectra1d.read.DnoraNc)
    # def import_spectra(
    #     self,
    #     spectral_reader: SpectraReader | None = None,
    #     point_picker: PointPicker | None = None,
    #     name: str | None = None,
    #     dry_run: bool = False,
    #     source: str | DataSource = DataSource.UNDEFINED,
    #     folder: str = None,
    #     **kwargs,
    # ):
    #     """Imports spectra. Which spectra to choose spatically
    #     are determined by the point_picker.

    #     source = 'remote' (default) / '<folder>' / 'met'

    #     The implementation of this is up to the SpectralReader, and all options might not be functional.
    #     'met' options will only work in MET Norway internal networks.

    #     To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
    #     """
    #     spectral_reader, name, source, folder = self._setup_import(
    #         DnoraDataType.Spectra, name, dry_run, spectral_reader, source, folder
    #     )

    #     point_picker = self._setup_point_picker(point_picker)

    #     inds = self._pick_points(
    #         spectral_reader,
    #         point_picker,
    #         self.grid().boundary_mask(),
    #         source,
    #         folder,
    #         **kwargs,
    #     )

    #     if self.dry_run():
    #         msg.info("Dry run! No spectra will be imported.")
    #         return

    #     msg.header(spectral_reader, "Loading omnidirectional spectra...")
    #     time, freq, spec, mdir, spr, lon, lat, x, y, metadata = spectral_reader(
    #         grid=self.grid(),
    #         start_time=self.start_time(),
    #         end_time=self.end_time(),
    #         inds=inds,
    #         source=source,
    #         folder=folder,
    #         **kwargs,
    #     )

    #     self[DnoraDataType.Spectra] = Spectra(
    #         x=x, y=y, lon=lon, lat=lat, time=time, freq=freq, name=name
    #     )

    #     self.spectra().set_spec(spec)
    #     self.spectra().set_mdir(mdir)
    #     self.spectra().set_spr(spr)

    #     self.spectra().set_metadata(metadata)

    #     # E.g. are the spectra oceanic convention etc.
    #     self.spectra()._mark_convention(spectral_reader.convention())

    # @cached_reader(DnoraDataType.WaveSeries, waveseries.read.DnoraNc)
    # def import_waveseries(
    #     self,
    #     waveseries_reader: WaveSeriesReader | None = None,
    #     point_picker: PointPicker | None = None,
    #     name: str | None = None,
    #     dry_run: bool = False,
    #     source: str | DataSource = DataSource.UNDEFINED,
    #     folder: str = None,
    #     **kwargs,
    # ):
    #     waveseries_reader, name, source, folder = self._setup_import(
    #         DnoraDataType.WaveSeries, name, dry_run, waveseries_reader, source, folder
    #     )

    #     point_picker = self._setup_point_picker(point_picker)

    #     inds = self._pick_points(
    #         waveseries_reader,
    #         point_picker,
    #         self.grid().boundary_mask(),
    #         source,
    #         folder,
    #         **kwargs,
    #     )

    #     if self.dry_run():
    #         msg.info("Dry run! No waveseries will be imported.")
    #         return

    #     msg.header(waveseries_reader, "Loading wave series data...")
    #     time, data_dict, lon, lat, x, y, metadata = waveseries_reader(
    #         grid=self.grid(),
    #         start_time=self.start_time(),
    #         end_time=self.end_time(),
    #         inds=inds,
    #         source=source,
    #         folder=folder,
    #         **kwargs,
    #     )

    #     self[DnoraDataType.WaveSeries] = WaveSeries(
    #         x, y, lon, lat, time=time, name=name
    #     )

    #     for wp, data in data_dict.items():
    #         self._waveseries = add_datavar(wp.name(), append=True)(
    #             self.waveseries()
    #         )  # Creates .hs() etc. methods
    #         self.waveseries()._update_datavar(wp.name(), data)
    #         self.waveseries().set_metadata(
    #             {
    #                 "name": wp.name(),
    #                 "unit": f"{wp.unit()}",
    #                 "standard_name": wp.standard_name(),
    #             },
    #             data_array_name=wp.name(),
    #         )

    #     self.waveseries().set_metadata(metadata)  # Global attributes

    # @cached_reader(DnoraDataType.WaterLevel, waterlevel.read.DnoraNc)
    # def import_waterlevel(
    #     self,
    #     waterlevel_reader: WaterLevelReader | None = None,
    #     name: str | None = None,
    #     dry_run: bool = False,
    #     source: str | DataSource = DataSource.UNDEFINED,
    #     folder: str = None,
    #     **kwargs,
    # ) -> None:
    #     """Imports waterlevel.

    #     source = 'remote' (default) / '<folder>' / 'met'

    #     The implementation of this is up to the WaterLevelReader, and all options might not be functional.
    #     'met' options will only work in MET Norway internal networks.

    #     To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
    #     """
    #     waterlevel_reader, name, source, folder = self._setup_import(
    #         DnoraDataType.WaterLevel, name, dry_run, waterlevel_reader, source, folder
    #     )

    #     if self.dry_run():
    #         msg.info("Dry run! No water level data will be imported.")
    #         return

    #     time, waterlevel, lon, lat, x, y, attributes = waterlevel_reader(
    #         grid=self.grid(),
    #         start_time=self.start_time(),
    #         end_time=self.end_time(),
    #         source=source,
    #         folder=folder,
    #         **kwargs,
    #     )
    #     self[DnoraDataType.WaterLevel] = WaterLevel(
    #         lon=lon, lat=lat, x=x, y=y, time=time, name=name
    #     )
    #     self.waterlevel().set_spacing(nx=len(x or lon), ny=len(y or lat))

    #     self.waterlevel().name = name
    #     self.waterlevel().set_waterlevel(waterlevel)
    #     self.waterlevel().set_metadata(attributes)

    # @cached_reader(DnoraDataType.OceanCurrent, current.read.DnoraNc)
    # def import_oceancurrent(
    #     self,
    #     oceancurrent_reader: OceanCurrentReader | None = None,
    #     name: str | None = None,
    #     dry_run: bool = False,
    #     source: str | DataSource = DataSource.UNDEFINED,
    #     folder: str = None,
    #     **kwargs,
    # ) -> None:
    #     """Imports waterlevel.

    #     source = 'remote' (default) / '<folder>' / 'met'

    #     The implementation of this is up to the OceanCurrentReader, and all options might not be functional.
    #     'met' options will only work in MET Norway internal networks.

    #     To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
    #     """
    #     oceancurrent_reader, name, source, folder = self._setup_import(
    #         DnoraDataType.OceanCurrent,
    #         name,
    #         dry_run,
    #         oceancurrent_reader,
    #         source,
    #         folder,
    #     )

    #     if self.dry_run():
    #         msg.info("Dry run! No ocean current data will be imported.")
    #         return

    #     time, u, v, lon, lat, x, y, attributes = oceancurrent_reader(
    #         grid=self.grid(),
    #         start_time=self.start_time(),
    #         end_time=self.end_time(),
    #         source=source,
    #         folder=folder,
    #         **kwargs,
    #     )
    #     self[DnoraDataType.OceanCurrent] = OceanCurrent(
    #         lon=lon, lat=lat, x=x, y=y, time=time, name=name
    #     )
    #     self.oceancurrent().set_spacing(nx=len(x or lon), ny=len(y or lat))

    #     self.oceancurrent().name = name
    #     self.oceancurrent().set_u(u)
    #     self.oceancurrent().set_v(v)
    #     self.oceancurrent().set_metadata(attributes)

    # @cached_reader(DnoraDataType.IceForcing, ice.read.DnoraNc)
    # def import_iceforcing(
    #     self,
    #     iceforcing_reader: IceForcingReader | None = None,
    #     name: str | None = None,
    #     dry_run: bool = False,
    #     source: str | DataSource = DataSource.UNDEFINED,
    #     folder: str = None,
    #     **kwargs,
    # ) -> None:
    #     """Imports waterlevel.

    #     source = 'remote' (default) / '<folder>' / 'met'

    #     The implementation of this is up to the WaterLevelReader, and all options might not be functional.
    #     'met' options will only work in MET Norway internal networks.

    #     To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
    #     """
    #     iceforcing_reader, name, source, folder = self._setup_import(
    #         DnoraDataType.IceForcing, name, dry_run, iceforcing_reader, source, folder
    #     )
    #     if self.dry_run():
    #         msg.info("Dry run! No ice forcing data will be imported.")
    #         return
    #     (
    #         time,
    #         concentration,
    #         thickness,
    #         lon,
    #         lat,
    #         x,
    #         y,
    #         attributes,
    #     ) = iceforcing_reader(
    #         grid=self.grid(),
    #         start_time=self.start_time(),
    #         end_time=self.end_time(),
    #         source=source,
    #         folder=folder,
    #         **kwargs,
    #     )
    #     self[DnoraDataType.IceForcing] = IceForcing(
    #         lon=lon, lat=lat, x=x, y=y, time=time, name=name
    #     )

    #     self.iceforcing().set_spacing(nx=len(x or lon), ny=len(y or lat))

    #     self.iceforcing().name = name
    #     self.iceforcing().set_concentration(concentration)
    #     self.iceforcing().set_thickness(thickness)
    #     self.iceforcing().set_metadata(attributes)

    # def boundary_to_spectra(
    #     self,
    #     dry_run: bool = False,
    #     name: str | None = None,
    #     write_cache=False,
    #     **kwargs,
    # ):
    #     self._dry_run = dry_run
    #     if self.boundary() is None:
    #         msg.warning("No Boundary to convert to Spectra!")
    #         return

    #     spectral_reader = spectra1d.read.BoundaryToSpectra(self.boundary())
    #     msg.header(
    #         spectral_reader,
    #         "Converting the boundary spectra to omnidirectional spectra...",
    #     )

    #     name = self.boundary().name

    #     if self.dry_run():
    #         msg.info("Dry run! No boundary will not be converted to spectra.")
    #         return

    #     self.import_spectra(
    #         spectral_reader=spectral_reader,
    #         point_picker=pick.TrivialPicker(),
    #         name=name,
    #         write_cache=write_cache,
    #         **kwargs,
    #     )

    # def spectra_to_waveseries(
    #     self,
    #     dry_run: bool = False,
    #     write_cache=False,
    #     freq: tuple = (0, 10_000),
    #     **kwargs,
    # ):
    #     self._dry_run = dry_run
    #     if self.spectra() is None:
    #         msg.warning("No Spectra to convert to WaveSeries!")
    #         return

    #     name = self.spectra().name

    #     if self.dry_run():
    #         msg.info("Dry run! No boundary will not be converted to spectra.")
    #         return

    #     waveseries_reader = SpectraToWaveSeries(self.spectra(), freq)
    #     msg.header(waveseries_reader, "Converting the spectra to wave series data...")

    #     self.import_waveseries(
    #         waveseries_reader=waveseries_reader,
    #         point_picker=pick.TrivialPicker(),
    #         name=name,
    #         write_cache=write_cache,
    #         **kwargs,
    #     )

    # def boundary_to_waveseries(
    #     self,
    #     dry_run: bool = False,
    #     write_cache=False,
    #     freq: tuple = (0, 10_000),
    #     **kwargs,
    # ):
    #     self.boundary_to_spectra(dry_run=dry_run, write_cache=write_cache, **kwargs)
    #     self.spectra_to_waveseries(
    #         dry_run=dry_run, write_cache=write_cache, freq=freq, **kwargs
    #     )

    # def set_spectral_grid_from_boundary(self, **kwargs):
    #     if self.boundary() is None:
    #         msg.warning("No Boundary exists. Can't set spectral grid.")
    #         return
    #     self.set_spectral_grid(
    #         freq=self.boundary().freq(), dirs=self.boundary().dirs(), **kwargs
    #     )

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
            freq = np.array(
                [freq0 * finc**n for n in np.linspace(0, nfreq - 1, nfreq)]
            )
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

    def run_model(
        self,
        model_executer: ModelExecuter | None = None,
        input_file: str | None = None,
        folder: str | None = None,
        dateformat: str | None = None,
        dry_run: bool = False,
        mat_to_nc: bool = False,
    ) -> None:
        """Run the model."""
        self._dry_run = dry_run
        model_executer = model_executer or self._get_model_executer()
        if model_executer is None:
            raise Exception("Define a ModelExecuter!")

        # We always assume that the model is located in the folder the input
        # file was written to

        # Option 1) Use user provided
        # Option 2) Use knowledge of where has been exported
        # Option 3) Use default values to guess where is has previously been exported
        exported_path = Path(self.exported_to(DnoraDataType.InputFile)[0])
        primary_file = input_file or exported_path.name
        primary_folder = folder  # or str(exported_path.parent)

        # if hasattr(self, "_input_file_writer"):
        #     extension = input_file_extension or self._input_file_writer._extension()
        # else:
        #     extension = input_file_extension or "swn"

        file_object = FileNames(
            model=self,
            filename=primary_file,
            folder=primary_folder,
            dateformat=dateformat,
            obj_type=DnoraDataType.InputFile,
            edge_object=DnoraDataType.Grid,
        )

        msg.header(model_executer, "Running model...")
        msg.plain(f"Using input file: {file_object.get_filepath()}")
        if not self.dry_run():
            model_executer(
                input_file=file_object.filename(), model_folder=file_object.folder()
            )
        else:
            msg.info("Dry run! Model will not run.")
        if mat_to_nc:
            input_file = f"{file_object.folder()}/{self.grid().name}.mat"
            output_file = f"{file_object.folder()}/{self.grid().name}.nc"
            convert_swash_mat_to_netcdf(
                input_file=input_file,
                output_file=output_file,
                lon=self.grid().lon_edges(),
                lat=self.grid().lat_edges(),
                dt=1,
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

    def input_file(self) -> None:
        """Only defined to have method for all objects"""
        return None

    def list_of_objects(
        self,
    ) -> list[DnoraObject]:
        """[ModelRun, Boundary] etc."""
        return [self[x] for x in DnoraDataType if self[x] is not None]

    def dict_of_object_names(self) -> dict[str:str]:
        """{'Boundary': 'NORA3'} etc."""
        dict_of_object_names = {}
        for obj_type in DnoraDataType:
            if self[obj_type] is None:
                dict_of_object_names[obj_type] = None
            else:
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

    def time(self, crop_with: list[str] = None):
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
                dnora_obj = self[obj_str]
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

    def __getitem__(self, dnora_obj: DnoraDataType):
        # dnora_str = camel_to_snake(dnora_str)

        if dnora_obj is None:
            return self
        if dnora_obj == DnoraDataType.TRIGRID:
            if type(self.grid()) == TriGrid:
                return self.grid()
            else:
                return None

        return self._dnora_objects.get(dnora_obj)

    def __setitem__(self, key: DnoraDataType, value: DnoraObject):
        self._dnora_objects[key] = value

    def _get_reader(self, obj_type: DnoraDataType):
        return self._reader_dict.get(obj_type)

    def _get_point_picker(self) -> PointPicker:
        return self._point_picker


# def camel_to_snake(string: str) -> str:
#     return re.sub(r"(?<!^)(?=[A-Z])", "_", string).lower()
