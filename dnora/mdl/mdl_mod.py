from __future__ import annotations  # For TYPE_CHECKING
from copy import copy
import pandas as pd
import numpy as np
from geo_skeletons import PointSkeleton
import re
from typing import Union

# Import objects
from ..grd.grd_mod import Grid, TriGrid

from ..file_module import FileNames

# Import abstract classes and needed instances of them
from ..wnd import Forcing
from ..wnd.read import ForcingReader
from .. import wnd

from ..bnd import Boundary
from ..bnd.read import BoundaryReader
from ..pick.point_pickers import PointPicker
from .. import bnd
from .. import pick


from ..spc import Spectra
from ..spc.read import SpectraReader
from .. import spc

from ..wsr import WaveSeries
from ..wsr.read import WaveSeriesReader, SpectraToWaveSeries
from .. import wsr

from ..wlv import WaterLevel
from ..wlv.read import WaterLevelReader
from .. import wlv

from ..ocr import OceanCurrent
from ..ocr.read import OceanCurrentReader
from .. import ocr

from ..ice import IceForcing
from ..ice.read import IceForcingReader
from .. import ice

from ..run.model_executers import ModelExecuter
from ..spectral_grid import SpectralGrid

from geo_skeletons.decorators import add_datavar

from .. import msg
from ..cacher.cache_decorator import cached_reader
from ..converters import convert_swash_mat_to_netcdf
from pathlib import Path

ReaderFunction = Union[
    ForcingReader,
    BoundaryReader,
    SpectraReader,
    WaveSeriesReader,
    WaterLevelReader,
    OceanCurrentReader,
    IceForcingReader,
]
DnoraObject = Union[
    Grid, TriGrid, Forcing, Boundary, WaveSeries, WaterLevel, OceanCurrent, IceForcing
]

OBJECT_STRINGS = [
    "ModelRun",
    "Grid",
    "Forcing",
    "Boundary",
    "Spectra",
    "WaveSeries",
    "WaterLevel",
    "OceanCurrent",
    "IceForcing",
    "SpectralGrid",
]


class ModelRun:
    def __init__(
        self,
        grid: Grid = None,
        start_time: str = "1970-01-01T00:00",
        end_time: str = "2030-12-31T23:59",
        name: str = "AnonymousModelRun",
        dry_run: bool = False,
    ):
        self.name = copy(name)
        self._grid = copy(grid)
        self._time = pd.date_range(start_time, end_time, freq="H")
        self._exported_to = {}
        self._global_dry_run = dry_run
        self._dry_run = False  # Set by methods
        self._reader_dict = {}
        self._point_picker = None
        self._consistency_check(
            objects_to_ignore_get=["ModelRun"],
            objects_to_ignore_import=["ModelRun", "SpectralGrid", "Grid"],
        )

    def _consistency_check(
        self, objects_to_ignore_get: list[str], objects_to_ignore_import: list[str]
    ):
        """Checks that the class contains the proper import and getter methods. This is a safety feature in case more object types are added."""
        for obj_type in OBJECT_STRINGS:
            if obj_type not in objects_to_ignore_get:
                if not hasattr(self, obj_type.lower()) and not hasattr(
                    self, camel_to_snake(obj_type)
                ):
                    raise SyntaxError(
                        f"No getter method self.{obj_type.lower()}() or self.{camel_to_snake(obj_type)}() defined for object {obj_type}!"
                    )

            if obj_type not in objects_to_ignore_import:
                if not hasattr(self, f"import_{obj_type.lower()}"):
                    raise SyntaxError(
                        f"No import method self.import_{obj_type.lower()}() defined for object {obj_type}!"
                    )

    def _setup_import(
        self,
        obj_type: str,
        name: str,
        dry_run: bool,
        reader: ReaderFunction,
    ) -> tuple[ReaderFunction, str]:
        """Sets up readers, names and dry runs porperties for import of object."""

        if obj_type not in OBJECT_STRINGS:
            raise ValueError(f"Dnora object {obj_type} not listed in OBJECT_STRINGS!")

        self._dry_run = dry_run

        reader = reader or self._get_reader(obj_type)
        if reader is None:
            raise Exception(f"Define a {obj_type}Reader!")

        name = name or reader.name()

        if name is None:
            raise ValueError(
                f"Provide either a name or a {obj_type}Reader that will then define the name!"
            )

        msg.header(reader, f"Importing {obj_type}...")
        return reader, name

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
        **kwargs,
    ):
        """Gets the indeces of the point defined by the logical mask with respect to all points available from the reader function."""
        if not self.dry_run():
            lon_all, lat_all, x_all, y_all = reader.get_coordinates(
                grid=self.grid(), start_time=self.start_time(), source=source
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
                **kwargs,
            )
            if len(inds) < 1:
                msg.warning(
                    "PointPicker didn't find any points. Aborting import of boundary."
                )
                return
            return inds

    @cached_reader("Forcing", wnd.read.DnoraNc)
    def import_forcing(
        self,
        forcing_reader: ForcingReader = None,
        name: str = None,
        dry_run: bool = False,
        source: str = "remote",
        **kwargs,
    ) -> None:
        """Imports wind forcing.

        source = 'remote' (default) / '<folder>' / 'met'

        The implementation of this is up to the ForcingReader, and all options might not be functional.
        'met' options will only work in MET Norway internal networks.

        To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
        """

        forcing_reader, name = self._setup_import(
            "Forcing", name, dry_run, forcing_reader
        )

        if not self.dry_run():
            time, u, v, lon, lat, x, y, attributes = forcing_reader(
                grid=self.grid(),
                start_time=self.start_time(),
                end_time=self.end_time(),
                source=source,
                **kwargs,
            )

            self._forcing = Forcing(lon=lon, lat=lat, x=x, y=y, time=time, name=name)
            x = x or lon
            y = y or lat
            self.forcing().set_spacing(nx=len(x), ny=len(y))

            self.forcing().name = name
            self.forcing().set_u(u)
            self.forcing().set_v(v)
            self.forcing().set_metadata(attributes)
        else:
            msg.info("Dry run! No forcing will be imported.")

    @cached_reader("Boundary", bnd.read.DnoraNc)
    def import_boundary(
        self,
        boundary_reader: BoundaryReader = None,
        point_picker: PointPicker = None,
        name: str = None,
        dry_run: bool = False,
        source: str = "remote",
        **kwargs,
    ):
        """Imports boundary spectra. Which spectra to choose spatically
        are determined by the point_picker.

        source = 'remote' (default) / '<folder>' / 'met'

        The implementation of this is up to the BoundaryReader, and all options might not be functional.
        'met' options will only work in MET Norway internal networks.

        To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
        """

        boundary_reader, name = self._setup_import(
            "Boundary", name, dry_run, boundary_reader
        )

        point_picker = self._setup_point_picker(point_picker)

        inds = self._pick_points(
            boundary_reader, point_picker, self.grid().boundary_mask(), source, **kwargs
        )

        if not self.dry_run():
            msg.header(boundary_reader, "Loading boundary spectra...")

            time, freq, dirs, spec, lon, lat, x, y, metadata = boundary_reader(
                grid=self.grid(),
                start_time=self.start_time(),
                end_time=self.end_time(),
                inds=inds,
                source=source,
                **kwargs,
            )
            self._boundary = Boundary(
                x=x, y=y, lon=lon, lat=lat, time=time, freq=freq, dirs=dirs, name=name
            )

            self.boundary().set_spec(spec)
            self.boundary().set_metadata(metadata)
            # E.g. are the spectra oceanic convention etc.
            self.boundary()._convention = boundary_reader.convention()

            self.boundary().set_metadata(
                {"spectral_convention": self.boundary().convention().value}, append=True
            )

            if boundary_reader.post_processing() is not None:
                self.boundary().process_boundary(boundary_reader.post_processing())
            return
        else:
            msg.info("Dry run! No boundary spectra will be imported.")

    @cached_reader("Spectra", spc.read.DnoraNc)
    def import_spectra(
        self,
        spectral_reader: SpectraReader = None,
        point_picker: PointPicker = None,
        name: str = None,
        dry_run: bool = False,
        source: str = "remote",
        **kwargs,
    ):
        """Imports spectra. Which spectra to choose spatically
        are determined by the point_picker.

        source = 'remote' (default) / '<folder>' / 'met'

        The implementation of this is up to the SpectralReader, and all options might not be functional.
        'met' options will only work in MET Norway internal networks.

        To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
        """
        spectral_reader, name = self._setup_import(
            "Spectra", name, dry_run, spectral_reader
        )

        point_picker = self._setup_point_picker(point_picker)

        inds = self._pick_points(
            spectral_reader, point_picker, self.grid().boundary_mask(), source, **kwargs
        )

        if not self.dry_run():
            msg.header(spectral_reader, "Loading omnidirectional spectra...")
            time, freq, spec, mdir, spr, lon, lat, x, y, metadata = spectral_reader(
                grid=self.grid(),
                start_time=self.start_time(),
                end_time=self.end_time(),
                inds=inds,
                source=source,
                **kwargs,
            )

            self._spectra = Spectra(
                x=x, y=y, lon=lon, lat=lat, time=time, freq=freq, name=name
            )

            self.spectra().set_spec(spec)
            self.spectra().set_mdir(mdir)
            self.spectra().set_spr(spr)

            self.spectra().set_metadata(metadata)

            # E.g. are the spectra oceanic convention etc.
            self.spectra()._convention = spectral_reader.convention()
            self.spectra().set_metadata(
                {"spectral_convention": self.spectra().convention().value}, append=True
            )
        else:
            msg.info("Dry run! No spectra will be imported.")

    @cached_reader("WaveSeries", wsr.read.DnoraNc)
    def import_waveseries(
        self,
        waveseries_reader: WaveSeriesReader = None,
        point_picker: PointPicker = None,
        name: str = None,
        dry_run: bool = False,
        source: str = "remote",
        **kwargs,
    ):
        waveseries_reader, name = self._setup_import(
            "WaveSeries", name, dry_run, waveseries_reader
        )

        point_picker = self._setup_point_picker(point_picker)

        inds = self._pick_points(
            waveseries_reader,
            point_picker,
            self.grid().boundary_mask(),
            source,
            **kwargs,
        )

        if not self.dry_run():
            msg.header(waveseries_reader, "Loading wave series data...")
            time, data_dict, lon, lat, x, y, metadata = waveseries_reader(
                grid=self.grid(),
                start_time=self.start_time(),
                end_time=self.end_time(),
                inds=inds,
                source=source,
                **kwargs,
            )

            self._waveseries = WaveSeries(x, y, lon, lat, time=time, name=name)

            for wp, data in data_dict.items():
                self._waveseries = add_datavar(wp.name(), append=True)(
                    self.waveseries()
                )  # Creates .hs() etc. methods
                self.waveseries()._update_datavar(wp.name(), data)
                self.waveseries().set_metadata(
                    {
                        "name": wp.name(),
                        "unit": f"{wp.unit()}",
                        "standard_name": wp.standard_name(),
                    },
                    data_array_name=wp.name(),
                )

            self.waveseries().set_metadata(metadata)  # Global attributes
        else:
            msg.info("Dry run! No waveseries will be imported.")

    @cached_reader("WaterLevel", wlv.read.DnoraNc)
    def import_waterlevel(
        self,
        waterlevel_reader: WaterLevelReader = None,
        name: str = None,
        dry_run: bool = False,
        source: str = "remote",
        **kwargs,
    ) -> None:
        """Imports waterlevel.

        source = 'remote' (default) / '<folder>' / 'met'

        The implementation of this is up to the WaterLevelReader, and all options might not be functional.
        'met' options will only work in MET Norway internal networks.

        To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
        """
        waterlevel_reader, name = self._setup_import(
            "WaterLevel", name, dry_run, waterlevel_reader
        )

        if not self.dry_run():
            time, waterlevel, lon, lat, x, y, attributes = waterlevel_reader(
                grid=self.grid(),
                start_time=self.start_time(),
                end_time=self.end_time(),
                source=source,
                **kwargs,
            )
            self._waterlevel = WaterLevel(
                lon=lon, lat=lat, x=x, y=y, time=time, name=name
            )
            self.waterlevel().set_spacing(nx=len(x or lon), ny=len(y or lat))

            self.waterlevel().name = name
            self.waterlevel().set_waterlevel(waterlevel)
            self.waterlevel().set_metadata(attributes)
        else:
            msg.info("Dry run! No water level data will be imported.")

    @cached_reader("OceanCurrent", ocr.read.DnoraNc)
    def import_oceancurrent(
        self,
        oceancurrent_reader: OceanCurrentReader = None,
        name: str = None,
        dry_run: bool = False,
        source: str = "remote",
        **kwargs,
    ) -> None:
        """Imports waterlevel.

        source = 'remote' (default) / '<folder>' / 'met'

        The implementation of this is up to the OceanCurrentReader, and all options might not be functional.
        'met' options will only work in MET Norway internal networks.

        To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
        """
        oceancurrent_reader, name = self._setup_import(
            "OceanCurrent", name, dry_run, oceancurrent_reader
        )

        if not self.dry_run():
            time, u, v, lon, lat, x, y, attributes = oceancurrent_reader(
                grid=self.grid(),
                start_time=self.start_time(),
                end_time=self.end_time(),
                source=source,
                **kwargs,
            )
            self._oceancurrent = OceanCurrent(
                lon=lon, lat=lat, x=x, y=y, time=time, name=name
            )
            x = x or lon
            y = y or lat
            self.oceancurrent().set_spacing(nx=len(x), ny=len(y))

            self.oceancurrent().name = name
            self.oceancurrent().set_u(u)
            self.oceancurrent().set_v(v)
            self.oceancurrent().set_metadata(attributes)
        else:
            msg.info("Dry run! No water level data will be imported.")

    @cached_reader("IceForcing", ice.read.DnoraNc)
    def import_iceforcing(
        self,
        iceforcing_reader: IceForcingReader = None,
        name: str = None,
        dry_run: bool = False,
        source: str = "remote",
        **kwargs,
    ) -> None:
        """Imports waterlevel.

        source = 'remote' (default) / '<folder>' / 'met'

        The implementation of this is up to the WaterLevelReader, and all options might not be functional.
        'met' options will only work in MET Norway internal networks.

        To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
        """
        iceforcing_reader, name = self._setup_import(
            "IceForcing", name, dry_run, iceforcing_reader
        )
        if not self.dry_run():
            (
                time,
                concentration,
                thickness,
                lon,
                lat,
                x,
                y,
                attributes,
            ) = iceforcing_reader(
                grid=self.grid(),
                start_time=self.start_time(),
                end_time=self.end_time(),
                source=source,
                **kwargs,
            )
            self._iceforcing = IceForcing(
                lon=lon, lat=lat, x=x, y=y, time=time, name=name
            )
            x = x or lon
            y = y or lat
            self.iceforcing().set_spacing(nx=len(x), ny=len(y))

            self.iceforcing().name = name
            self.iceforcing().set_concentration(concentration)
            self.iceforcing().set_thickness(thickness)
            self.iceforcing().set_metadata(attributes)
        else:
            msg.info("Dry run! No ice forcing data will be imported.")

    def boundary_to_spectra(
        self, dry_run: bool = False, name: str = None, write_cache=False, **kwargs
    ):
        self._dry_run = dry_run
        if self.boundary() is None:
            msg.warning("No Boundary to convert to Spectra!")

        spectral_reader = spc.read.BoundaryToSpectra(self.boundary())
        msg.header(
            spectral_reader,
            "Converting the boundary spectra to omnidirectional spectra...",
        )
        name = self.boundary().name

        if not self.dry_run():
            self.import_spectra(
                spectral_reader=spectral_reader,
                point_picker=pick.TrivialPicker(),
                name=name,
                write_cache=write_cache,
                **kwargs,
            )
        else:
            msg.info("Dry run! No boundary will not be converted to spectra.")

    def spectra_to_waveseries(
        self,
        dry_run: bool = False,
        write_cache=False,
        freq: tuple = (0, 10_000),
        **kwargs,
    ):
        self._dry_run = dry_run
        if self.spectra() is None:
            msg.warning("No Spectra to convert to WaveSeries!")
            return
        name = self.spectra().name

        if not self.dry_run():
            waveseries_reader = SpectraToWaveSeries(self.spectra(), freq)
            msg.header(
                waveseries_reader, "Converting the spectra to wave series data..."
            )

            self.import_waveseries(
                waveseries_reader=waveseries_reader,
                point_picker=pick.TrivialPicker(),
                name=name,
                write_cache=write_cache,
                **kwargs,
            )
        else:
            msg.info("Dry run! No boundary will not be converted to spectra.")

    def boundary_to_waveseries(
        self,
        dry_run: bool = False,
        write_cache=False,
        freq: tuple = (0, 10_000),
        **kwargs,
    ):
        self.boundary_to_spectra(dry_run=dry_run, write_cache=write_cache, **kwargs)
        self.spectra_to_waveseries(
            dry_run=dry_run, write_cache=write_cache, freq=freq, **kwargs
        )

    def set_spectral_grid_from_boundary(self, **kwargs):
        if self.boundary() is None:
            msg.warning("No Boundary exists. Can't set spectral grid.")
            return
        self.set_spectral_grid(
            freq=self.boundary().freq(), dirs=self.boundary().dirs(), **kwargs
        )

    def set_spectral_grid_from_spectra(self, **kwargs):
        if self.spectra() is None:
            msg.warning("No Spectra exists. Can't set spectral grid.")
            return
        self.set_spectral_grid(freq=self.spectra().freq(), **kwargs)

    def set_spectral_grid(
        self,
        freq: np.ndarray = None,
        dirs: np.ndarray = None,
        freq0: float = 0.04118,
        nfreq: int = 32,
        ndir: int = 36,
        finc: float = 1.1,
        dirshift: float = None,
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

        self._spectral_grid = SpectralGrid(name="spectral_grid", freq=freq, dirs=dirs)

    def run_model(
        self,
        model_executer: ModelExecuter = None,
        input_file: str = None,
        folder: str = None,
        dateformat: str = None,
        input_file_extension: str = None,
        dry_run: bool = False,
        mat_to_nc: bool = False,
    ) -> None:
        """Run the model."""
        self._dry_run = dry_run
        self._model_executer = model_executer or self._get_model_executer()
        if self._model_executer is None:
            raise Exception("Define a ModelExecuter!")

        # We always assume that the model is located in the folder the input
        # file was written to

        # Option 1) Use user provided
        # Option 2) Use knowledge of where has been exported
        # Option 3) Use default values to guess where is has previously been exported
        exported_path = Path(self.exported_to("input_file")[0])
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
            obj_type="model_executer",
            edge_object="Grid",
        )

        msg.header(self._model_executer, "Running model...")
        msg.plain(f"Using input file: {file_object.get_filepath()}")
        if not self.dry_run():
            self._model_executer(
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
        return self._grid

    def forcing(self) -> Forcing:
        """Returns the forcing object if exists."""
        if hasattr(self, "_forcing"):
            return self._forcing
        else:
            return None

    def boundary(self) -> Boundary:
        """Returns the boundary object if exists."""
        if hasattr(self, "_boundary"):
            return self._boundary
        else:
            return None

    def spectra(self) -> Spectra:
        """Returns the spectral object if exists."""
        if hasattr(self, "_spectra"):
            return self._spectra
        else:
            return None

    def waveseries(self) -> WaveSeries:
        """Returns the wave series object if exists."""
        if hasattr(self, "_waveseries"):
            return self._waveseries
        else:
            return None

    def waterlevel(self) -> WaterLevel:
        """Returns the water level object if exists."""
        if hasattr(self, "_waterlevel"):
            return self._waterlevel
        else:
            return None

    def oceancurrent(self) -> OceanCurrent:
        """Returns the ocean current object if exists."""
        if hasattr(self, "_oceancurrent"):
            return self._oceancurrent
        else:
            return None

    def iceforcing(self) -> IceForcing:
        """Returns the ocean current object if exists."""
        if hasattr(self, "_iceforcing"):
            return self._iceforcing
        else:
            return None

    def spectral_grid(self) -> SpectralGrid:
        """Returns the spectral grid object if exists."""
        if hasattr(self, "_spectral_grid"):
            return self._spectral_grid
        else:
            return None

    def input_file(self) -> None:
        """Only defined to have method for all objects"""
        return None

    def list_of_objects(
        self,
    ) -> list[DnoraObject]:
        """[ModelRun, Boundary] etc."""
        return [x for x in OBJECT_STRINGS if self[x] is not None]

    def dict_of_object_names(self) -> dict[str:str]:
        """{'Boundary': 'NORA3'} etc."""
        dict_of_object_names = {}
        for obj_type in OBJECT_STRINGS:
            if self[obj_type] is None:
                dict_of_object_names[obj_type] = None
            else:
                dict_of_object_names[obj_type] = self[obj_type].name
        return dict_of_object_names

    def exported_to(self, obj_str: str) -> str:
        """Returns the path the object (e.g. grid) was exported to.

        If object has not been exported, the default filename is returned as
        a best guess
        """

        if self[obj_str] is None:
            return [""]

        default_name = FileNames(
            model=self, format=self._get_default_format(), obj_type=obj_str
        ).get_filepath()
        return self._exported_to.get(obj_str.lower(), default_name)

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

    def __getitem__(self, dnora_str: str):
        dnora_str = camel_to_snake(dnora_str)

        if dnora_str == "model_run":
            return self

        if dnora_str == "tri_grid":
            if type(self.grid()) == TriGrid:
                return self.grid()
            else:
                return None
        try:
            return eval(f"self.{dnora_str}()")
        except AttributeError:
            return eval(f"self.{dnora_str.replace('_', '')}()")
        except:
            return None

    def _get_reader(self, obj_type: str):
        if obj_type not in OBJECT_STRINGS:
            raise Warning(
                f"Trying to access a reader function for object {obj_type}, which is not defined in the list of OBJECT_STRINGS!"
            )
        return self._reader_dict.get(obj_type)

    def _get_point_picker(self) -> PointPicker:
        return self._point_picker


def camel_to_snake(string: str) -> str:
    return re.sub(r"(?<!^)(?=[A-Z])", "_", string).lower()
