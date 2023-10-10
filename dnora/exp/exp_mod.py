from ..mdl import ModelRun
from .. import msg
from ..file_module import FileNames
from ..bnd.conventions import SpectralConvention
from typing import Union
from .general.general_writing_functions import GeneralWritingFunction, DnoraNc, DumpToNc
from .grd.grid_writers import GridWriter
from .wnd.forcing_writers import ForcingWriter
from .bnd.boundary_writers import BoundaryWriter
from .spc.spectral_writers import SpectralWriter
from .wsr.waveseries_writers import WaveSeriesWriter
from .wlv.waterlevel_writers import WaterLevelWriter
from .ocr.oceancurrent_writers import OceanCurrentWriter
from .ice.iceforcing_writers import IceForcingWriter
from .inp.input_file_writers import InputFileWriter

from .decorators import add_export_method

WriterFunction = Union[
    GeneralWritingFunction,
    GridWriter,
    ForcingWriter,
    BoundaryWriter,
    SpectralWriter,
    WaveSeriesWriter,
    WaterLevelWriter,
    OceanCurrentWriter,
    IceForcingWriter,
    InputFileWriter,
]


@add_export_method("Grid")
@add_export_method("Forcing")
@add_export_method("Boundary")
@add_export_method("Spectra")
@add_export_method("WaveSeries")
@add_export_method("WaterLevel")
@add_export_method("OceanCurrent")
@add_export_method("IceForcing")
class DataExporter:
    def __init__(self, model: ModelRun):
        self.model = model
        self.exported_to = {}

    def _setup_export(
        self, obj_type: str, writer_function, dry_run: bool
    ) -> WriterFunction:
        self._dry_run = dry_run

        if self.model[obj_type] is None:
            return None

        writer_function = writer_function or self[f"{obj_type}_writer"]

        if writer_function is None:
            raise Exception(f"Define a {obj_type}Writer!")

        return writer_function

    def _export_object(
        self,
        obj_type,
        filename: str,
        folder: str,
        dateformat: str,
        writer_function: WriterFunction,
        format: str,
        **kwargs,
    ) -> list[str]:
        # Controls generation of file names using the proper defaults etc.
        if writer_function is None:
            msg.info(f"No {obj_type} data exists. Won't export anything.")
            return

        msg.header(
            writer_function,
            f"Writing {obj_type} data from {self.model[obj_type].name}",
        )

        format = format or self._get_default_format()
        file_object = FileNames(
            format=format,
            obj_type=obj_type,
            model=self.model,
            filename=filename,
            folder=folder,
            dateformat=dateformat,
        )
        if self.dry_run():
            msg.info("Dry run! No files will be written.")
            output_files = [file_object.get_filepath()]
        else:
            # Write the object using the WriterFunction
            file_object.create_folder()
            output_files = writer_function(self.model, file_object, obj_type, **kwargs)
            if type(output_files) is not list:
                output_files = [output_files]

        # Store name and location where file was written
        self.exported_to[obj_type.lower()] = output_files

        for file in output_files:
            msg.to_file(file)

    def write_input_file(
        self,
        input_file_writer: InputFileWriter = None,
        filename=None,
        folder=None,
        dateformat=None,
        format: str = None,
        grid_path: str = None,
        forcing_path: str = None,
        boundary_path: str = None,
        start_time: str = None,
        end_time: str = None,
        dry_run=False,
        **kwargs,
    ) -> None:
        """Writes the grid data in the Grid-object to an external source,
        e.g. a file."""
        self._dry_run = dry_run
        input_file_writer = input_file_writer or self._get_input_file_writer()

        if input_file_writer is None:
            msg.info("No InputFileWriter defines. Won't do anything.")
            return

        msg.header(input_file_writer, "Writing model input file...")

        # Controls generation of file names using the proper defaults etc.
        format = format or self._get_default_format()

        file_object = FileNames(
            model=self.model,
            filename=filename,
            folder=folder,
            format=format,
            dateformat=dateformat,
            obj_type="input_file",
            edge_object="Grid",
        )
        file_object.create_folder()

        if self.dry_run():
            msg.info("Dry run! No files will be written.")
            output_files = [file_object.get_filepath()]
        else:
            # Write the grid using the InputFileWriter object
            output_files = input_file_writer(
                self.model, file_object, self.exported_to, **kwargs
            )
            if type(output_files) is not list:
                output_files = [output_files]

        if self.dry_run():
            for file in output_files:
                msg.to_file(file)

        return

    def __getitem__(self, writer: str):
        """writer = 'grid_writer, forcing_writer etc."""
        return eval(f"self._get_{writer.lower()}()")

    def dry_run(self) -> bool:
        return self._dry_run or self.model.dry_run()

    def _get_default_format(self) -> str:
        return "ModelRun"

    def _get_boundary_writer(self) -> BoundaryWriter:
        return DnoraNc()

    def _get_forcing_writer(self) -> ForcingWriter:
        return DnoraNc()

    def _get_spectra_writer(self) -> SpectralWriter:
        return DnoraNc()

    def _get_waveseries_writer(self) -> WaveSeriesWriter:
        return DnoraNc()

    def _get_waterlevel_writer(self) -> WaterLevelWriter:
        return DnoraNc()

    def _get_grid_writer(self) -> GridWriter:
        return DumpToNc()

    def _get_trigrid_writer(self) -> GridWriter:
        return DnoraNc()

    def _get_input_file_writer(self) -> InputFileWriter:
        return None

    def _get_spectral_convention(self) -> SpectralConvention:
        """Used only if method is not defined, such as for GeneralWritingFunctions that just dump everything to montly netcdf-files."""
        return SpectralConvention.OCEAN
