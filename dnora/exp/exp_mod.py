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
from ..dnora_object_type import DnoraObjectType, object_type_from_string
from ..model_formats import ModelFormat

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


@add_export_method(DnoraObjectType.Grid)
@add_export_method(DnoraObjectType.Forcing)
@add_export_method(DnoraObjectType.Boundary)
@add_export_method(DnoraObjectType.Spectra)
@add_export_method(DnoraObjectType.WaveSeries)
@add_export_method(DnoraObjectType.WaterLevel)
@add_export_method(DnoraObjectType.OceanCurrent)
@add_export_method(DnoraObjectType.IceForcing)
class DataExporter:
    _writer_dict = {DnoraObjectType.InputFile: None}

    def _get_default_writer(self) -> WriterFunction:
        return DumpToNc()

    def _get_default_format(self) -> str:
        return ModelFormat.MODELRUN

    def _get_writer(self, obj_type: DnoraObjectType) -> WriterFunction:
        return self._writer_dict.get(obj_type, self._get_default_writer())

    def __init__(self, model):
        self.model = model
        self.exported_to = {}

    def export(
        self,
        obj_type: DnoraObjectType | str,
        writer: str = None,
        filename: str = None,
        folder: str = None,
        dateformat: str = None,
        format: str = None,
        dry_run=False,
        **kwargs,
    ) -> None:
        obj_type = object_type_from_string(obj_type)
        writer_function = self._setup_export(obj_type, writer, dry_run)

        if not self.dry_run():
            try:  # GeneralWritingFunction might not have this method defined
                wanted_convention = writer_function.convention()
            except AttributeError:
                wanted_convention = self._get_spectral_convention()

            try:
                self.model[obj_type]._set_convention(wanted_convention)
            except AttributeError:  # Can only be done for spectra
                pass

        self._export_object(
            obj_type,
            filename,
            folder,
            dateformat,
            writer_function=writer_function,
            format=format,
            **kwargs,
        )

    def _setup_export(
        self, obj_type: DnoraObjectType, writer_function, dry_run: bool
    ) -> WriterFunction:
        self._dry_run = dry_run

        if self.model[obj_type] is None:
            return None

        writer_function = writer_function or self._get_writer(obj_type)

        if writer_function is None:
            raise Exception(f"Define a {obj_type.name}Writer!")

        return writer_function

    def _export_object(
        self,
        obj_type: DnoraObjectType,
        filename: str,
        folder: str,
        dateformat: str,
        writer_function: WriterFunction,
        format: str,
        **kwargs,
    ) -> list[str]:
        # Controls generation of file names using the proper defaults etc.
        if writer_function is None:
            msg.info(f"No {obj_type.name} data exists. Won't export anything.")
            return

        msg.header(
            writer_function,
            f"Writing {obj_type.name} data from {self.model[obj_type].name}",
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
        self.exported_to[obj_type] = output_files

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
            obj_type=DnoraObjectType.InputFile,
            edge_object=DnoraObjectType.Grid,
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

    # def __getitem__(self, obj_type: DnoraObjectType):
    #     """writer = 'grid_writer, forcing_writer etc."""
    #     return self._get_writer(obj_type)

    def dry_run(self) -> bool:
        return self._dry_run or self.model.dry_run()

    def _get_spectral_convention(self) -> SpectralConvention:
        """Used only if method is not defined, such as for GeneralWritingFunctions that just dump everything to montly netcdf-files."""
        return SpectralConvention.OCEAN
