from dnora import msg
from dnora.file_module import FileNames
from dnora.type_manager.spectral_conventions import SpectralConvention
from typing import Union
from .data_writers import DataWriter, Netcdf
from .spectra_writers import SpectraWriter

from .decorators import add_export_method
from dnora.type_manager.dnora_types import (
    DnoraDataType,
    data_type_from_string,
    DnoraFileType,
)
from dnora.type_manager.model_formats import ModelFormat

WriterFunction = Union[
    DataWriter,
    SpectraWriter,
]


@add_export_method(DnoraDataType.GRID)
@add_export_method(DnoraDataType.TRIGRID)
@add_export_method(DnoraDataType.WIND)
@add_export_method(DnoraDataType.SPECTRA)
@add_export_method(DnoraDataType.SPECTRA1D)
@add_export_method(DnoraDataType.WAVESERIES)
@add_export_method(DnoraDataType.WATERLEVEL)
@add_export_method(DnoraDataType.CURRENT)
@add_export_method(DnoraDataType.ICE)
class DataExporter:
    _writer_dict = {}
    _silent = False

    def _get_default_writer(self) -> WriterFunction:
        return Netcdf()

    def _get_default_format(self) -> str:
        return ModelFormat.MODELRUN

    def _get_writer(self, obj_type: DnoraDataType | DnoraFileType) -> WriterFunction:
        return self._writer_dict.get(obj_type, self._get_default_writer())

    def _get_spectral_convention(self) -> SpectralConvention:
        """Used only if method is not defined, such as for GeneralWritingFunctions that just dump everything to montly netcdf-files."""
        return SpectralConvention.OCEAN

    def __init__(self, model):
        self.model = model

    def export(
        self,
        obj_type: DnoraDataType | str,
        writer: WriterFunction = None,
        filename: str = None,
        folder: str = None,
        dateformat: str = None,
        format: str = None,
        dry_run=False,
        **kwargs,
    ) -> None:
        obj_type = data_type_from_string(obj_type)
        writer_function = self._setup_export(obj_type, writer, dry_run)

        if not self.dry_run():
            if self.model.get(obj_type) is None:
                msg.info(f"No {obj_type.name} data exists. Won't export anything.")
                return

            if not self._silent:
                msg.header(
                    writer_function,
                    f"Writing {obj_type.name} data from {self.model[obj_type].name}",
                )

            try:  # GeneralWritingFunction might not have this method defined
                wanted_convention = writer_function.convention()
            except AttributeError:
                wanted_convention = self._get_spectral_convention()

            if obj_type in [DnoraDataType.SPECTRA, DnoraDataType.SPECTRA1D]:
                self.model[obj_type].set_convention(wanted_convention)
        else:
            if not self._silent:
                if self.model.get(obj_type) is None:
                    msg.info(f"Writing {obj_type.name} data")
                else:
                    msg.header(
                        writer_function,
                        f"Writing {obj_type.name} data from {self.model[obj_type].name}",
                    )

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
        self, obj_type: DnoraDataType, writer_function, dry_run: bool
    ) -> WriterFunction:
        self._dry_run = dry_run

        writer_function = writer_function or self._get_writer(obj_type)

        if writer_function is None:
            raise Exception(f"Define a {obj_type.name}Writer!")

        return writer_function

    def _export_object(
        self,
        obj_type: DnoraDataType,
        filename: str,
        folder: str,
        dateformat: str,
        writer_function: WriterFunction,
        format: str,
        **kwargs,
    ) -> list[str]:
        # Controls generation of file names using the proper defaults etc.
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
            if not self._silent:
                msg.info("Dry run! No files will be written.")
            output_files = [file_object.get_filepath()]
        else:
            # Write the object using the WriterFunction
            file_object.create_folder()
            output_files = writer_function(self.model, file_object, obj_type, **kwargs)
            if type(output_files) is not list:
                output_files = [output_files]

        # Store name and location where file was written
        self.model._data_exported_to[data_type_from_string(obj_type.name)] = (
            output_files
        )
        if not self._silent:
            msg.to_multifile(output_files)

    def dry_run(self) -> bool:
        return self._dry_run or self.model.dry_run()
