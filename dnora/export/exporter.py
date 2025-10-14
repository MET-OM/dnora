from dnora import msg
from dnora.file_module import FileNames
from dnora.type_manager.spectral_conventions import (
    SpectralConvention,
    spectral_convention_from_string,
)
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

    def _get_writer(
        self, obj_type: Union[DnoraDataType, DnoraFileType]
    ) -> WriterFunction:
        return self._writer_dict.get(obj_type, self._get_default_writer())

    # def _get_spectral_convention(self) -> SpectralConvention:
    #     """Used only if method is not defined, such as for GeneralWritingFunctions that just dump everything to montly netcdf-files."""
    #     return SpectralConvention.OCEAN

    def __init__(self, model, include_nest: bool = True):
        self.model = model
        self._nest = {}
        if not self.model.parent():
            msg.header(self, "Initializing model exporter...")
            msg.plain(f"Exporting data from '{model.grid().name}'")
        if include_nest and model.nest():
            for name, nest in self.model.nest(get_dict=True).items():
                msg.plain(
                    f"Exporting data from '{name}' nested inside '{model.grid().name}'"
                )
                self._nest[name] = self.__class__(nest)

    def export(
        self,
        obj_type: Union[DnoraDataType, str],
        writer: WriterFunction = None,
        spectral_convention: Union[SpectralConvention, str] = None,
        filename: str = None,
        folder: str = None,
        dateformat: str = None,
        dateformat_folder: str = None,
        dry_run=False,
        **kwargs,
    ) -> None:
        obj_type = data_type_from_string(obj_type)
        spectral_convention = spectral_convention_from_string(spectral_convention)
        writer_function = self._setup_export(obj_type, writer, dry_run)

        if not self.dry_run():
            if not self._silent:
                if self.model.get(obj_type) is None:
                    return
                    # if self.model.parent() is None or (
                    #     self.model.parent() is not None
                    #     and self.model.parent().get(obj_type) is not None
                    #     and self.model.parent().get(obj_type).name
                    #     == self.model[obj_type].name
                    # ):
                msg.header(
                    writer_function,
                    f"Writing {obj_type.name} data from {self.model[obj_type].name}",
                )

            # if self.model.get(obj_type) is None:
            #     msg.info(f"No {obj_type.name} data exists. Won't export anything.")
            #     return

            if spectral_convention is None:
                try:  # GeneralWritingFunction might not have this method defined
                    spectral_convention = writer_function.convention()
                except AttributeError:
                    pass
                    # wanted_convention = self._get_spectral_convention()
            if obj_type in [DnoraDataType.SPECTRA]:
                self.model[obj_type].set_convention(spectral_convention)
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
            dateformat_folder,
            writer_function=writer_function,
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
        dateformat_folder: str,
        writer_function: WriterFunction,
        **kwargs,
    ) -> list[str]:
        # Controls generation of file names using the proper defaults etc.

        edge_object = kwargs.get("edge_object")
        file_object = FileNames(
            format=self._get_default_format(),
            obj_type=obj_type,
            model=self.model,
            filename=filename,
            folder=folder,
            dateformat=dateformat,
            dateformat_folder=dateformat_folder,
            edge_object=edge_object,
        )
        if self.dry_run():
            if not self._silent:
                msg.info("Dry run! No files will be written.")
            output_files = [file_object.get_filepath()]
        else:
            # Write the object using the WriterFunction
            file_object.create_folder()
            output_files = writer_function(self.model, file_object, obj_type, **kwargs)
            if not isinstance(output_files, list):
                output_files = [output_files]

        # Store name and location where file was written
        old_files = self.model._data_exported_to.get(obj_type, [])
        self.model._data_exported_to[obj_type] = old_files + output_files

        if not self._silent:
            msg.to_multifile(output_files)

    def dry_run(self) -> bool:
        return self._dry_run or self.model.dry_run()
