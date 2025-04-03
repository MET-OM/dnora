from .inputfile.inputfile_writers import InputFileWriter
from dnora.type_manager.dnora_types import DnoraFileType, DnoraDataType
from dnora import msg
from dnora.file_module import FileNames
from dnora.type_manager.model_formats import ModelFormat

from .post_processors import PostProcessor
from .model_runners import ModelRunner
from pathlib import Path
from dnora.type_manager.dnora_types import file_type_from_string
from .decorators import add_write_method, add_run_method
from dnora.defaults import read_environment_variable

@add_run_method(DnoraFileType.WIND)
@add_run_method(DnoraFileType.GRID)
@add_write_method(DnoraFileType.INPUT)
@add_write_method(DnoraFileType.GRID)
@add_write_method(DnoraFileType.TRIGRID)
@add_write_method(DnoraFileType.WIND)
@add_write_method(DnoraFileType.SPECTRA)
@add_write_method(DnoraFileType.SPECTRA1D)
@add_write_method(DnoraFileType.WAVESERIES)
@add_write_method(DnoraFileType.WATERLEVEL)
@add_write_method(DnoraFileType.CURRENT)
@add_write_method(DnoraFileType.ICE)
class ModelExecuter:
    _input_file_writers = {}
    _model_runners = {}

    def _get_default_format(self) -> str:
        return ModelFormat.MODELRUN

    def __init__(self, model):
        self.model = model

    def dry_run(self) -> bool:
        return self._dry_run or self.model.dry_run()

    def _write(self,
        file_type: DnoraFileType | str,
        input_file_writer: InputFileWriter = None,
        filename: str = None,
        folder: str = None,
        dateformat: str = None,
        format: ModelFormat = None,
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
        file_type = file_type_from_string(file_type)
        input_file_writer = input_file_writer or self._input_file_writers.get(file_type)

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
            obj_type=file_type,
            edge_object=DnoraDataType.GRID,
        )
        file_object.create_folder()

        if self.dry_run():
            msg.info("Dry run! No files will be written.")
            output_files = [file_object.get_filepath()]
            for file in output_files:
                msg.to_file(file)
        else:
            # Write the grid using the InputFileWriter object
            output_files = input_file_writer(
                self.model,
                file_object,
                exported_files=self.model.exported_files(),
                **kwargs,
            )
            if not isinstance(output_files, list):
                output_files = [output_files]

        msg.to_multifile(output_files)
        self.model._input_file_exported_to[file_type] = output_files

        return

    #def run_grid(self, model_runner: ModelRunner | None = None, model_folder: str=''):
    #    self._run_model(file_type=DnoraFileType.GRID, model_runner=model_runner, model_folder=model_folder)
        
    def _run(
        self,
        file_type: DnoraFileType | str,
        model_runner: ModelRunner | None = None,
        model_folder: str = '',
        input_file: str | None = None,
        folder: str | None = None,
        dateformat: str | None = None,
        post_processors: list[PostProcessor] | None = None,
        dry_run: bool = False,
        **kwargs,
    ) -> None:
        """Run the model."""
        self._dry_run = dry_run
        file_type = file_type_from_string(file_type)
        model_runner = model_runner or self._model_runners.get(file_type)
        if model_runner is None:
            raise Exception("Define a ModelRunner!")

        # Find location of model executable
        # E.g. For writing GRID and a preferred format of WW3 search for DNORA_GRID_WW3_PATH and DNORA_WW3_PATH
        model_folder = model_folder or read_environment_variable(file_type, model_runner.preferred_format())

        # Option 1) Use user provided
        # Option 2) Use knowledge of where has been exported
        # Option 3) Use default values to guess where is has previously been exported
        exported_path = Path(self.model.input_file_exported_to(file_type)[0])
        primary_file = input_file or exported_path.stem
        primary_folder = folder or str(exported_path.parent)
        file_object = FileNames(
            model=self.model,
            filename=primary_file,
            folder=primary_folder,
            dateformat=dateformat,
            obj_type=file_type,
            format=self._get_default_format(),
            edge_object=DnoraDataType.GRID,
        )

        msg.header(model_runner, "Running model...")
        msg.plain(f"Using input file: {file_object.get_filepath()}")
        if not self.dry_run():
            model_runner(
                file_object=file_object,
                model_folder=model_folder, 
                **kwargs,
            )

            post_processors = post_processors or model_runner.post_processors()

            if post_processors:
                self.post_process(post_processors, file_object, **kwargs)
        else:
            msg.info("Dry run! Model will not run.")

    def post_process(
        self, post_processors: list[PostProcessor], file_object: FileNames, **kwargs
    ) -> None:
        """Post processes model run output, e.g. convert to netcdf or move files"""
        for post_processor in post_processors:
            print(post_processor)
            post_processor(self.model, file_object, **kwargs)
