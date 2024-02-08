from .inputfile.inputfile_writers import InputFileWriter
from dnora.dnora_types import DnoraFileType, DnoraDataType
from dnora import msg
from dnora.file_module import FileNames
from dnora.model_formats import ModelFormat

from .post_processors import PostProcessor
from .model_runners import ModelRunner
from pathlib import Path


class ModelExecuter:
    def _get_default_format(self) -> str:
        return ModelFormat.MODELRUN

    def _get_input_file_writer(self):
        return None

    def _get_model_runner(self):
        return None

    def __init__(self, model):
        self.model = model
        self.exported_to = {}

    def dry_run(self) -> bool:
        return self._dry_run or self.model.dry_run()

    def write_input_file(
        self,
        input_file_writer: InputFileWriter = None,
        filename: str = None,
        folder: str = None,
        dateformat: str = None,
        format: ModelFormat = ModelFormat.MODELRUN,
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
            obj_type=DnoraFileType.INPUT,
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
                self.model, file_object, self.exported_to, **kwargs
            )
            if type(output_files) is not list:
                output_files = [output_files]

        self.exported_to[DnoraFileType.INPUT] = output_files

        return

    def run_model(
        self,
        model_runner: ModelRunner | None = None,
        input_file: str | None = None,
        folder: str | None = None,
        dateformat: str | None = None,
        post_processors: list[PostProcessor] | None = None,
        dry_run: bool = False,
        **kwargs,
    ) -> None:
        """Run the model."""
        self._dry_run = dry_run
        model_runner = model_runner or self._get_model_runner()
        if model_runner is None:
            raise Exception("Define a ModelRunner!")

        # We always assume that the model is located in the folder the input
        # file was written to

        # Option 1) Use user provided
        # Option 2) Use knowledge of where has been exported
        # Option 3) Use default values to guess where is has previously been exported
        exported_path = Path(self.model.exported_to(DnoraFileType.INPUT)[0])
        primary_file = input_file or exported_path.name
        primary_folder = folder  # or str(exported_path.parent)

        file_object = FileNames(
            model=self.model,
            filename=primary_file,
            folder=primary_folder,
            dateformat=dateformat,
            obj_type=DnoraFileType.INPUT,
            edge_object=DnoraDataType.GRID,
        )

        msg.header(model_runner, "Running model...")
        msg.plain(f"Using input file: {file_object.get_filepath()}")
        if not self.dry_run():
            model_runner(
                file_object=file_object,
                model_folder=folder,
                **kwargs,
            )

            post_processors = post_processors or model_runner.post_processors()

            if post_processors:
                for post_processor in post_processors:
                    print(post_processor)
                    post_processor(self.model, file_object, **kwargs)
        else:
            msg.info("Dry run! Model will not run.")
