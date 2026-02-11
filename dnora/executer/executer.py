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
from typing import Union, Optional

from geo_parameters.metaparameter import MetaParameter
import geo_parameters as gp


@add_run_method(DnoraFileType.ATMOSPHERE)
@add_run_method(DnoraFileType.WAVESERIES)
@add_run_method(DnoraFileType.WAVEGRID)
@add_run_method(DnoraFileType.OCEAN)
@add_run_method(DnoraFileType.INPUT)
@add_run_method(DnoraFileType.SPECTRA)
@add_run_method(DnoraFileType.ICE)
@add_run_method(DnoraFileType.WATERLEVEL)
@add_run_method(DnoraFileType.CURRENT)
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
@add_write_method(DnoraFileType.OCEAN)
@add_write_method(DnoraFileType.WAVEGRID)
@add_write_method(DnoraFileType.ATMOSPHERE)
class ModelExecuter:
    _input_file_writers = {}
    _model_runners = {}
    _output_var_aliases = {}
    _output_vars = []

    def _get_default_format(self) -> str:
        return ModelFormat.MODELRUN

    def __init__(self, model, include_nest: bool = True):
        self.model = model
        self.model._input_file_export_format["general"] = self._get_default_format()
        self._model_output_files: dict[DnoraDataType, list[str]] = {}
        self._nest = {}
        if not self.model.parent():
            msg.header(self, "Initializing model executer...")
            msg.plain(f"Will execute '{model.grid().name}'")
        if include_nest and model.nest():
            for name, nest in self.model.nest(get_dict=True).items():
                msg.plain(f"Will execute '{name}' inside '{model.grid().name}'")
                self._nest[name] = self.__class__(nest)
        self.set_output_vars(self._output_vars)

    def dry_run(self) -> bool:
        return self._dry_run or self.model.dry_run()

    def set_output_vars(
        self,
        output_vars: list[Union[MetaParameter, str]],
        verbose: bool = False,
    ) -> None:
        """Set the output variables that will be used when running the model"""
        self._output_vars = []
        self.add_output_vars(output_vars=output_vars, verbose=verbose)

    def add_output_vars(
        self,
        output_vars: list[Union[MetaParameter, str]],
        verbose: bool = False,
    ) -> None:
        """Adds output variables for the wave model while keeping the old ones"""
        if not isinstance(output_vars, list):
            output_vars = [output_vars]

        for var in output_vars:
            if isinstance(var, str):
                if verbose:
                    msg.plain(f"Mapping '{var}' >> '{var}'")
                self._output_vars.append(var)
            elif gp.is_gp(var):
                key = var.find_me_in(self._output_var_aliases.keys(), return_first=True)
                if key is None:
                    msg.info(
                        f"Cannot find variable {var} in the alias list! Skipping..."
                    )
                else:
                    model_var = self._output_var_aliases.get(key)
                    if verbose:
                        msg.plain(f"Mapping {key} >> '{model_var}'")
                    self._output_vars.append(model_var)
            else:
                raise TypeError(
                    f"Variables need to be of type 'str' or geo-parameters, not {var}!"
                )

    def write_pre_process_files(self):
        """Calls the write file methods for all the objects (i.e. write_grid_file(), wrtie_wind_file() etc."""
        for obj_type in DnoraDataType:
            if self.model.get(obj_type) is not None:
                # Use this instead of calling self.export to get export of possible nested grids right
                exec(f"self.write_{obj_type.name.lower()}_file()")

    def pre_process(self):
        """Calls the run methods for all the objects (i.e. run_grid(), run_wind() etc."""
        for obj_type in DnoraDataType:
            if self.model.get(obj_type) is not None:
                # Use this instead of calling self.export to get export of possible nested grids right
                exec(f"self.run_{obj_type.name.lower()}()")

    def _write(
        self,
        file_type: Union[DnoraFileType, str],
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
            msg.info("No InputFileWriter defined. Won't do anything.")
            return

        msg.header(input_file_writer, f"Writing model input file for '{self.model.grid().name}'...")

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
                output_vars=self._output_vars,
                stride = self.model._stride,
                **kwargs,
            )
            if not isinstance(output_files, list):
                output_files = [output_files]
        if any(fn for fn in output_files):
            msg.to_multifile(output_files)
        self.model._input_file_exported_to[file_type] = output_files
        self.model._input_file_export_format[file_type] = self._get_default_format()

        return

    def run_model(
        self,
        model_runner: Optional[ModelRunner] = None,
        model_folder: str = "",
        post_process: bool = True,
        dry_run: bool = False,
        **kwargs,
    ):
        """Run the main model. Set post_process=False to disable any post-processing that might be defined."""

        # Use the method generated by the decorateor, since that will automatically go through all nested grids if present

        
        self.run_input(
            model_runner=model_runner,
            model_folder=model_folder,
            post_process=post_process,
            dry_run=dry_run,
            **kwargs,
        )

    def _run(
        self,
        file_type: Union[DnoraFileType, str],
        model_runner: Optional[ModelRunner] = None,
        model_folder: str = "",
        input_file: Optional[str] = None,
        folder: Optional[str] = None,
        dateformat: Optional[str] = None,
        post_process: bool = True,
        post_processors: Optional[list[PostProcessor]] = None,
        parent_folder: Optional[str] = None,
        dry_run: bool = False,
        **kwargs,
    ) -> None:
        """Run the model."""
       
        if self.model.parent() is not None:
            parent_folder = parent_folder or str(Path(self.model.parent().input_file_exported_to(DnoraFileType.INPUT)[0]).parent)
        else:
            parent_folder = ''        
        self._dry_run = dry_run
        file_type = file_type_from_string(file_type)
        model_runner = model_runner or self._model_runners.get(file_type)
        if model_runner is None:
            raise Exception("Define a ModelRunner!")
        msg.header(model_runner, f"Running model '{self.model.grid().name}'...")
        # Find location of model executable
        # E.g. For writing GRID and a preferred format of WW3 search for DNORA_GRID_WW3_PATH and DNORA_WW3_PATH
        model_folder = model_folder or read_environment_variable(
            file_type, model_runner.preferred_format()
        )

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

        
        msg.plain(f"Using input file: {file_object.get_filepath()}")
        if not self.dry_run():
            outfile = model_runner(
                file_object=file_object,
                model_folder=model_folder,
                parent_folder=parent_folder,
                **kwargs,
            )
            if outfile is not None:
                if not isinstance(outfile,list):
                    outfile = [outfile]
                self._model_output_files[file_type] = outfile
        else:
            msg.info("Dry run! Model will not run.")

        post_processors = post_processors or model_runner.post_processors()

        if post_processors and post_process:
            self.post_process(post_processors, file_object, model_folder, parent_folder,file_type, **kwargs)

    def post_process(
        self,
        post_processors: list[PostProcessor],
        file_object: FileNames,
        model_folder,
        parent_folder: str,
        file_type,
        **kwargs,
    ) -> None:
        """Post processes model run output, e.g. convert to netcdf or move files"""
        for post_processor in post_processors:
            msg.header(post_processor, "Post processing...")
            if post_processor.for_nest is not None:
                model = self.model.nest(get_dict=True)[post_processor.for_nest]
            else:
                model = self.model
            file_objects = []
            for fn in model.input_file_exported_to(post_processor.for_file_type or file_type):
                
                exported_path = Path(fn)
                primary_file = exported_path.stem
                primary_folder = str(exported_path.parent)
                file_objects.append(FileNames(
                    model=model,
                    filename=primary_file,
                    folder=primary_folder,
                    obj_type=post_processor.for_file_type or file_type,
                    format=self._get_default_format(),
                    edge_object=DnoraDataType.GRID,
                ))

            post_processor(
                file_object=file_objects,
                model_folder=model_folder,
                parent_folder=parent_folder,
                **kwargs,
            )

    def _output_file(self, obj_type: Union[DnoraFileType, str]) -> str:
        """Returns the path the object (e.g. grid) was exported to.

        If object has not been exported, the default filename is returned as
        a best guess
        """
        obj_type = file_type_from_string(obj_type)
        export_format = (
            self._data_export_format.get(obj_type)
            or self._data_export_format.get("general")
            or self._get_default_format()
        )
        default_name = [
            FileNames(
                model=self, format=export_format, obj_type=obj_type
            ).get_filepath()
        ]  # Want a list of strings
        return self._data_exported_to.get(obj_type, default_name)

    def output_files(self) -> dict:
        """Gives a dict of the exported files"""
        files = {}
        for dnora_type in DnoraFileType:
            dnora_type = file_type_from_string(dnora_type)
            output_files = self._model_output_files.get(dnora_type)
            if output_files is not None:
                files[dnora_type.name.lower()] = output_files
        return files