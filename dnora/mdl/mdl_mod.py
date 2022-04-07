from __future__ import annotations # For TYPE_CHECKING
from copy import copy
import re
from pathlib import Path
import yaml
# Import objects
from ..grd.grd_mod import Grid
from ..bnd.bnd_mod import Boundary
from ..wnd.wnd_mod import Forcing
from ..spc.spc_mod import Spectra

# Import abstract classes and needed instances of them
from ..bnd.read import BoundaryReader
from ..bnd.write import BoundaryWriter
from ..bnd.pick import PointPicker

from ..wnd.read import ForcingReader
from ..wnd.write import ForcingWriter

from ..spc.read import SpectralReader, BoundaryToSpectra
from ..spc.write import SpectralWriter

from ..grd.write import GridWriter
from ..grd.process import GridProcessor, TrivialFilter
from ..trg.write import TrGridWriter

from ..dnplot import GridPlotter, TopoPlotter
from ..inp import InputFileWriter
from ..run import ModelExecuter

from ..file_module import FileNames
from typing import Union
# Import default values and auxiliry functions
from .. import msg
from ..bnd.process import processor_for_convention_change

WritingFunction = Union[GridWriter, BoundaryWriter, SpectralWriter, ForcingWriter]

class ModelRun:
    def __init__(self, grid: Grid, start_time: str='1970-01-01T00:00', end_time: str='2030-12-31T23:59', name: str='AnonymousModelRun'):
        self._name = copy(name)
        self._grid = copy(grid)
        self.start_time = copy(start_time)
        self.end_time = copy(end_time)
        self._exported_to = {}

    def import_boundary(self, boundary_reader: BoundaryReader=None, point_picker: PointPicker=None, name: str=None) -> None:
        """Creates a Boundary-object and imports boundary spectra."""

        self._boundary_reader = boundary_reader or self._get_boundary_reader()
        self._point_picker = point_picker or self._get_point_picker()

        if self._boundary_reader is None :
            raise Exception('Define a BoundaryReader!')
        elif self._point_picker is None:
            raise Exception('Define a PointPicker!')

        # Create boundary object
        name = name or type(self._boundary_reader).__name__
        self._boundary = Boundary(grid=self.grid(), name=name)

        # Import the boundary spectra into the Boundary-object
        self.boundary().import_boundary(start_time=self.start_time, end_time=self.end_time, boundary_reader=self._boundary_reader, point_picker=self._point_picker)

        return

    def import_forcing(self, forcing_reader: ForcingReader=None, name: str=None) -> None:
        """Creates a Forcing-objects and imports forcing data."""

        self._forcing_reader = forcing_reader or self._get_forcing_reader()

        if self._forcing_reader is None:
            raise Exception('Define a ForcingReader!')

        # Create forcing object
        name = name or type(self._forcing_reader).__name__
        self._forcing = Forcing(grid=self.grid(), name=name)

        # Import the forcing data into the Forcing-object
        self.forcing().import_forcing(start_time=self.start_time, end_time=self.end_time, forcing_reader=self._forcing_reader)

        return

    def import_spectra(self, spectral_reader: SpectralReader=None, name: str=None) -> None:
        """Creates a Spectra-object and import omnidirectional spectral data."""
        self._spectral_reader = spectral_reader or self._get_spectral_reader()

        if self._spectral_reader is None:
            raise Exception('Define a SpectralReader!')

        # Create forcing object
        name = name or type(self._spectral_reader).__name__
        self._spectra = Spectra(grid=self.grid(), name=name)

        # Import the forcing data into the Forcing-object
        self.spectra().import_spectra(start_time=self.start_time, end_time=self.end_time, spectral_reader=self._spectral_reader)

    def boundary_to_spectra(self):
        if self.boundary() is None:
            msg.warning('No Boundary to convert to Spectra!')

        spectral_reader = BoundaryToSpectra(self.boundary())
        msg.header(spectral_reader, 'Converting the boundary spectra to omnidirectional spectra...')
        name = self.boundary().name()
        self.import_spectra(spectral_reader, name)


    def export_grid(self, grid_writer: Union[GridWriter, TrGridWriter]=None,
                    filename: str=None, folder: str=None, dateformat: str=None,
                    dry_run=False) -> None:
        """Writes the grid data in the Grid-object to an external source,
        e.g. a file."""

        if len(self.grid().topo())==0:
            msg.warning('Grid not meshed so nothing to export!')
            return

        # Try to use defaul grid writer if not provided
        self._grid_writer = grid_writer or self._get_grid_writer()
        if self._grid_writer is None:
            raise Exception('Define a GridWriter!')

        msg.header(self._grid_writer, f"Writing grid topography from {self.grid().name()}")

        output_files = self._export_object(filename, folder, dateformat,
                            writer_function=self._grid_writer,
                            dnora_obj='grid', dry_run=dry_run)

        # Write status file
        infofilename = str(Path(output_files[0]).with_suffix('')) + '_dnora_info.txt'
        self.grid().write_status(filename=infofilename)

    def export_boundary(self, boundary_writer: BoundaryWriter=None,
                        filename: str=None, folder: str=None,
                        dateformat: str=None, dry_run=False) -> None:
        """Writes the spectra in the Boundary-object to an external source, e.g.
        a file."""

        if self.boundary() is None:
            raise Exception('Import boundary before exporting!')

        self._boundary_writer = boundary_writer or self._get_boundary_writer()
        if self._boundary_writer is None:
            raise Exception('Define a BoundaryWriter!')

        msg.header(self._boundary_writer, f"Writing boundary spectra from {self.boundary().name()}")

        boundary_processor = processor_for_convention_change(
                            current_convention = self.boundary().convention(),
                            wanted_convention = self._boundary_writer._convention())

        if boundary_processor is None:
            msg.info(f"Convention ({self.boundary().convention()}) already equals wanted convention ({self._boundary_writer._convention()}).")
        else:
            self.boundary().process_boundary(boundary_processor)

        __ = self._export_object(filename, folder, dateformat,
                            writer_function=self._boundary_writer,
                            dnora_obj='boundary', dry_run=dry_run)

    def export_spectra(self, spectral_writer: SpectralWriter=None,
                        filename: str=None, folder: str=None,
                        dateformat: str=None, dry_run=False) -> None:
        """Writes the spectra in the Spectra-object to an external source, e.g.
        a file."""

        if self.spectra() is None:
            raise Exception('Import spectra before exporting!')

        self._spectral_writer = spectral_writer or self._get_spectral_writer()
        if self._spectral_writer is None:
            raise Exception('Define a SpectralWriter!')

        msg.header(self._spectral_writer, f"Writing omnidirectional spectra from {self.spectra().name()}")

        # NOT IMPPELMENTED YET
        # spectral_processor = processor_for_convention_change(current_convention = self.spectra().convention(), wanted_convention = self._spectral_writer._convention_in())
        # if spectral_processor is None:
        #     msg.info(f"Convention ({self.spectra().convention()}) already equals wanted convention ({self._spectral_writer._convention_in()}).")
        # else:
        #     self.spectra().process_spectra(spectral_processor)

        # Replace #Spectra etc and add file extension

        __ = self._export_object(filename, folder, dateformat,
                            writer_function=self._spectral_writer,
                            dnora_obj='spectra', dry_run=dry_run)

    def export_forcing(self, forcing_writer: ForcingWriter=None,
                        filename: str=None, folder: str=None,
                         dateformat: str=None, dry_run=False) -> None:
        """Writes the forcing data in the Forcing-object to an external source,
        e.g. a file."""

        if self.forcing() is None:
            raise Exception('Import forcing before exporting!')

        self._forcing_writer = forcing_writer or self._get_forcing_writer()

        if self._forcing_writer is None:
            raise Exception('Define a ForcingWriter!')

        msg.header(self._forcing_writer, f"Writing wind forcing from {self.forcing().name()}")

        __ = self._export_object(filename, folder, dateformat,
                            writer_function=self._forcing_writer,
                            dnora_obj='forcing', dry_run=dry_run)

    def write_input_file(self, input_file_writer: InputFileWriter=None,
                        filename=None, folder=None, dateformat=None,
                        grid_path: str=None, forcing_path: str=None,
                        boundary_path: str=None, start_time: str=None,
                        end_time: str=None, dry_run=True) -> None:
        """Writes the grid data in the Grid-object to an external source,
        e.g. a file."""

        self._input_file_writer = input_file_writer or self._get_input_file_writer()
        if self._input_file_writer is None:
            raise Exception('Define an InputFileWriter!')


        msg.header(self._input_file_writer, "Writing model input file...")

        # Controls generation of file names using the proper defaults etc.
        file_object = FileNames(format=self._get_default_format(),
                                clean_names=self._input_file_writer._clean_filename(),
                                list_of_objects=self.list_of_objects(),
                                _filename=filename,
                                _folder=folder,
                                _dateformat=dateformat,
                                extension=self._input_file_writer._extension(),
                                module='wnd')

        # These are used to write info to the input file where to find the
        # forcing etc. data that has previously been exported
        grid_path = grid_path or self.exported_to('grid')[-1]
        forcing_path = forcing_path or self.exported_to('forcing')[-1]
        boundary_path = boundary_path or self.exported_to('boundary')[-1]

        start_time = start_time or self.start_time
        end_time = end_time or self.end_time


        if dry_run:
            msg.info('Dry run! No files will be written.')
            output_files = [file_object.filepath()]
        else:
            # Write the grid using the InputFileWriter object
            output_files = self._input_file_writer(grid=self.grid(),
                            forcing=self.forcing(), boundary=self.boundary(),
                            start_time=start_time, end_time=end_time,
                            filename=file_object.filepaht(),
                            grid_path=grid_path, forcing_path=forcing_path,
                            boundary_path=boundary_path)
            if type(output_files) is not list:
                output_files = [output_files]


        # Store name and location where file was written
        self._exported_to['input_file'] = []
        for file in output_files:
            self._exported_to['input_file'].append(file)

        if self._input_file_writer._im_silent() or dry_run:
            for file in output_files:
                msg.to_file(file)

        return

    def run_model(self, model_executer: ModelExecuter=None,
                input_file: str=None, folder: str=None,
                dateformat: str=None) -> None:
        """Run the model."""

        self._model_executer = model_executedr or self._get_model_executer()
        if self._model_executer is None:
            raise Exception('Define a ModelExecuter!')

        filename, folder = self._prepare_file_output(module='run',
                                filename=filename, folder=folder,
                                extension=self._model_executer._extension(),
                                dateformat=dateformat)

        if self._model_execture._clean_filename():
            filename = clean_filename(filename, self._defaults['list_of_placeholders'])

        # We always assume that the model is located in the folder the input
        # file was written to
        input_path = Path(self.exported_to('input_file'))
        input_file = input_path.name
        model_folder = input_path.parent

        msg.header(self._model_executer, "Running model...")
        msg.plain(f"Using input file: {str(input_path)}")
        self._model_executer(input_file=input_file, model_folder=model_folder)

        return

    def _export_object(self, filename: str, folder: str, dateformat: str,
                    writer_function: WritingFunction, dnora_obj: str,
                    dry_run: bool) -> list[str]:

        # Controls generation of file names using the proper defaults etc.
        file_object = FileNames(format=self._get_default_format(),
                                clean_names=writer_function._clean_filename(),
                                list_of_objects=self.list_of_objects(),
                                _filename=filename,
                                _folder=folder,
                                _dateformat=dateformat,
                                extension=writer_function._extension(),
                                dnora_obj=dnora_obj)

        if dry_run:
            msg.info('Dry run! No files will be written.')
            output_files = [file_object.filepath()]
        else:
            # Write the object using the WriterFunction
            obj = eval(f'self.{dnora_obj}()')
            output_files = writer_function(obj, file_object.filepath())
            if type(output_files) is not list:
                output_files = [output_files]

        # Store name and location where file was written
        self._exported_to[object] = []
        for file in output_files:
            self._exported_to[object].append(file)

        if writer_function._im_silent() or dry_run:
            for file in output_files:
                msg.to_file(file)

        return output_files

    def plot_grid(self, grid_plotter: GridPlotter=None, grid_processor: GridProcessor=None,
                plain: bool=False, save_fig: bool=False, show_fig: bool=True) -> None:
        """Plot the data in the Grid-object, possibly overlaying data from the
        Boundary- and Forcing-objects."""

        if len(self.grid().topo())==0:
            msg.warning('Grid not meshed and nothing to plot!')
            return

        if grid_plotter is None:
            self._grid_plotter = self._get_grid_plotter()
        else:
            self._grid_plotter = grid_plotter

        if self._grid_plotter is None:
            raise Exception('Define a GridPlotter!')


        grid_plot = copy(self.grid())
        if grid_processor is not None:
            grid_plot.process_grid.grid(grid_processor)


        filename = create_filename_obj(filestring=filestring, objects=[self, self.grid(), self.forcing(), self.boundary()])

        fig, filename_out = self._grid_plotter.grid(grid_plot, forcing=self.forcing(), boundary=self.boundary(), filename=filename, plain=plain)

        # Cleans out e.g. #T0 or "__" if they were in the filename
        if filename_out is not None:
            filename_out = clean_filename(filename_out, list_of_placeholders)

        if fig is not None:
            if save_fig:
                fig.savefig(filename_out, dpi=300)
                msg.to_file(filename_out)
            if show_fig:
                fig.show()

        return

    def plot_topo(self, grid_plotter: GridPlotter=None, grid_processor: GridProcessor=None,
                plain: bool=True, save_fig: bool=False, show_fig: bool=True) -> None:
        """Plot the raw data in the Grid-object, possibly overlaying data from the
        Boundary- and Forcing-objects."""


        if len(self.grid().raw_topo()) == 0:
            msg.warning('No topography imported so nothing to plot!')
            return

        if grid_plotter is None:
            self._grid_plotter = self._get_topo_plotter()
        else:
            self._grid_plotter = grid_plotter

        if self._grid_plotter is None:
            raise Exception('Define a GridPlotter!')


        grid_plot = copy(self.grid())
        if grid_processor is not None:
            grid_plot.process_grid.topo(grid_processor)


        filename = create_filename_obj(filestring=filestring, objects=[self, self.grid(), self.forcing(), self.boundary()])

        fig, filename_out = self._grid_plotter.topo(grid_plot, forcing=self.forcing(), boundary=self.boundary(), filename=filename, plain=plain)

        # Cleans out e.g. #T0 or "__" if they were in the filename
        if filename_out is not None:
            filename_out = clean_filename(filename_out, list_of_placeholders)

        if fig is not None:
            if save_fig:
                fig.savefig(filename_out, dpi=300)
                msg.to_file(filename_out)
            if show_fig:
                fig.show()

        return

    def name(self) -> str:
        return self._name

    def grid(self) -> str:
        """Returns the grid object."""
        return self._grid

    def forcing(self) -> Forcing:
        """Returns the forcing object if exists."""
        if hasattr(self, '_forcing'):
            return self._forcing
        else:
            return None

    def boundary(self) -> Boundary:
        """Returns the boundary object if exists."""
        if hasattr(self, '_boundary'):
            return self._boundary
        else:
            return None

    def spectra(self) -> Spectra:
        """Returns the spectral object if exists."""
        if hasattr(self, '_spectra'):
            return self._spectra
        else:
            return None

    def list_of_objects(self) -> list[ModelRun, Grid, Forcing, Boundary, Spectra]:
        return [self, self.grid(), self.forcing(), self.boundary(), self.spectra()]

    def exported_to(self, object: str) -> str:
        """Returns the path the object (e.g. grid) was exported to.

        If object has not been exported, the default filename is returned as
        a best guess
        """

        if eval(f'self.{object}()') is None:
            return ['']

        if self._exported_to.get(object) is not None:
            return self._exported_to[object]

        return ['']


    def _get_default_format(self) -> str:
        return 'ModelRun'

    def _get_boundary_reader(self) -> BoundaryReader:
        return None

    def _get_boundary_writer(self) -> BoundaryWriter:
        return None

    def _get_forcing_reader(self) -> ForcingReader:
        return None

    def _get_forcing_writer(self) -> ForcingWriter:
        return None

    def _get_spectral_reader(self) -> SpectralReader:
        return None

    def _get_spectral_writer(self) -> SpectralWriter:
        return None

    def _get_grid_writer(self) -> GridWriter:
        return None

    def _get_point_picker(self) -> PointPicker:
        return None

    def _get_input_file_writer(self) -> InputFileWriter:
        return None

    def _get_model_executer(self) -> ModelExecuter:
        return None

    def _get_grid_plotter(self) -> GridPlotter:
        return TopoPlotter()

    def _get_topo_plotter(self) -> GridPlotter:
        return TopoPlotter()

    def __repr__(self):
        lines = [f"<dnora ModelRun object> ({type(self).__name__})", f"  Name: {self.name()}"]

        lines.append(f'  Covering time: {self.start_time} - {self.end_time}')

        lines.append(f"  Contains:")
        if self.grid().structured():
            gridtype = 'structured'
        else:
            gridtype = 'unstructured'

        lines.append(f'\tgrid object ({gridtype}): {self.grid().name()} {self.grid().size()}')
        if self.forcing() is not None:
            lines.append(f'\tforcing object: {self.forcing().name()} {self.forcing().size()}')
        else:
            lines.append(f'\tforcing: use .import_forcing()')
        if self.boundary() is not None:
            lines.append(f'\tboundary object: {self.boundary().name()} {self.boundary().size()}')
        else:
            lines.append(f'\tboundary: use .import_boundary()')


        return "\n".join(lines)
