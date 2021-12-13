from copy import copy
import re

# Import objects
from ..grd.grd_mod import Grid
from ..bnd.bnd_mod import Boundary
from ..wnd.wnd_mod import Forcing

# Import abstract classes and needed instances of them
from ..bnd.read import BoundaryReader
from ..bnd.write import BoundaryWriter
from ..bnd.pick import PointPicker

from ..wnd.read import ForcingReader
from ..wnd.write import ForcingWriter

from ..grd.write import GridWriter
from ..grd.process import GridProcessor, TrivialFilter
from ..trg.write import TrGridWriter

from ..dnplot import GridPlotter, TopoPlotter
from ..inp import InputFileWriter
from ..run import ModelExecuter

from typing import Union
# Import default values and auxiliry functions
from .. import msg
from ..defaults import dflt_mdl, dflt_inp, dflt_bnd, dflt_frc, dflt_grd, dflt_plt, list_of_placeholders
from ..aux import create_filename_obj, create_filename_time, add_folder_to_filename, check_if_folder, clean_filename, add_extension

class ModelRun:
    def __init__(self, grid: Grid, start_time: str, end_time: str, name='AnonymousModelRun'):
        self._name = copy(name)
        self._grid = copy(grid)
        self.start_time = copy(start_time)
        self.end_time = copy(end_time)

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

    def export_boundary(self, boundary_writer: BoundaryWriter=None, out_format: str=None, filestring: str=None, datestring: str=None, folder: str=None) -> None:
        """Writes the spectra in the Boundary-object to an external source, e.g.
        a file."""

        if self.boundary() is None:
            raise Exception('Import boundary before exporting!')

        self._boundary_writer = boundary_writer or self._get_boundary_writer()
        if self._boundary_writer is None:
            raise Exception('Define a BoundaryWriter!')

        msg.header(self._boundary_writer, f"Writing boundary spectra from {self.boundary().name()}")

        # Get formats for file names
        out_format = out_format or self._boundary_writer._preferred_format()
        extension = self._boundary_writer._preferred_extension()

        filestring = filestring or dflt_bnd['fs'][out_format]
        datestring = datestring or dflt_bnd['ds'][out_format]
        folderstring = folder or dflt_bnd['fldr'][out_format]

        # Create filenames
        filestring = create_filename_obj(filestring=filestring, objects=[self.grid(), self.forcing(), self.boundary()])
        filestring = create_filename_time(filestring=filestring, times=[self.start_time, self.end_time], datestring=datestring)
        filestring = add_extension(filestring, extension)


        # Folder
        folderstring = create_filename_obj(filestring=folderstring, objects=[self.grid(),self.boundary(), self.forcing()])
        folderstring = create_filename_time(filestring=folderstring, times=[self.start_time, self.end_time], datestring=datestring)
        existed = check_if_folder(folder=folderstring, create=True)
        if not existed:
            msg.plain(f"Creating folder {folderstring}")


        output_files, output_folder = self._boundary_writer(self.boundary(), filename=filestring, folder=folderstring)

        self._boundary_exported_as = output_files
        self._boundary_exported_to = output_folder

        return

    def export_forcing(self, forcing_writer: ForcingWriter=None, out_format: str=None, filestring: str=None, datestring: str=None, folder: str=None) -> None:
        """Writes the forcing data in the Forcing-object to an external source,
        e.g. a file."""

        if self.forcing() is None:
            raise Exception('Import forcing before exporting!')

        self._forcing_writer = forcing_writer or self._get_forcing_writer()

        if self._forcing_writer is None:
            raise Exception('Define a ForcingWriter!')

        msg.header(self._forcing_writer, f"Writing wind forcing from {self.forcing().name()}")

        # Get formats for file names
        out_format = out_format or self._forcing_writer._preferred_format()
        extension = self._forcing_writer._preferred_extension()

        filestring = filestring or dflt_frc['fs'][out_format]
        datestring = datestring or dflt_frc['ds'][out_format]
        folderstring = folder or dflt_frc['fldr'][out_format]

        # Create filenames
        filestring = create_filename_obj(filestring=filestring, objects=[self.grid(), self.forcing(), self.boundary()])
        filestring = create_filename_time(filestring=filestring, times=[self.start_time, self.end_time], datestring=datestring)
        filestring = add_extension(filestring, extension)


        # Folder
        folderstring = create_filename_obj(filestring=folderstring, objects=[self.grid(),self.boundary(), self.forcing()])
        folderstring = create_filename_time(filestring=folderstring, times=[self.start_time, self.end_time], datestring=datestring)
        existed = check_if_folder(folder=folderstring, create=True)
        if not existed:
            msg.plain(f"Creating folder {folderstring}")


        output_files, output_folder = self._forcing_writer(self.forcing(), filename=filestring, folder=folderstring)

        self._forcing_exported_as = output_files
        self._forcing_exported_to = output_folder


        return

    def export_grid(self, grid_writer: Union[GridWriter, TrGridWriter]=None, out_format: str=None, filestring: str=None, infofilestring: str=None, folder: str=None) -> None:
        """Writes the grid data in the Grid-object to an external source,
        e.g. a file."""

        # Try to use defaul grid writer if not provided
        self._grid_writer = grid_writer or self._get_grid_writer()
        if self._grid_writer is None:
            raise Exception('Define a GridWriter!')

        msg.header(self._grid_writer, f"Writing grid topography from {self.grid().name()}")

        # Get formats for file names
        out_format = out_format or self._grid_writer._preferred_format()
        extension = self._grid_writer._preferred_extension()

        filestring = filestring or dflt_grd['fs'][out_format]
        infofilestring = infofilestring or dflt_grd['info'][out_format]
        folderstring = folder or dflt_grd['fldr'][out_format]

        # Create filenames
        filestring = create_filename_obj(filestring=filestring, objects=[self.grid()])
        filestring = add_extension(filestring, extension)

        infofilestring = create_filename_obj(filestring=infofilestring, objects=[self.grid()])
        infofilestring = add_extension(infofilestring, extension)

        # Folder
        folderstring = create_filename_obj(filestring=folderstring, objects=[self.grid(),self.boundary(), self.forcing()])
        existed = check_if_folder(folder=folderstring, create=True)
        if not existed:
            msg.plain(f"Creating folder {folderstring}")

        # Write the grid using the GridWriter object
        output_files, output_folder = self._grid_writer(self.grid(), filename=filestring, infofilename=infofilestring, folder=folderstring)

        self._grid_exported_as = output_files
        self._grid_exported_to = output_folder

        return

    def write_input_file(self, input_file_writer: InputFileWriter=None, out_format=None, filestring=None, datestring=None, folder=None,
                        grid_path: str=None, forcing_path: str=None, boundary_path: str=None,
                        start_time: str=None, end_time: str=None) -> None:
        """Writes the grid data in the Grid-object to an external source,
        e.g. a file."""

        self._input_file_writer = input_file_writer or self._get_input_file_writer()
        if self._input_file_writer is None:
            raise Exception('Define an InputFileWriter!')

        out_format = out_format or self._input_file_writer._preferred_format()
        extension = self._input_file_writer._preferred_extension()

        filestring = filestring or dflt_inp['fs'][out_format]
        datestring = datestring or dflt_inp['ds'][out_format]
        folderstring = folder or dflt_inp['fldr'][out_format]


        start_time = start_time or self.start_time
        end_time = end_time or self.end_time

        # Filename and folder for the input file
        filestring = create_filename_obj(filestring=filestring, objects=[self, self.grid(), self.forcing(), self.boundary()])
        filestring = create_filename_time(filestring=filestring, times=[start_time, end_time], datestring=datestring)
        filestring = add_extension(filestring, extension)

        folderstring = create_filename_obj(filestring=folderstring, objects=[self, self.grid(), self.forcing(), self.boundary()])
        folderstring = create_filename_time(filestring=folderstring, times=[start_time, end_time], datestring=datestring)

        existed = check_if_folder(folder=folderstring, create=True)
        if not existed:
            msg.plain(f"Creating folder {folderstring}")


        grid_path = grid_path or add_folder_to_filename(filename=self.grid_exported_as(out_format), folder=self.grid_exported_to(out_format))
        forcing_path = forcing_path or add_folder_to_filename(filename=self.forcing_exported_as(out_format), folder=self.forcing_exported_to(out_format))
        boundary_path = boundary_path or add_folder_to_filename(filename=self.boundary_exported_as(out_format), folder=self.boundary_exported_to(out_format))

        msg.header(self._input_file_writer, "Writing model input file...")

        output_files, output_folder = self._input_file_writer(grid=self.grid(), forcing=self.forcing(), boundary=self.boundary(), start_time=start_time, end_time=end_time,
                        filename=filestring, folder=folderstring,
                        grid_path=grid_path, forcing_path=forcing_path, boundary_path=boundary_path)


        msg.to_file(add_folder_to_filename(output_files, output_folder))


        self._input_file_written_as = output_files
        self._input_file_written_to = output_folder

        return

    def run_model(self, model_executer: ModelExecuter=None, out_format=None, filestring=None, datestring=None, folder=None) -> None:
        """Run the model."""

        if model_executer is None:
            self._model_executer = self._get_model_executer()
        else:
            self._model_executer = model_executer


        file_not_provided = filestring is None and datestring is None and out_format is None
        folder_not_provided = folder is None

        if out_format is None:
            out_format = self._model_executer._preferred_format()

        # Set possible file- and datestrings
        if filestring is None:
            filestring = dflt_inp['fs'][out_format]

        if datestring is None:
            datestring = dflt_inp['ds'][out_format]

        if folder is None:
            folder = dflt_inp['fldr'][out_format]

        # If no filename information was provided by the user, then use info about where input file was written (if this info exists)
        if file_not_provided and hasattr(self, '_input_file_written_as'):
            input_file = self._input_file_written_as
        else: # Otherwise use default values and hope for the best
            input_file = create_filename_obj(filestring=filestring, objects=[self, self.grid(), self.forcing(), self.boundary()])
            input_file = create_filename_time(filestring=input_file, times=[self.start_time, self.end_time], datestring=datestring)

        input_file = clean_filename(input_file, list_of_placeholders)


        # If no folder information was provided by the user, then use info about where the input file was written to (if this info exists)
        if file_not_provided and hasattr(self, '_input_file_written_to'):
            model_folder = self._input_file_written_to
        else: # Otherwise use default values and hope for the best
            model_folder = create_filename_obj(filestring=folder, objects=[self, self.grid(), self.forcing(), self.boundary()])
            model_folder = create_filename_time(filestring=model_folder, times=[self.start_time, self.end_time], datestring=datestring)

        model_folder = clean_filename(model_folder, list_of_placeholders)
        msg.header(self._model_executer, "Running model...")
        msg.plain(f"Using input file: {add_folder_to_filename(input_file, model_folder)}")
        self._model_executer(input_file=input_file, model_folder=model_folder)

        return

    def plot_grid(self, grid_plotter: GridPlotter=None, grid_processor: GridProcessor=None,
                plain: bool=False, save_fig: bool=False, show_fig: bool=True,
                filestring: str=dflt_plt['fs']['Grid']) -> None:
        """Plot the data in the Grid-object, possibly overlaying data from the
        Boundary- and Forcing-objects."""

        if grid_plotter is None:
            self._grid_plotter = self._get_grid_plotter()
        else:
            self._grid_plotter = grid_plotter

        if self._grid_plotter is None:
            raise Exception('Define a GridPlotter!')


        grid_plot = copy(self.grid())
        if grid_processor is not None:
            grid_plot.process_grid(grid_processor)


        filename = create_filename_obj(filestring=filestring, objects=[self, self.grid(), self.forcing(), self.boundary()])

        fig, filename_out = self._grid_plotter(grid_plot, forcing=self.forcing(), boundary=self.boundary(), filename=filename, plain=plain)

        # Cleans out e.g. #T0 or "__" if they were in the filename
        filename_out = clean_filename(filename_out, list_of_placeholders)

        if save_fig:
            fig.savefig(filename_out, dpi=300)
            msg.to_file(filename_out)
        if show_fig:
            fig.show()

        return

    def name(self) -> str:
        return self._name

    def grid(self) -> str:
        return self._grid

    def forcing(self) -> Forcing:
        if hasattr(self, '_forcing'):
            return self._forcing
        else:
            return None

    def boundary(self) -> Boundary:
        if hasattr(self, '_boundary'):
            return self._boundary
        else:
            return None

    def grid_exported_to(self, out_format: str='General') -> str:
        if hasattr(self, '_grid_exported_to'):
            return self._grid_exported_to
        else:
            return self.grid().folder(folderstring=dflt_grd['fldr'][out_format])

    def forcing_exported_to(self, out_format: str='General') -> str:
        if hasattr(self, '_forcing_exported_to'):
            return self._forcing_exported_to
        elif self.forcing() is not None:
            return self.forcing().folder(folderstring=dflt_frc['fldr'][out_format], datestring=dflt_frc['ds'][out_format])
        else:
            return ''

    def boundary_exported_to(self, out_format: str='General') -> str:
        if hasattr(self, '_boundary_exported_to'):
            return self._boundary_exported_to
        elif self.boundary() is not None:
            return self.boundary().folder(folderstring=dflt_bnd['fldr'][out_format], datestring=dflt_bnd['ds'][out_format])
        else:
            return ''

    def grid_exported_as(self, out_format: str='General') -> str:
        if hasattr(self, '_grid_exported_tas'):
            return self._grid_exported_as
        else:
            return self.grid().filename(filestring=dflt_grd['fs'][out_format])

    def forcing_exported_as(self, out_format: str='General') -> str:
        if hasattr(self, '_forcing_exported_as'):
            return self._forcing_exported_as
        elif self.forcing() is not None:
            return self.forcing().filename(filestring=dflt_frc['fs'][out_format], datestring=dflt_frc['ds'][out_format])
        else:
            return ''

    def boundary_exported_as(self, out_format: str='General') -> str:
        if hasattr(self, '_boundary_exported_as'):
            return self._boundary_exported_as
        elif self.boundary() is not None:
            return self.boundary().filename(filestring=dflt_bnd['fs'][out_format], datestring=dflt_bnd['ds'][out_format])
        else:
            return ''

    def _get_boundary_reader(self) -> BoundaryReader:
        return None

    def _get_boundary_writer(self) -> BoundaryWriter:
        return None

    def _get_forcing_reader(self) -> ForcingReader:
        return None

    def _get_forcing_writer(self) -> ForcingWriter:
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
