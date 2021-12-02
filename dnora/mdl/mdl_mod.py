from copy import copy
import re
from ..grd.grd_mod import Grid
from ..bnd.read import BoundaryReader
from ..bnd.write import BoundaryWriter
from ..wnd.read import ForcingReader
from ..wnd.write import ForcingWriter
from ..grd.write import GridWriter

from ..bnd.pick import PointPicker

from ..bnd.bnd_mod import Boundary
from ..wnd.wnd_mod import Forcing

from ..inp import InputFileWriter
from .. import msg
from ..defaults import dflt_mdl, dflt_inp, dflt_bnd, dflt_frc, dflt_grd, dflt_plt, list_of_placeholders

from ..aux import create_filename_obj, create_filename_time, add_folder_to_filename, check_if_folder, clean_filename
from ..grd.process import GridProcessor, TrivialFilter
from ..dnplot import GridPlotter, grd_topo
from ..run import ModelExecuter # Abstract class

class ModelRun:
    def __init__(self, grid: Grid, start_time: str, end_time: str,
                name='AnonymousModelRun'): #, folder=dflt_mdl['fldr']['General']):
        self._name = copy(name)

        self._grid = copy(grid)
        self.start_time = copy(start_time)
        self.end_time = copy(end_time)
        #self._folder = copy(folder)
        #self._foldername_generated = False


    # def set_run_time(self, start_time: str, end_time: str):
    #     self.start_time = start_time
    #     self.end_time = end_time
    #     return
    #
    # def set_grid(self, grid: Grid):
    #     self.grid = copy(grid)
    #     return

    def import_boundary(self, boundary_reader: BoundaryReader=None, point_picker: PointPicker=None, name: str=''):
        if boundary_reader is None:
            self._boundary_reader = self._get_boundary_reader()
        else:
            self._boundary_reader = boundary_reader

        if point_picker is None:
            self._point_picker = self._get_point_picker()
        else:
            self._point_picker = point_picker

        if self._boundary_reader is None :
            raise Exception('Define a BoundaryReader!')
        elif self._point_picker is None:
            raise Exception('Define a PointPicker!')
        else:
            if name:
                self._boundary = Boundary(grid=self.grid(), name=name)
            else:
                self._boundary = Boundary(grid=self.grid(), name=type(self._boundary_reader).__name__)
            self.boundary().import_boundary(start_time=self.start_time, end_time=self.end_time, boundary_reader=self._boundary_reader, point_picker=self._point_picker)

        return None

    def import_forcing(self, forcing_reader: ForcingReader=None, name: str=''):
        if forcing_reader is None:
            self._forcing_reader = self._get_forcing_reader()
        else:
            self._forcing_reader = forcing_reader

        if self._forcing_reader is None:
            raise Exception('Define a ForcingReader!')
        else:
            if name:
                self._forcing = Forcing(grid=self.grid(), name=name)
            else:
                self._forcing = Forcing(grid=self.grid(), name=type(self._forcing_reader).__name__)
            self.forcing().import_forcing(start_time=self.start_time, end_time=self.end_time, forcing_reader=self._forcing_reader)

        return None

    def export_boundary(self, boundary_writer: BoundaryWriter=None, out_format: str=None, filestring: str=None, datestring: str=None, folder: str=None):
        # Make sure we have set the foldername, since the name given
        # might contain regular expressions
        #self.set_foldername()

        if hasattr(self, 'boundary'):
            if boundary_writer is None:
                 self._boundary_writer = self._get_boundary_writer()
            else:
                self._boundary_writer = boundary_writer

            if self._boundary_writer is None:
                raise Exception('Define a BoundaryWriter!')
            else:
                output_files, output_folder = self.boundary().export_boundary(boundary_writer=self._boundary_writer, out_format=out_format, filestring=filestring, datestring=datestring, folder=folder)

            self._boundary_exported_as = output_files
            self._boundary_exported_to = output_folder
        else:
            raise Exception('Import boundary before exporting!')

        return

    def export_forcing(self, forcing_writer: ForcingWriter=None, out_format: str=None, filestring: str=None, datestring: str=None, folder: str=None):
        # Make sure we have set the foldername, since the name given
        # might contain regular expressions
        #self.set_foldername()
        if hasattr(self, 'forcing'):
            if forcing_writer is None:
                 self._forcing_writer = self._get_forcing_writer()
            else:
                self._forcing_writer = forcing_writer

            if self._forcing_writer is None:
                raise Exception('Define a ForcingWriter!')

            output_files, output_folder = self.forcing().export_forcing(forcing_writer=self._forcing_writer, out_format=out_format, filestring=filestring, datestring=datestring, folder=folder)

            self._forcing_exported_as = output_files
            self._forcing_exported_to = output_folder
        else:
            raise Exception('Import forcing before exporting!')
        return

    def export_grid(self, grid_writer: GridWriter=None, out_format: str=None, filestring: str=None, infofilestring: str=None, folder: str=None):
        if grid_writer is None:
             self._grid_writer = self._get_grid_writer()
        else:
            self._grid_writer = grid_writer

        if self._grid_writer is None:
            raise Exception('Define a GridWriter!')

        output_files, output_folder = self.grid().export_grid(grid_writer=self._grid_writer, out_format=out_format, filestring=filestring, infofilestring=infofilestring, folder=folder)

        self._grid_exported_as = output_files
        self._grid_exported_to = output_folder

        return

    def write_input_file(self, input_file_writer: InputFileWriter=None, out_format=None, filestring=None, datestring=None, folder=None,
                        grid_path: str=None, forcing_path: str=None, boundary_path: str=None,
                        start_time: str=None, end_time: str=None):
        if input_file_writer is None:
            self._input_file_writer = self._get_input_file_writer()
        else:
            self._input_file_writer = input_file_writer

        if self._input_file_writer is None:
            raise Exception('Define an InputFileWriter!')

        if out_format is None:
            out_format = self._input_file_writer._preferred_format()

        if filestring is None:
            filestring = dflt_inp['fs'][out_format]

        if datestring is None:
            datestring = dflt_inp['ds'][out_format]

        if folder is None:
            folder = dflt_inp['fldr'][out_format]

        if start_time is None:
            start_time = self.start_time

        if end_time is None:
            end_time = self.end_time

        # Filename and folder for the input file
        filename = create_filename_obj(filestring=filestring, objects=[self, self.grid(), self.forcing(), self.boundary()])
        filename = create_filename_time(filestring=filename, times=[start_time, end_time], datestring=datestring)

        folder = create_filename_obj(filestring=folder, objects=[self, self.grid(), self.forcing(), self.boundary()])
        folder = create_filename_time(filestring=folder, times=[start_time, end_time], datestring=datestring)

        existed = check_if_folder(folder=folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {folder}")


        if grid_path is None:
            grid_path = add_folder_to_filename(filename=self.grid_exported_as(out_format), folder=self.grid_exported_to(out_format))
        if forcing_path is None:
            forcing_path = add_folder_to_filename(filename=self.forcing_exported_as(out_format), folder=self.forcing_exported_to(out_format))
        if boundary_path is None:
            boundary_path = add_folder_to_filename(filename=self.boundary_exported_as(out_format), folder=self.boundary_exported_to(out_format))

        msg.header(self._input_file_writer, "Writing model input file...")

        output_files, output_folder = self._input_file_writer(grid=self.grid(), forcing=self.forcing(), boundary=self.boundary(), start_time=start_time, end_time=end_time,
                        filename=filename, folder=folder,
                        grid_path=grid_path, forcing_path=forcing_path, boundary_path=boundary_path)


        msg.to_file(add_folder_to_filename(output_files, folder))


        self._input_file_written_as = output_files
        self._input_file_written_to = output_folder

        return

    def run_model(self, model_executer: ModelExecuter=None, out_format=None, filestring=None, datestring=None, folder=None):
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

    def plot_grid(self, grid_plotter: GridPlotter=None, grid_processor: GridProcessor=None, plain: bool=False, save_fig: bool=False, show_fig: bool=True, filestring: str=dflt_plt['fs']['Grid']):
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

    # def set_foldername(self, folder: str=dflt_mdl['fldr']['General'], datestring: str=dflt_mdl['ds']['General']):
    #     """Creates a foldername for the object.
    #
    #     The foldername can be based on e.g. the name of the ModelRun object
    #     itself, or the containing Grid or Boundary/Forcing objects (if imported)
    #     ,or the start and end times.
    #
    #     This is typically called by the objects when using the
    #     .export_XXXXXXX() method.
    #     """
    #     if not self._foldername_generated:
    #         # Substitute placeholders for objects ($Grid etc.)
    #         objects=[self, self.grid()]
    #         if hasattr(self, 'forcing'):
    #             objects.append(self.forcing)
    #         if hasattr(self, 'boundary'):
    #             objects.append(self.boundary())
    #         folder = create_filename_obj(filestring=folder, objects=objects)
    #         # Substitute placeholders for times ($T0 etc.)
    #         folder = create_filename_time(filestring=folder, times=[self.start_time, self.end_time], datestring=datestring)
    #
    #         # Possible clean up
    #
    #
    #         self._folder = folder
    #         self._foldername_generated = True
    #     return
    #
    # def folder(self):
    #     return self._folder

    def name(self):
        return self._name

    def grid(self):
        return self._grid

    def forcing(self):
        if hasattr(self, '_forcing'):
            return self._forcing
        else:
            return None

    def boundary(self):
        if hasattr(self, '_boundary'):
            return self._boundary
        else:
            return None

    def grid_exported_to(self, out_format: str='General'):
        if hasattr(self, '_grid_exported_to'):
            return self._grid_exported_to
        else:
            return self.grid().folder(folderstring=dflt_grd['fldr'][out_format])

    def forcing_exported_to(self, out_format: str='General'):
        if hasattr(self, '_forcing_exported_to'):
            return self._forcing_exported_to
        elif self.forcing() is not None:
            return self.forcing().folder(folderstring=dflt_frc['fldr'][out_format], datestring=dflt_frc['ds'][out_format])
        else:
            return ''

    def boundary_exported_to(self, out_format: str='General'):
        if hasattr(self, '_boundary_exported_to'):
            return self._boundary_exported_to
        elif self.boundary() is not None:
            return self.boundary().folder(folderstring=dflt_bnd['fldr'][out_format], datestring=dflt_bnd['ds'][out_format])
        else:
            return ''

    def grid_exported_as(self, out_format: str='General'):
        if hasattr(self, '_grid_exported_tas'):
            return self._grid_exported_as
        else:
            return self.grid().filename(filestring=dflt_grd['fs'][out_format])

    def forcing_exported_as(self, out_format: str='General'):
        if hasattr(self, '_forcing_exported_as'):
            return self._forcing_exported_as
        elif self.forcing() is not None:
            return self.forcing().filename(filestring=dflt_frc['fs'][out_format], datestring=dflt_frc['ds'][out_format])
        else:
            return ''

    def boundary_exported_as(self, out_format: str='General'):
        if hasattr(self, '_boundary_exported_as'):
            return self._boundary_exported_as
        elif self.boundary() is not None:
            return self.boundary().filename(filestring=dflt_bnd['fs'][out_format], datestring=dflt_bnd['ds'][out_format])
        else:
            return ''

    def _get_boundary_reader(self):
        return None

    def _get_boundary_writer(self):
        return None

    def _get_forcing_reader(self):
        return None

    def _get_forcing_writer(self):
        return None

    def _get_grid_writer(self):
        return None

    def _get_point_picker(self):
        return None

    def _get_input_file_writer(self):
        return None

    def _get_model_executer(self):
        return None

    def _get_grid_plotter(self):
        return grd_topo()
