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
from ..wlv.wlv_mod import WaterLevel
from ..ocr.ocr_mod import OceanCurrent

# Import abstract classes and needed instances of them
from ..bnd.read import BoundaryReader
from ..bnd.write import BoundaryWriter
from ..bnd.pick import PointPicker

from ..wnd.read import ForcingReader
from ..wnd.write import ForcingWriter

from ..wlv.read import WaterLevelReader
from ..wlv.write import WaterLevelWriter

from ..ocr.read import OceanCurrentReader
from ..ocr.write import OceanCurrentWriter

from ..spc.read import SpectralReader, BoundaryToSpectra
from ..spc.write import SpectralWriter

from ..grd.write import GridWriter
from ..grd.process import GridProcessor, TrivialFilter
from ..trg.write import TrGridWriter

from ..dnplot.dnplot import GridPlotter, TopoPlotter, ForcingPlotter, OceanCurrentPlotter
from ..inp import InputFileWriter
from ..run import ModelExecuter

from ..file_module import FileNames
from typing import Union
# Import default values and aux_funcsiliry functions
from .. import msg
from ..bnd.process import processor_for_convention_change

from .. import file_module
from ..converters import convert_swash_mat_to_netcdf
WritingFunction = Union[GridWriter, BoundaryWriter, SpectralWriter, ForcingWriter, WaterLevelWriter, OceanCurrentWriter]
PlottingFunction = Union[GridPlotter]

class ModelRun:
    def __init__(self, grid: Grid, start_time: str='1970-01-01T00:00',
    end_time: str='2030-12-31T23:59', name: str='AnonymousModelRun',
    dry_run: bool=False):
        self._name = copy(name)
        self._grid = copy(grid)
        self.start_time = copy(start_time)
        self.end_time = copy(end_time)
        self._exported_to = {}
        self._global_dry_run = dry_run
        self._dry_run = False  # Set by methods

    def import_boundary(self, boundary_reader: BoundaryReader=None,
                        point_picker: PointPicker=None, name: str=None,
                        dry_run: bool=False,
                        write_cache: bool=False,
                        read_cache: bool=False,
                        cache_name: str='#Grid_#Lon0_#Lon1_#Lat0_#Lat1') -> None:
        """Creates a Boundary-object and imports boundary spectra."""

        self._dry_run = dry_run
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
        if not self.dry_run():
            self.boundary().import_boundary(start_time=self.start_time,
                                            end_time=self.end_time,
                                            boundary_reader=self._boundary_reader,
                                            point_picker=self._point_picker,
                                            write_cache=write_cache,
                                            read_cache=read_cache,
                                            cache_name=cache_name)
        else:
            msg.info('Dry run! No boundary spectra will be imported.')

    def import_forcing(self, forcing_reader: ForcingReader=None,
                        name: str=None, dry_run: bool=False,
                        expansion_factor: float=1.2,
                        write_cache: bool=False,
                        read_cache: bool=False,
                        cache_name: str='#Grid_#Lon0_#Lon1_#Lat0_#Lat1') -> None:
        """Creates a Forcing-objects and imports forcing data."""
        self._dry_run = dry_run

        self._forcing_reader = forcing_reader or self._get_forcing_reader()

        if self._forcing_reader is None:
            raise Exception('Define a ForcingReader!')

        # Create forcing object
        name = name or type(self._forcing_reader).__name__
        self._forcing = Forcing(grid=self.grid(), name=name)

        # Import the forcing data into the Forcing-object
        if not self.dry_run():
            self.forcing().import_forcing(start_time=self.start_time,
                                        end_time=self.end_time,
                                        forcing_reader=self._forcing_reader,
                                        expansion_factor=expansion_factor,
                                        read_cache=read_cache,
                                        write_cache=write_cache,
                                        cache_name=cache_name)
        else:
            msg.info('Dry run! No forcing will be imported.')


    def import_oceancurrent(self, oceancurrent_reader: OceanCurrentReader=None,
                        name: str=None, dry_run: bool=False,
                        expansion_factor: float=1.2,
                        write_cache: bool=False,
                        read_cache: bool=False,
                        cache_name: str='#Grid_#Lon0_#Lon1_#Lat0_#Lat1') -> None:
        """Creates an OceanCurrent-object and imports oceancurrent data."""
        self._dry_run = dry_run

        self._oceancurrent_reader = oceancurrent_reader or self._get_oceancurrent_reader()

        if self._oceancurrent_reader is None:
            raise Exception('Define an OceanCurrent Reader!')

        # Create oceancurrent object
        name = name or type(self._oceancurrent_reader).__name__
        self._oceancurrent = OceanCurrent(grid=self.grid(), name=name)

        # Import the oceancurrent data into the Ocean Current-object
        if not self.dry_run():
            self.oceancurrent().import_oceancurrent(start_time=self.start_time,
                                        end_time=self.end_time,
                                        oceancurrent_reader=self._oceancurrent_reader,
                                        expansion_factor=expansion_factor,
                                        read_cache=read_cache,
                                        write_cache=write_cache,
                                        cache_name=cache_name)
        else:
            msg.info('Dry run! No OceanCurrent oceancurrent will be imported.')



    def import_waterlevel(self, waterlevel_reader: WaterLevelReader=None,
                        name: str=None, dry_run: bool=False,
                        expansion_factor: float=1.2,
                        write_cache: bool=False,
                        read_cache: bool=False,
                        cache_name: str='#Grid_#Lon0_#Lon1_#Lat0_#Lat1') -> None:
        """Creates a Forcing-objects and imports forcing data."""
        self._dry_run = dry_run

        self._waterlevel_reader = waterlevel_reader or self._get_waterlevel_reader()

        if self._waterlevel_reader is None:
            raise Exception('Define a WaterLevelReader!')

        # Create waterlevel object
        name = name or type(self._waterlevel_reader).__name__
        self._waterlevel = WaterLevel(grid=self.grid(), name=name)

        # Import the forcing data into the Forcing-object
        if not self.dry_run():
            self.waterlevel().import_waterlevel(start_time=self.start_time,
                                        end_time=self.end_time,
                                        waterlevel_reader=self._waterlevel_reader,
                                        expansion_factor=expansion_factor,
                                        read_cache=read_cache,
                                        write_cache=write_cache,
                                        cache_name=cache_name)
        else:
            msg.info('Constant water level! No time dependent water level data specified.')

    def import_spectra(self, spectral_reader: SpectralReader=None,
                        name: str=None, dry_run: bool=False) -> None:
        """Creates a Spectra-object and import omnidirectional spectral data."""
        self._dry_run = dry_run
        self._spectral_reader = spectral_reader or self._get_spectral_reader()

        if self._spectral_reader is None:
            raise Exception('Define a SpectralReader!')

        # Create forcing object
        name = name or type(self._spectral_reader).__name__
        self._spectra = Spectra(grid=self.grid(), name=name)

        # Import the forcing data into the Forcing-object
        if not self.dry_run():
            self.spectra().import_spectra(start_time=self.start_time,
                                        end_time=self.end_time,
                                        spectral_reader=self._spectral_reader)
        else:
            msg.info('Dry run! No omnidirectional spectra will be imported.')


    def boundary_to_spectra(self, dry_run: bool=False):
        self._dry_run = dry_run
        if self.boundary() is None:
            msg.warning('No Boundary to convert to Spectra!')

        spectral_reader = BoundaryToSpectra(self.boundary())
        msg.header(spectral_reader, 'Converting the boundary spectra to omnidirectional spectra...')
        name = self.boundary().name()
        if not self.dry_run():
            self.import_spectra(spectral_reader, name)
        else:
            msg.info('Dry run! No boundary will not be converted to spectra.')

    def export_grid(self, grid_writer: Union[GridWriter, TrGridWriter]=None,
                    filename: str=None, folder: str=None, dateformat: str=None,
                    dry_run=False) -> None:
        """Writes the grid data in the Grid-object to an external source,
        e.g. a file."""
        self._dry_run = dry_run
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
                            dnora_obj='grid')

        # Write status file
        infofilename = str(Path(output_files[0]).with_suffix('')) + '_dnora_info.txt'
        self.grid().write_status(filename=infofilename)

    def export_boundary(self, boundary_writer: BoundaryWriter=None,
                        filename: str=None, folder: str=None,
                        dateformat: str=None, dry_run=False) -> None:
        """Writes the spectra in the Boundary-object to an external source, e.g.
        a file."""
        self._dry_run = dry_run
        if self.boundary() is None and not self.dry_run():
            raise Exception('Import boundary before exporting!')

        self._boundary_writer = boundary_writer or self._get_boundary_writer()
        if self._boundary_writer is None:
            raise Exception('Define a BoundaryWriter!')

        if self.boundary() is None:
            msg.header(self._boundary_writer, f"Writing boundary spectra from DryRunBoundary")
        else:
            msg.header(self._boundary_writer, f"Writing boundary spectra from {self.boundary().name()}")

        if self.dry_run():
            boundary_processor = None
        else:
            boundary_processor = processor_for_convention_change(
                                current_convention = self.boundary().convention(),
                                wanted_convention = self._boundary_writer._convention())

        if boundary_processor is None:
            msg.info(f"Convention ({self.boundary().convention()}) already equals wanted convention ({self._boundary_writer._convention()}).")
        else:
            self.boundary().process_boundary(boundary_processor)

        __ = self._export_object(filename, folder, dateformat,
                            writer_function=self._boundary_writer,
                            dnora_obj='boundary')

    def export_spectra(self, spectral_writer: SpectralWriter=None,
                        filename: str=None, folder: str=None,
                        dateformat: str=None, dry_run=False) -> None:
        """Writes the spectra in the Spectra-object to an external source, e.g.
        a file."""
        self._dry_run = dry_run
        if self.spectra() is None and not self.dry_run():
            raise Exception('Import spectra before exporting!')

        self._spectral_writer = spectral_writer or self._get_spectral_writer()
        if self._spectral_writer is None:
            raise Exception('Define a SpectralWriter!')

        if self.spectra() is None:
            msg.header(self._spectral_writer, f"Writing omnidirectional spectra from DryRunSpectra")
        else:
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
                            dnora_obj='spectra')

    def export_forcing(self, forcing_writer: ForcingWriter=None,
                        filename: str=None, folder: str=None,
                         dateformat: str=None, dry_run=False) -> None:
        """Writes the forcing data in the Forcing-object to an external source,
        e.g. a file."""
        self._dry_run = dry_run
        if self.forcing() is None and not self.dry_run():
            raise Exception('Import forcing before exporting!')

        self._forcing_writer = forcing_writer or self._get_forcing_writer()

        if self._forcing_writer is None:
            raise Exception('Define a ForcingWriter!')

        if self.forcing() is not None:
            msg.header(self._forcing_writer, f"Writing wind forcing from {self.forcing().name()}")
        else:
            msg.header(self._forcing_writer, f"Writing wind forcing from DryRunForcing")


        __ = self._export_object(filename, folder, dateformat,
                            writer_function=self._forcing_writer,
                            dnora_obj='forcing')

    def export_oceancurrent(self, oceancurrent_writer: OceanCurrentWriter=None,
                        filename: str=None, folder: str=None,
                         dateformat: str=None, dry_run=False) -> None:
        """Writes the OceanCurrent data in the OceanCurrent-object to an external source,
        e.g. a file."""
        self._dry_run = dry_run
        if self.oceancurrent() is None and not self.dry_run():
            raise Exception('Import OceanCurrent before exporting!')

        self._oceancurrent_writer = oceancurrent_writer or self._get_oceancurrent_writer()

        if self._oceancurrent_writer is None:
            raise Exception('Define a OceanCurrentWriter!')

        if self.oceancurrent() is not None:
            msg.header(self._oceancurrent_writer, f"Writing oceancurrent from {self.oceancurrent().name()}")
        else:
            msg.header(self._oceancurrent_writer, f"Writing oceancurrent from DryRunOceanCurrent")


        __ = self._export_object(filename, folder, dateformat,
                            writer_function=self._oceancurrent_writer,
                            dnora_obj='oceancurrent')

    def export_waterlevel(self, waterlevel_writer: WaterLevelWriter=None,
                        filename: str=None, folder: str=None,
                         dateformat: str=None, dry_run=False) -> None:
        """Writes the forcing data in the Forcing-object to an external source,
        e.g. a file."""
        self._dry_run = dry_run
        if self.waterlevel() is None and not self.dry_run():
            raise Exception('Import waterlevel before exporting!')

        self._waterlevel_writer = waterlevel_writer or self._get_waterlevel_writer()

        if self._waterlevel_writer is None:
            raise Exception('Define a WaterLevelWriter!')

        if self.waterlevel() is not None:
            msg.header(self._waterlevel_writer, f"Writing water level from {self.waterlevel().name()}")
        else:
            msg.header(self._waterlevel_writer, f"Static water level.")


        __ = self._export_object(filename, folder, dateformat,
                            writer_function=self._waterlevel_writer,
                            dnora_obj='waterlevel')

    def write_input_file(self, input_file_writer: InputFileWriter=None,
                        filename=None, folder=None, dateformat=None,
                        grid_path: str=None, forcing_path: str=None,
                        boundary_path: str=None, waterlevel_path: str=None,
                        oceancurrent_path: str=None,
                        start_time: str=None, end_time: str=None,
                        dry_run=False) -> None:
        """Writes the grid data in the Grid-object to an external source,
        e.g. a file."""
        self._dry_run = dry_run
        self._input_file_writer = input_file_writer or self._get_input_file_writer()
        if self._input_file_writer is None:
            raise Exception('Define an InputFileWriter!')


        msg.header(self._input_file_writer, "Writing model input file...")

        # Controls generation of file names using the proper defaults etc.
        file_object = FileNames(format=self._get_default_format(),
                                clean_names=self._input_file_writer._clean_filename(),
                                dict_of_object_names=self.dict_of_object_names(),
                                start_time=self.start_time,
                                end_time=self.end_time,
                                _filename=filename,
                                _folder=folder,
                                _dateformat=dateformat,
                                extension=self._input_file_writer._extension(),
                                dnora_obj='input_file')

        file_object.create_folder()

        # These are used to write info to the input file where to find the
        # forcing etc. data that has previously been exported
        grid_path = grid_path or self.exported_to('grid')[-1]
        forcing_path = forcing_path or self.exported_to('forcing')[-1]
        boundary_path = boundary_path or self.exported_to('boundary')[-1]
        waterlevel_path = waterlevel_path or self.exported_to('waterlevel')[-1]
        oceancurrent_path = oceancurrent_path or self.exported_to('oceancurrent')[-1]

        start_time = start_time or self.start_time
        end_time = end_time or self.end_time


        if self.dry_run():
            msg.info('Dry run! No files will be written.')
            output_files = [file_object.filepath()]
        else:
            # Write the grid using the InputFileWriter object
            output_files = self._input_file_writer(grid=self.grid(),
                            forcing=self.forcing(), boundary=self.boundary(), waterlevel=self.waterlevel(),  oceancurrent=self.oceancurrent(),
                            start_time=start_time, end_time=end_time,
                            filename=file_object.filepath(),
                            grid_path=grid_path, forcing_path=forcing_path,
                            boundary_path=boundary_path, waterlevel_path=waterlevel_path, oceancurrent_path=oceancurrent_path)
            if type(output_files) is not list:
                output_files = [output_files]


        # Store name and location where file was written
        self._exported_to['input_file'] = []
        for file in output_files:
            self._exported_to['input_file'].append(file)

        if self._input_file_writer._im_silent() or self.dry_run():
            for file in output_files:
                msg.to_file(file)

        return

    def run_model(self, model_executer: ModelExecuter=None,
                input_file: str=None, folder: str=None,
                dateformat: str=None, input_file_extension: str=None,
                dry_run: bool=False, mat_to_nc: bool=False) -> None:
        """Run the model."""
        self._dry_run = dry_run
        self._model_executer = model_executer or self._get_model_executer()
        if self._model_executer is None:
            raise Exception('Define a ModelExecuter!')

        # We always assume that the model is located in the folder the input
        # file was written to

        # Option 1) Use user provided
        # Option 2) Use knowledge of where has been exported
        # Option 3) Use default values to guess where is has previously been exported
        exported_path = Path(self.exported_to('input_file')[0])
        primary_file = input_file or exported_path.name
        primary_folder = folder #or str(exported_path.parent)

        if hasattr(self, '_input_file_writer'):
            extension = input_file_extension or self._input_file_writer._extension()
        else:
            extension = input_file_extension or 'swn'

        file_object = FileNames(format=self._get_default_format(),
                                clean_names=True,
                                dict_of_object_names=self.dict_of_object_names(),
                                start_time=self.start_time,
                                end_time=self.end_time,
                                _filename=primary_file,
                                _folder=primary_folder,
                                _dateformat=dateformat,
                                extension=extension,
                                dnora_obj='input_file')


        msg.header(self._model_executer, "Running model...")
        msg.plain(f"Using input file: {file_object.filepath()}")
        if not self.dry_run():
            self._model_executer(input_file=file_object.filename(), model_folder=file_object.folder())
        else:
            msg.info('Dry run! Model will not run.')
        if mat_to_nc:
            input_file = f'{file_object.folder()}/{self.grid().name()}.mat'
            output_file = f'{file_object.folder()}/{self.grid().name()}.nc'
            convert_swash_mat_to_netcdf(input_file=input_file,output_file=output_file, lon=self.grid().lon_edges(), lat=self.grid().lat_edges(), dt=1)

    def dry_run(self):
        """Checks if method or global ModelRun dryrun is True.
        """
        return self._dry_run or self._global_dry_run

    def _export_object(self, filename: str, folder: str, dateformat: str,
                    writer_function: WritingFunction,
                    dnora_obj: str) -> list[str]:

        # Controls generation of file names using the proper defaults etc.
        file_object = FileNames(format=self._get_default_format(),
                                clean_names=writer_function._clean_filename(),
                                dict_of_object_names=self.dict_of_object_names(),
                                start_time=self.start_time,
                                end_time=self.end_time,
                                _filename=filename,
                                _folder=folder,
                                _dateformat=dateformat,
                                extension=writer_function._extension(),
                                dnora_obj=dnora_obj)
        if self.dry_run():
            msg.info('Dry run! No files will be written.')
            output_files = [file_object.filepath()]
        else:
            # Write the object using the WriterFunction
            file_object.create_folder()
            obj = eval(f'self.{dnora_obj}()')
            output_files = writer_function(obj, file_object.filepath())
            if type(output_files) is not list:
                output_files = [output_files]

        # Store name and location where file was written
        self._exported_to[dnora_obj] = []
        for file in output_files:
            self._exported_to[dnora_obj].append(file)

        if writer_function._im_silent() or self.dry_run():
            for file in output_files:
                msg.to_file(file)

        return output_files

    def plot_spectra(self):
        from dnora.dnplot import dnplot
        __ = dnplot.SpecPlotter().spectra(dict_of_objects=self.dict_of_objects(), plain=False)

    def plot_boundary(self):
        from dnora.dnplot import dnplot
        __ = dnplot.SpecPlotter().boundary(dict_of_objects=self.dict_of_objects(), plain=False)


    def plot_grid(self, grid_plotter: GridPlotter=None, filename: str=None,
                    folder: str=None, dateformat: str=None, plain: bool=False,
                    save_fig: bool=False, show_fig: bool=True) -> dict:
        """Plot the data in the Grid-object, possibly overlaying data from the
        Boundary- and Forcing-objects."""

        if len(self.grid().topo())==0:
            msg.warning('Grid not meshed and nothing to plot!')
            return

        self._grid_plotter = grid_plotter or self._get_grid_plotter()

        if self._grid_plotter is None:
            raise Exception('Define a GridPlotter!')

        figure_dict = self._plot_object(filename=filename, folder=folder,
                            dateformat=dateformat,
                            plotting_function=self._grid_plotter,
                            plain=plain, save_fig=save_fig,
                            show_fig=show_fig, dnora_obj='dnplot_grid')

        return figure_dict

    def plot_topo(self, grid_plotter: GridPlotter=None, filename: str=None,
                folder: str=None, dateformat: str=None, plain: bool=True,
                save_fig: bool=False, show_fig: bool=True) -> dict:
        """Plot the raw data in the Grid-object, possibly overlaying data from the
        Boundary- and Forcing-objects."""


        if len(self.grid().raw_topo()) == 0:
            msg.warning('No topography imported so nothing to plot!')
            return

        self._grid_plotter = grid_plotter or self._get_topo_plotter()

        if self._grid_plotter is None:
            raise Exception('Define a GridPlotter!')

        figure_dict = self._plot_object(filename=filename, folder=folder,
                            dateformat=dateformat,
                            plotting_function=self._grid_plotter,
                            plain=plain, save_fig=save_fig,
                            show_fig=show_fig, dnora_obj='dnplot_topo')
        return figure_dict

    def plot_forcing(self, forcing_plotter: GridPlotter=None, filename: str=None,
                    folder: str=None, dateformat: str=None, plain: bool=False,
                    save_fig: bool=False, show_fig: bool=True) -> dict:
        """Plot the data in the Forcing-object."""

        if self.forcing() is None:
            msg.warning('No forcing data to plot!')
            return

        self._forcing_plotter = forcing_plotter or self._get_forcing_plotter()

        if self._forcing_plotter is None:
            raise Exception('Define a GridPlotter!')

        figure_dict = self._plot_object(filename=filename, folder=folder,
                            dateformat=dateformat,
                            plotting_function=self._forcing_plotter,
                            plain=plain, save_fig=save_fig,
                            show_fig=show_fig, dnora_obj='dnplot_forcing')

        return figure_dict

    def plot_oceancurrent(self, oceancurrent_plotter: GridPlotter=None, filename: str=None,
                    folder: str=None, dateformat: str=None, plain: bool=False,
                    save_fig: bool=False, show_fig: bool=True) -> dict:
        """Plot the data in the OceanCurrent-object."""

        if self.oceancurrent() is None:
            msg.warning('No oceancurrent data to plot!')
            return

        self._oceancurrent_plotter = oceancurrent_plotter or self._get_oceancurrent_plotter()

        if self._oceancurrent_plotter is None:
            raise Exception('Define a GridPlotter!')

        figure_dict = self._plot_object(filename=filename, folder=folder,
                            dateformat=dateformat,
                            plotting_function=self._oceancurrent_plotter,
                            plain=plain, save_fig=save_fig,
                            show_fig=show_fig, dnora_obj='dnplot_oceancurrent')

        return figure_dict


    def _plot_object(self, filename: str, folder: str, dateformat: str,
                    plotting_function: PlottingFunction, plain: bool,
                    save_fig: bool, show_fig: bool, dnora_obj: str) -> dict:
        """Plots a dnora object, e.g. a grid"""

        if filename is not None:
            save_fig = True

        # Controls generation of file names using the proper defaults etc.
        file_object = FileNames(format=self._get_default_format(),
                                clean_names=True,
                                dict_of_object_names=self.dict_of_object_names(),
                                start_time=self.start_time,
                                end_time=self.end_time,
                                _filename=filename,
                                _folder=folder,
                                _dateformat=dateformat,
                                extension=plotting_function._extension(),
                                dnora_obj=dnora_obj)

        file_object.create_folder()

        if dnora_obj in ['dnplot_grid', 'dnplot_forcing', 'dnplot_oceancurrent']:
            figure_dict = plotting_function.grid(dict_of_objects=self.dict_of_objects(), plain=plain)
        elif dnora_obj == 'dnplot_topo':
            figure_dict = plotting_function.topo(dict_of_objects=self.dict_of_objects(), plain=plain)

        if figure_dict is not None:
            fig = figure_dict.get('fig')
        else:
            fig = None
        if fig is not None:
            if save_fig:
                fig.savefig(file_object.filepath(),bbox_inches='tight', dpi=300)
                msg.to_file(file_object.filepath())
            if show_fig:
                fig.show()
        return figure_dict

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

    def oceancurrent(self) -> OceanCurrent:
        """Returns the OceanCurrent object if exists."""
        if hasattr(self, '_oceancurrent'):
            return self._oceancurrent
        else:
            return None

    def waterlevel(self) -> WaterLevel:
        """Returns the waterlevel object if exists."""
        if hasattr(self, '_waterlevel'):
            return self._waterlevel
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

    def input_file(self) -> None:
        """Only defined to have method for all objects"""
        return None


    def dict_of_objects(self) -> dict[str: ModelRun, str: Grid, str: Forcing, str: Boundary, str: Spectra]:
        return {'ModelRun': self, 'Grid': self.grid(), 'Forcing': self.forcing(), 'Boundary': self.boundary(),
        'OceanCurrent': self.oceancurrent() ,  'Spectra': self.spectra()}

    def list_of_objects(self) -> list[ModelRun, Grid, Forcing, Boundary, OceanCurrent, Spectra]:
        """[ModelRun, Boundary] etc."""
        return list(self.dict_of_objects().values())

    def list_of_object_strings(self) -> list[str]:
        """['ModelRun', 'Boundary'] etc."""
        return list(self.dict_of_objects().keys())

    def dict_of_object_names(self) -> dict[str: str]:
        """ {'Boundary': 'NORA3'} etc."""
        d = {}
        for a,b in self.dict_of_objects().items():
            if b is None:
                d[a] = None
            else:
                d[a] = b.name()
        return d

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

    def _get_waterlevel_reader(self) -> WaterLevelReader:
        return None

    def _get_waterlevel_writer(self) -> WaterLevelWriter:
        return None

    def _get_oceancurrent_reader(self) -> OceanCurrentReader:
        return None

    def _get_oceancurrent_writer(self) -> OceanCurrentWriter:
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

    def _get_forcing_plotter(self) -> GridPlotter:
        return ForcingPlotter()

    def _get_oceancurrent_plotter(self) -> GridPlotter:
        return OceanCurrentPlotter()

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
        if self.oceancurrent() is not None:
            lines.append(f'\toceancurrent object: {self.oceancurrent().name()} {self.oceancurrent().size()}')
        else:
            lines.append(f'\toceancurrent: use .import_oceancurrent()')

        return "\n".join(lines)
