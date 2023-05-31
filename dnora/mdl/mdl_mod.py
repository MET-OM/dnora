from __future__ import annotations # For TYPE_CHECKING
from copy import copy
import re
from pathlib import Path
import yaml
import numpy as np
import pandas as pd
# Import objects
from ..grd.grd_mod import Grid

from ..spc.spc_mod import Spectra
from ..wsr.wsr_mod import WaveSeries
# Import abstract classes and needed instances of them
from .. import bnd, wnd, spc, wsr

from ..bnd.conventions import SpectralConvention
from ..wsr.read import WaveSeriesReader, SpectraToWaveSeries

from ..grd.write import GridWriter
from ..grd.process import GridProcessor, TrivialFilter
#from ..trg.write import TrGridWriter

from ..dnplot.dnplot import GridPlotter, TopoPlotter, ForcingPlotter
from ..inp.inp import InputFileWriter
from ..run import ModelExecuter

from ..file_module import FileNames
from typing import Union
# Import default values and aux_funcsiliry functions
from .. import msg
from ..cacher.cacher import Cacher
from ..cacher.cache_decorator import cached_reader
import os, glob
from .. import file_module
from ..converters import convert_swash_mat_to_netcdf
WritingFunction = Union[GridWriter, bnd.write.BoundaryWriter, spc.write.SpectralWriter, wnd.write.ForcingWriter]
#PlottingFunction = Union[GridPlotter]

class ModelRun:
    def __init__(self, grid: Grid, start_time: str='1970-01-01T00:00',
    end_time: str='2030-12-31T23:59', name: str='AnonymousModelRun',
    dry_run: bool=False):
        self.name = copy(name)
        self._grid = copy(grid)
        self._grid.exported_to = [None]
        self._time = pd.date_range(start_time, end_time, freq='H')
        #self.start_time = copy(start_time)
        #self.end_time = copy(end_time)
        self._global_dry_run = dry_run
        self._dry_run = False  # Set by methods

    @cached_reader('Boundary', bnd.read.DnoraNc)
    def import_boundary(self, boundary_reader: bnd.read.BoundaryReader=None,
                        point_picker: PointPicker=None, name: str=None,
                        dry_run: bool=False,
                        source: str='remote',
                        **kwargs) -> None:
        """Imports boundary spectra.

        source = 'remote' (default) / '<folder>' / 'met'

        The implementation of this is up to the BoundaryReader, and all options might not be functional.
        'met' options will only work in MET Norway internal networks.

        To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
        """

        self._dry_run = dry_run
        boundary_reader = boundary_reader or self._get_boundary_reader()
        point_picker = point_picker or self._get_point_picker()

        # This is to allow importing from cache using only a name
        if boundary_reader is None:
            raise Exception('Define a BoundaryReader!')
        if point_picker is None:
            raise Exception('Define a PointPicker!')

        name = name or boundary_reader.name()
        boundary_reader.set_source(source)

        if name is None:
            raise ValueError('Provide either a name or a BoundaryReader that will then define the name!')

        self._boundary = bnd.Boundary(grid=self.grid(), name=name)
        self._boundary.exported_to = [None]

        # Import the boundary spectra into the Boundary-object
        if not self.dry_run():
            self.boundary().import_boundary(start_time=self.time()[0],
                                            end_time=self.time()[-1],
                                            boundary_reader=boundary_reader,
                                            point_picker=point_picker,
                                            **kwargs)
        else:
            msg.info('Dry run! No boundary spectra will be imported.')

    @cached_reader('Forcing', wnd.read.DnoraNc)
    def import_forcing(self, forcing_reader: wnd.read.ForcingReader=None,
                        name: str=None, dry_run: bool=False,
                        source: str='remote',
                        **kwargs) -> None:
        """Imports wind forcing.

        source = 'remote' (default) / '<folder>' / 'met'

        The implementation of this is up to the ForcingReader, and all options might not be functional.
        'met' options will only work in MET Norway internal networks.

        To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
        """
        self._dry_run = dry_run

        forcing_reader = forcing_reader or self._get_forcing_reader()

        # This is to allow importing from cache using only a name
        if forcing_reader is None:
            raise Exception('Define a ForcingReader!')

        name = name or forcing_reader.name()
        forcing_reader.set_source(source)

        if name is None:
            raise ValueError('Provide either a name or a ForcingReader that will then define the name!')

        self._forcing = wnd.Forcing(grid=self.grid(), name=name)
        self._forcing.exported_to = [None]

        # Import the forcing data into the Forcing-object
        if not self.dry_run():
            self.forcing().import_forcing(start_time=self.time()[0],
                                        end_time=self.time()[-1],
                                        forcing_reader=forcing_reader,
                                        **kwargs)
        else:
            msg.info('Dry run! No forcing will be imported.')

    @cached_reader('Spectra', spc.read.DnoraNc)
    def import_spectra(self, spectral_reader: SpectralReader=None,
                        point_picker: PointPicker=None,
                        name: str=None, dry_run: bool=False,
                        source: str='remote', files: list[str] = None, **kwargs) -> None:
        """Imports omnidirectional spectra.

        source = 'remote' (default) / '<folder>' / 'met'

        The implementation of this is up to the SpectralReader, and all options might not be functional.
        'met' options will only work in MET Norway internal networks.

        To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
        """

        self._dry_run = dry_run

        if files is not None:
            if not isinstance(files, list):
                files = [files]
            spectral_reader = spc.read.DnoraNc(files=files)
        spectral_reader = spectral_reader or self._get_spectra_reader()
        point_picker = point_picker or self._get_point_picker()

        # This is to allow importing from cache using only a name
        if spectral_reader is None:
            raise Exception('Define a SpectralReader!')
        if point_picker is None:
            raise Exception('Define a PointPicker!')

        name = name or spectral_reader.name()
        spectral_reader.set_source(source)

        if name is None:
            raise ValueError('Provide either a name or a SpectralReader that will then define the name!')

        # Create spectral object
        self._spectra = spc.Spectra(grid=self.grid(), name=name)
        self._spectra.exported_to = [None]

        # Import the forcing data into the Forcing-object
        if not self.dry_run():
            self.spectra().import_spectra(start_time=self.time()[0],
                                        end_time=self.time()[-1],
                                        spectral_reader=spectral_reader,
                                        point_picker=point_picker,
                                        **kwargs)
        else:
            msg.info('Dry run! No omnidirectional spectra will be imported.')

    @cached_reader('WaveSeries', wsr.read.DnoraNc)
    def import_waveseries(self, waveseries_reader: WavesSeriesReader=None,
                        point_picker: PointPicker=None,
                        name: str=None, dry_run: bool=False,
                        source: str='remote', files: list[str] = None, **kwargs) -> None:

        """Imports wave timeseries data.

        source = 'remote' (default) / '<folder>' / 'met'

        The implementation of this is up to the WaveSeriesReader, and all options might not be functional.
        'met' options will only work in MET Norway internal networks.

        To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
        """
        self._dry_run = dry_run
        if files is not None:
            if not isinstance(files, list):
                files = [files]
            waveseries_reader = wsr.read.DnoraNc(files=files)
        waveseries_reader = waveseries_reader or self._get_waveseries_reader()

        point_picker = point_picker or self._get_point_picker()
        # This is to allow importing from cache using only a name
        if waveseries_reader is None:
            raise Exception('Define a WaveSeriesReader!')
        if point_picker is None:
            raise Exception('Define a PointPicker!')

        name = name or waveseries_reader.name()
        waveseries_reader.set_source(source)

        if name is None:
            raise ValueError('Provide either a name or a WaveSeriesReader that will then define the name!')

        self._waveseries = WaveSeries(grid=self.grid(), name=name)
        self._waveseries.exported_to = [None]

        # Import the forcing data into the Forcing-object
        if not self.dry_run():
            self.waveseries().import_waveseries(start_time=self.time()[0],
                                        end_time=self.time()[-1],
                                        waveseries_reader=waveseries_reader,
                                        point_picker=point_picker,
                                        **kwargs)
        else:
            msg.info('Dry run! No wave data will be imported.')

    def cache_boundary(self):
        """Writes existing data to cached files."""
        self.export_boundary(boundary_writer=bnd.write.DnoraNc(), format='Cache')

    def cache_spectra(self):
        """Writes existing data to cached files."""
        self.export_spectra(spectral_writer=spc.write.DnoraNc(), format='Cache')

    def cache_forcing(self):
        """Writes existing data to cached files."""
        self.export_forcing(forcing_writer=wnd.write.DnoraNc(), format='Cache')

    def cache_waveseries(self):
        """Writes existing data to cached files."""
        self.export_waveseries(waveseries_writer=wsr.write.DnoraNc(), format='Cache')

    def boundary_to_spectra(self, dry_run: bool=False, name :str=None,
                            write_cache=False, **kwargs):
        self._dry_run = dry_run
        if self.boundary() is None:
            msg.warning('No Boundary to convert to Spectra!')

        spectral_reader = spc.read.BoundaryToSpectra(self.boundary())
        msg.header(spectral_reader, 'Converting the boundary spectra to omnidirectional spectra...')
        name = self.boundary().name

        if not self.dry_run():
            self.import_spectra(spectral_reader=spectral_reader,
                                point_picker=bnd.pick.TrivialPicker(),
                                name=name,
                                write_cache=write_cache, **kwargs)
        else:
            msg.info('Dry run! No boundary will not be converted to spectra.')

    def spectra_to_waveseries(self, dry_run: bool=False, write_cache=False,
                                freq: tuple=(0, 10_000), **kwargs):
        self._dry_run = dry_run
        if self.spectra() is None:
            msg.warning('No Spectra to convert to WaveSeries!')
            return

        waveseries_reader = SpectraToWaveSeries(self.spectra(), freq)
        msg.header(waveseries_reader, 'Converting the spectra to wave series data...')
        name = self.spectra().name
        if not self.dry_run():
            self.import_waveseries(waveseries_reader=waveseries_reader,
                                    point_picker=bnd.pick.TrivialPicker(),
                                    name=name,
                                    write_cache=write_cache, **kwargs)
        else:
            msg.info('Dry run! No boundary will not be converted to spectra.')

    def boundary_to_waveseries(self, dry_run: bool=False, write_cache=False,
                                freq: tuple=(0, 10_000), **kwargs):
        self.boundary_to_spectra(dry_run=dry_run, write_cache=write_cache, **kwargs)
        self.spectra_to_waveseries(dry_run=dry_run, write_cache=write_cache, freq=freq, **kwargs)

    def set_spectral_grid(self, freq0: float=0.04118, nfreq: int=32, ndir: int=36, finc: float=1.1, dirshift: float=0.):
        """Sets spectral grid for model run. Will be used to write input files."""
        self._spectral_grid = bnd.Boundary(grid=None, name='spectral_grid')
        freq = np.array([freq0*finc**n for n in np.linspace(0,nfreq-1,nfreq)])
        dirs = np.linspace(0,360,ndir+1)[0:-1]+dirshift
        time = self.time()
        
        self._spectral_grid._init_structure(x=0, y=0, lon=None, lat=None, time=time, freq=freq, dirs=dirs)
        

    def export_grid(self, grid_writer: GridWriter=None,
                    filename: str=None, folder: str=None, dateformat: str=None,
                    format: str=None, dry_run=False) -> None:
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

        msg.header(self._grid_writer, f"Writing grid topography from {self.grid().name}")

        output_files = self._export_object('Grid',filename, folder, dateformat,
                            writer_function=self._grid_writer, format=format)

        # Write status file
        #infofilename = str(Path(output_files[0]).with_suffix('')) + '_dnora_info.txt'
        #self.grid().write_status(filename=infofilename)

    def export_boundary(self, boundary_writer: BoundaryWriter=None,
                        filename: str=None, folder: str=None,
                        dateformat: str=None, format: str=None, dry_run=False) -> None:
        """Writes the spectra in the Boundary-object to an external source, e.g.
        a file."""
        self._dry_run = dry_run
        if self.boundary() is None and not self.dry_run(dry_run):
            raise Exception('Import boundary before exporting!')

        self._boundary_writer = boundary_writer or self._get_boundary_writer()
        if self._boundary_writer is None:
            raise Exception('Define a BoundaryWriter!')

        if self.boundary() is None:
            msg.header(self._boundary_writer, f"Writing boundary spectra from DryRunBoundary")
        else:
            msg.header(self._boundary_writer, f"Writing boundary spectra from {self.boundary().name}")

        if not self.dry_run():
            self.boundary()._set_convention(self._boundary_writer.convention())

        __ = self._export_object('Boundary', filename, folder, dateformat,
                            writer_function=self._boundary_writer, format=format)

    def export_spectra(self, spectral_writer: SpectralWriter=None,
                        filename: str=None, folder: str=None,
                        dateformat: str=None, format: str=None, dry_run=False) -> None:
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
            msg.header(self._spectral_writer, f"Writing omnidirectional spectra from {self.spectra().name}")

        if not self.dry_run():
            self.spectra()._set_convention(self._spectral_writer.convention())

        # Replace #Spectra etc and add file extension
        __ = self._export_object('Spectra', filename, folder, dateformat,
                            writer_function=self._spectral_writer, format=format)

    def export_waveseries(self, waveseries_writer: WaveSeriesWriter=None,
                        filename: str=None, folder: str=None,
                        dateformat: str=None, format: str=None, dry_run=False) -> None:
        """Writes the data of the WaveSeries-object to an external source, e.g.
        a file."""
        self._dry_run = dry_run
        if self.waveseries() is None and not self.dry_run():
            raise Exception('Import waveseries data before exporting!')

        self._waveseries_writer = waveseries_writer or self._get_waveseries_writer()
        if self._waveseries_writer is None:
            raise Exception('Define a WaveSeriesWriter!')

        if self.waveseries() is None:
            msg.header(self._waveseries_writer, f"Writing wave series data from DryRunSpectra")
        else:
            msg.header(self._waveseries_writer, f"Writing wave series data from {self.waveseries().name}")

        # Replace #Spectra etc and add file extension
        __ = self._export_object('WaveSeries', filename, folder, dateformat,
                            writer_function=self._waveseries_writer, format=format)

    def export_forcing(self, forcing_writer: ForcingWriter=None,
                        filename: str=None, folder: str=None,
                         dateformat: str=None, format: str=None, dry_run=False) -> None:
        """Writes the forcing data in the Forcing-object to an external source,
        e.g. a file."""
        self._dry_run = dry_run
        if self.forcing() is None and not self.dry_run():
            raise Exception('Import forcing before exporting!')

        self._forcing_writer = forcing_writer or self._get_forcing_writer()

        if self._forcing_writer is None:
            raise Exception('Define a ForcingWriter!')

        if self.forcing() is not None:
            msg.header(self._forcing_writer, f"Writing wind forcing from {self.forcing().name}")
        else:
            msg.header(self._forcing_writer, f"Writing wind forcing from DryRunForcing")


        __ = self._export_object('Forcing', filename, folder, dateformat,
                            writer_function=self._forcing_writer, format=format)

    def write_input_file(self, input_file_writer: InputFileWriter=None,
                        filename=None, folder=None, dateformat=None,
                        grid_path: str=None, forcing_path: str=None,
                        boundary_path: str=None, start_time: str=None,
                        end_time: str=None, dry_run=False) -> None:
        """Writes the grid data in the Grid-object to an external source,
        e.g. a file."""
        self._dry_run = dry_run
        self._input_file_writer = input_file_writer or self._get_input_file_writer()
        if self._input_file_writer is None:
            raise Exception('Define an InputFileWriter!')


        msg.header(self._input_file_writer, "Writing model input file...")

        # Controls generation of file names using the proper defaults etc.
        file_object = FileNames(dict_of_objects=self.dict_of_objects(),
                                filename=filename,
                                folder=folder,
                                dateformat=dateformat,
                                extension=self._input_file_writer._extension(),
                                obj_type='input_file',
                                edge_object='Grid')

        file_object.create_folder()


        if self.dry_run():
            msg.info('Dry run! No files will be written.')
            output_files = [file_object.get_filepath()]
        else:
            # Write the grid using the InputFileWriter object
            output_files = self._input_file_writer(dict_of_objects=self.dict_of_objects(),
                            filename=file_object.get_filepath())
            if type(output_files) is not list:
                output_files = [output_files]


        # # Store name and location where file was written
        # self._exported_to['input_file'] = []
        # for file in output_files:
        #     self._exported_to['input_file'].append(file)

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

        file_object = FileNames(dict_of_objects=self.dict_of_objects(),
                                filename=primary_file,
                                folder=primary_folder,
                                dateformat=dateformat,
                                extension=extension,
                                obj_type='model_executer',
                                edge_object='Grid')


        msg.header(self._model_executer, "Running model...")
        msg.plain(f"Using input file: {file_object.get_filepath()}")
        if not self.dry_run():
            self._model_executer(input_file=file_object.filename(), model_folder=file_object.folder())
        else:
            msg.info('Dry run! Model will not run.')
        if mat_to_nc:
            input_file = f'{file_object.folder()}/{self.grid().name}.mat'
            output_file = f'{file_object.folder()}/{self.grid().name}.nc'
            convert_swash_mat_to_netcdf(input_file=input_file,output_file=output_file, lon=self.grid().lon_edges(), lat=self.grid().lat_edges(), dt=1)

    def dry_run(self):
        """Checks if method or global ModelRun dryrun is True.
        """
        return self._dry_run or self._global_dry_run

    def _export_object(self, obj_type, filename: str, folder: str, dateformat: str,
                    writer_function: WritingFunction, format: str) -> list[str]:

        # Controls generation of file names using the proper defaults etc.
        format = format or self._get_default_format()
        file_object = FileNames(format=format,
                                obj_type=obj_type,
                                dict_of_objects=self.dict_of_objects(),
                                filename=filename,
                                folder=folder,
                                dateformat=dateformat,
                                extension=writer_function._extension())
        if self.dry_run():
            msg.info('Dry run! No files will be written.')
            output_files = [file_object.get_filepath()]
        else:
            # Write the object using the WriterFunction
            file_object.create_folder()
            output_files = writer_function(self.dict_of_objects(), file_object)
            if type(output_files) is not list:
                output_files = [output_files]

        # Store name and location where file was written
        self.dict_of_objects().get(obj_type).exported_to = output_files
        # for file in output_files:
        #     self._exported_to[obj_type].append(file)

        if writer_function._im_silent() or self.dry_run():
            for file in output_files:
                msg.to_file(file)

        return output_files

    # def plot_spectra(self):
    #     from dnora.dnplot import dnplot
    #     __ = dnplot.SpecPlotter().spectra(dict_of_objects=self.dict_of_objects(), plain=False)

    # def plot_boundary(self):
    #     from dnora.dnplot import dnplot
    #     __ = dnplot.SpecPlotter().boundary(dict_of_objects=self.dict_of_objects(), plain=False)


    # def plot_grid(self, grid_plotter: GridPlotter=None, filename: str=None,
    #                 folder: str=None, dateformat: str=None, plain: bool=False,
    #                 save_fig: bool=False, show_fig: bool=True) -> dict:
    #     """Plot the data in the Grid-object, possibly overlaying data from the
    #     Boundary- and Forcing-objects."""

    #     if len(self.grid().topo())==0:
    #         msg.warning('Grid not meshed and nothing to plot!')
    #         return

    #     self._grid_plotter = grid_plotter or self._get_grid_plotter()

    #     if self._grid_plotter is None:
    #         raise Exception('Define a GridPlotter!')

    #     figure_dict = self._plot_object(obj_type='Grid',
    #                         filename=filename, folder=folder,
    #                         dateformat=dateformat,
    #                         plotting_function=self._grid_plotter,
    #                         plain=plain, save_fig=save_fig,
    #                         show_fig=show_fig)

    #     return figure_dict

    # def plot_topo(self, grid_plotter: GridPlotter=None, filename: str=None,
    #             folder: str=None, dateformat: str=None, plain: bool=True,
    #             save_fig: bool=False, show_fig: bool=True) -> dict:
    #     """Plot the raw data in the Grid-object, possibly overlaying data from the
    #     Boundary- and Forcing-objects."""


    #     if len(self.grid().raw_topo()) == 0:
    #         msg.warning('No topography imported so nothing to plot!')
    #         return

    #     self._grid_plotter = grid_plotter or self._get_topo_plotter()

    #     if self._grid_plotter is None:
    #         raise Exception('Define a GridPlotter!')

    #     figure_dict = self._plot_object(filename=filename, folder=folder,
    #                         dateformat=dateformat,
    #                         plotting_function=self._grid_plotter,
    #                         plain=plain, save_fig=save_fig,
    #                         show_fig=show_fig, dnora_obj='Topo')
    #     return figure_dict

    # def plot_forcing(self, forcing_plotter: GridPlotter=None, filename: str=None,
    #                 folder: str=None, dateformat: str=None, plain: bool=False,
    #                 save_fig: bool=False, show_fig: bool=True) -> dict:
    #     """Plot the data in the Forcing-object."""

    #     if self.forcing() is None:
    #         msg.warning('No forcing data to plot!')
    #         return

    #     self._forcing_plotter = forcing_plotter or self._get_forcing_plotter()

    #     if self._forcing_plotter is None:
    #         raise Exception('Define a GridPlotter!')

    #     figure_dict = self._plot_object(obj_type='Forcing', filename=filename,
    #                         folder=folder,
    #                         dateformat=dateformat,
    #                         plotting_function=self._forcing_plotter,
    #                         plain=plain, save_fig=save_fig,
    #                         show_fig=show_fig)

    #     return figure_dict

    # def _plot_object(self, obj_type:str, filename: str, folder: str, dateformat: str,
    #                 plotting_function: PlottingFunction, plain: bool,
    #                 save_fig: bool, show_fig: bool) -> dict:
    #     """Plots a dnora object, e.g. a grid"""

    #     if filename is not None:
    #         save_fig = True

    #     # Controls generation of file names using the proper defaults etc.
    #     file_object = FileNames(dict_of_objects=self.dict_of_objects(),
    #                             filename=filename,
    #                             folder=folder,
    #                             dateformat=dateformat,
    #                             extension=plotting_function._extension(),
    #                             obj_type=obj_type,
    #                             edge_object='Grid')
    #     # file_object = FileNames(format=self._get_default_format(),
    #     #                         dnora_obj=dnora_obj,
    #     #                         clean_names=True,
    #     #                         dict_of_object_names=self.dict_of_object_names(),
    #     #                         start_time=self.start_time,
    #     #                         end_time=self.end_time,
    #     #                         _filename=filename,
    #     #                         _folder=folder,
    #     #                         _dateformat=dateformat,
    #     #                         extension=plotting_function._extension())

    #     file_object.create_folder(key='plotfolder')
    #     dnora_obj = self.dict_of_objects()[obj_type]
    #     if dnora_obj.is_gridded():
    #         figure_dict = plotting_function.gridded(dict_of_objects=self.dict_of_objects(), plain=plain)
    #     else:
    #         figure_dict = plotting_function.ungridded(dict_of_objects=self.dict_of_objects(), plain=plain)

    #     if figure_dict is not None:
    #         fig = figure_dict.get('fig')
    #     else:
    #         fig = None
    #     if fig is not None:
    #         if save_fig:
    #             fig.savefig(file_object.get_filepath(key='plotname'),bbox_inches='tight', dpi=300)
    #             msg.to_file(file_object.get_filepath(key='plotname'))
    #         if show_fig:
    #             fig.show()
    #     return figure_dict
    
    def time(self, crop: bool=False):
        """Returns start and end time of ModelRun
        crop = True: Give the period that is covered by all objects (Forcing etc.)"""
        t0 = self._time[0]
        t1 = self._time[-1]

        if crop:
            for dnora_obj in self.list_of_objects():
                time = dnora_obj.time()
                if time[0] is not None:
                    t0 = pd.to_datetime([t0, time[0]]).max()
                if time[-1] is not None:
                    t1 = pd.to_datetime([t1, time[-1]]).min()
        time = pd.date_range(t0, t1, freq='H')
        return time[::len(time)-1]
        

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

    def waveseries(self) -> WaveSeries:
        """Returns the wave series object if exists."""
        if hasattr(self, '_waveseries'):
            return self._waveseries
        else:
            return None

    def topo(self) -> Grid:
        """Returns the raw topography object if exists."""
        if hasattr(self.grid(), '_raw'):
            return self.grid().raw()
        else:
            return None

    def spectral_grid(self) -> Boundary:
        """Returns spectral grid if exists"""
        if hasattr(self, '_spectral_grid'):
            return self._spectral_grid
        else:
            return None

    def input_file(self) -> None:
        """Only defined to have method for all objects"""
        return None


    def dict_of_objects(self) -> dict[str: ModelRun, str: Grid, str: Forcing, str: Boundary, str: Spectra]:
        return {'ModelRun': self, 'Grid': self.grid(), 'Topo': self.topo(), 'Forcing': self.forcing(),
                'Boundary': self.boundary(), 'Spectra': self.spectra(), 'WaveSeries': self.waveseries(),
                'SpectralGrid': self.spectral_grid()}

    def list_of_objects(self) -> list[ModelRun, Grid, Forcing, Boundary, Spectra]:
        """[ModelRun, Boundary] etc."""
        return [x for x in list(self.dict_of_objects().values()) if x is not None]

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
                d[a] = b.name
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

    def _get_spectra_reader(self) -> SpectralReader:
        return None

    def _get_waveseries_reader(self) -> WaveSeriesReader:
        return None

    def _get_spectral_writer(self) -> SpectralWriter:
        return None

    def _get_waveseries_writer(self) -> WaveSeriesWriter:
        return None

    def _get_grid_writer(self) -> GridWriter:
        return None

    def _get_point_picker(self) -> PointPicker:
        return bnd.pick.TrivialPicker()

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

    # def __repr__(self):
    #     lines = [f"<dnora ModelRun object> ({type(self).__name__})", f"  Name: {self.name}"]
    #
    #     lines.append(f'  Covering time: {self.start_time} - {self.end_time}')
    #
    #     # lines.append(f"  Contains:")
    #     # if self.grid().structured():
    #     #     gridtype = 'structured'
    #     # else:
    #     #     gridtype = 'unstructured'
    #
    #     lines.append(f'\tgrid object ({gridtype}): {self.grid().name} {self.grid().size()}')
    #     if self.forcing() is not None:
    #         lines.append(f'\tforcing object: {self.forcing().name} {self.forcing().size()}')
    #     else:
    #         lines.append(f'\tforcing: use .import_forcing()')
    #     if self.boundary() is not None:
    #         lines.append(f'\tboundary object: {self.boundary().name} {self.boundary().size()}')
    #     else:
    #         lines.append(f'\tboundary: use .import_boundary()')
    #
    #
    #     return "\n".join(lines)
