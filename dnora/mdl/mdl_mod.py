from __future__ import annotations # For TYPE_CHECKING
from copy import copy
import pandas as pd
import numpy as np
from skeletons import PointSkeleton
# Import objects
from ..grd.grd_mod import Grid


# Import abstract classes and needed instances of them
from ..wnd import Forcing
from ..wnd.read import ForcingReader
from .. import wnd

from ..bnd import Boundary
from ..bnd.read import BoundaryReader
from ..bnd.pick import PointPicker
from .. import bnd


from ..spc import Spectra
from ..spc.read import SpectralReader
from .. import spc

from ..wsr import WaveSeries
from ..wsr.read import WaveSeriesReader, SpectraToWaveSeries
from .. import wsr

from skeletons.datavar_factory import add_datavar
# Import default values and aux_funcsiliry functions
from .. import msg
#from ..cacher.cache_decorator import cached_reader

class ModelRun:
    def __init__(self, grid: Grid, start_time: str='1970-01-01T00:00',
    end_time: str='2030-12-31T23:59', name: str='AnonymousModelRun',
    dry_run: bool=False):
        self.name = copy(name)
        self._grid = copy(grid)
        self._grid.exported_to = [None]
        self._time = pd.date_range(start_time, end_time, freq='H')
        self._global_dry_run = dry_run
        self._dry_run = False  # Set by methods

    #@cached_reader('Forcing', wnd.read.DnoraNc)
    def import_forcing(self, forcing_reader: ForcingReader=None,
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

        if not self.dry_run():
            time, u, v, lon, lat, x, y, attributes = forcing_reader(self.grid(), self.start_time(), self.end_time(), **kwargs)
            
            self._forcing = Forcing(lon=lon, lat=lat, x=x, y=y, time=time)
            x = x or lon
            y = y or lat
            self.forcing().set_spacing(nx=len(x), ny=len(y))
            
            self.forcing().name = name 
            self.forcing().set_u(u)
            self.forcing().set_v(v)
            self.forcing().set_metadata(attributes)
        else:
            msg.info('Dry run! No forcing will be imported.')


    def import_boundary(self, boundary_reader: BoundaryReader=None,
                        point_picker: PointPicker=None, name: str=None,
                        dry_run: bool=False, source: str='remote',
                       **kwargs):
        """Imports boundary spectra. Which spectra to choose spatically
        are determined by the point_picker.

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


        msg.header(boundary_reader, "Reading coordinates of spectra...")
        lon_all, lat_all, x_all, y_all = boundary_reader.get_coordinates(self.grid(), self.start_time())
        all_points = PointSkeleton(lon=lon_all, lat=lat_all, x=x_all, y=y_all)
        
        if np.all(np.logical_not(self.grid().boundary_mask())):
            boundary_points = None
        else:
            boundary_points = PointSkeleton.from_skeleton(self.grid(), mask=self.grid().boundary_mask())
        
        msg.header(point_picker, "Choosing boundary spectra...")
        inds = point_picker(self.grid(), all_points, selected_points=boundary_points, **kwargs)
        if len(inds) < 1:
            msg.warning("PointPicker didn't find any points. Aborting import of boundary.")
            return

        # Main reading happens here
        msg.header(boundary_reader, "Loading boundary spectra...")

        time, freq, dirs, spec, lon, lat, x, y, metadata = boundary_reader(self.grid(), self.start_time(), self.end_time(), inds)
        self._boundary = Boundary(x=x, y=y, lon=lon, lat=lat, time=time, freq=freq, dirs=dirs)

        self.boundary().set_spec(spec)
        self.boundary().set_metadata(metadata)
        # E.g. are the spectra oceanic convention etc.
        self.boundary()._convention = boundary_reader.convention()

        self.boundary().set_metadata({'spectral_convention': self.boundary().convention().value}, append=True)

        if boundary_reader.post_processing() is not None:
            self.boundary().process_boundary(boundary_reader.post_processing())
        return



    def import_spectra(self, spectral_reader: SpectralReader=None,
                        point_picker: PointPicker=None, name: str=None,
                        dry_run: bool=False, source: str='remote',
                       **kwargs):
        """Imports spectra. Which spectra to choose spatically
        are determined by the point_picker.

        source = 'remote' (default) / '<folder>' / 'met'

        The implementation of this is up to the SpectralReader, and all options might not be functional.
        'met' options will only work in MET Norway internal networks.
        
        To import local netcdf files saved in DNORA format (by write_cache=True), use read_cache=True.
        """

        self._dry_run = dry_run
        spectral_reader = spectral_reader or self._get_spectral_reader()
        point_picker = point_picker or self._get_point_picker()

        # This is to allow importing from cache using only a name
        if spectral_reader is None:
            raise Exception('Define a SpectralReader!')
        if point_picker is None:
            raise Exception('Define a PointPicker!')

        msg.header(spectral_reader, "Reading coordinates of spectra...")
        lon_all, lat_all, x_all, y_all = spectral_reader.get_coordinates(self.grid(), self.start_time())
        all_points = PointSkeleton(lon=lon_all, lat=lat_all, x=x_all, y=y_all)
        
        if np.all(np.logical_not(self.grid().boundary_mask())):
            boundary_points = None
        else:
            boundary_points = PointSkeleton.from_skeleton(self.grid(), mask=self.grid().boundary_mask())
        
        msg.header(point_picker, "Choosing spectra...")
        inds = point_picker(self.grid(), all_points, selected_points=boundary_points, **kwargs)

        msg.header(spectral_reader, "Loading omnidirectional spectra...")
        time, freq, spec, mdir, spr, lon, lat, x, y, metadata = spectral_reader(self.grid(), self.start_time(), self.end_time(), inds)

        self._spectra = Spectra(x=x, y=y, lon=lon, lat=lat, time=time, freq=freq)

        self.spectra().set_spec(spec)
        self.spectra().set_mdir(mdir)
        self.spectra().set_spr(spr)
        
        self.spectra().set_metadata(metadata)

        # E.g. are the spectra oceanic convention etc.
        self.spectra()._convention = spectral_reader.convention()
        self.spectra().set_metadata({'spectral_convention': self.spectra().convention().value}, append=True)
   

    def import_waveseries(self, waveseries_reader: WaveSeriesReader=None,
                        point_picker: PointPicker=None, name: str=None,
                        dry_run: bool=False, source: str='remote',
                       **kwargs):

        msg.header(waveseries_reader, "Reading coordinates of WaveSeries...")
        lon_all, lat_all, x_all, y_all = waveseries_reader.get_coordinates(self.grid(), self.start_time())

        all_points = PointSkeleton(lon=lon_all, lat=lat_all, x=x_all, y=y_all)
        
        if np.all(np.logical_not(self.grid().boundary_mask())):
            boundary_points = None
        else:
            boundary_points = PointSkeleton.from_skeleton(self.grid(), mask=self.grid().boundary_mask())
        
        msg.header(point_picker, "Choosing wave series points...")
        inds = point_picker(self.grid(), all_points, selected_points=boundary_points, **kwargs)

        msg.header(waveseries_reader, "Loading wave series data...")
        time, data_dict, lon, lat, x, y, metadata = waveseries_reader(self.grid(), self.start_time(), self.end_time(), inds, **kwargs)

        self._waveseries = WaveSeries(x, y, lon, lat, time=time, name=name)

        for wp, data in data_dict.items():
            self._waveseries = add_datavar(wp.name(), append=True)(self.waveseries()) # Creates .hs() etc. methods
            self.waveseries()._update_datavar(wp.name(), data)
            self.waveseries().set_metadata({'name': wp.name(), 'unit': f'{wp.unit()}', 'standard_name': wp.standard_name()}, data_array_name=wp.name())
                       
        self.waveseries().set_metadata(metadata) # Global attributes

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

    def dry_run(self):
        """Checks if method or global ModelRun dryrun is True.
        """
        return self._dry_run or self._global_dry_run


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

    def time(self, crop: bool=False):
        """Returns times of ModelRun
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
    
    def start_time(self, crop: bool=False):
        """Returns start time of ModelRun
        crop = True: Give the period that is covered by all objects (Forcing etc.)"""
        return self.time(crop=crop)[0]

    def end_time(self, crop: bool=False):
        """Returns start time of ModelRun
        crop = True: Give the period that is covered by all objects (Forcing etc.)"""
        return self.time(crop=crop)[-1]

    def _get_forcing_reader(self) -> ForcingReader:
        return None
    
    def _get_boundary_reader(self) -> BoundaryReader:
        return None
    
    def _get_point_picker(self) -> PointPicker:
        return None

    def _get_spectral_reader(self) -> SpectralReader:
        return None
    
    def _get_waveseries_reader(self) -> WaveSeriesReader:
        return None