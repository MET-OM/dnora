from __future__ import annotations
from abc import ABC, abstractmethod

# Import objects

from dnora.type_manager.data_sources import DataSource
from dnora.type_manager.spectral_conventions import SpectralConvention

import re
from typing import TYPE_CHECKING, Union, Optional

if TYPE_CHECKING:
    from type_manager.dnora_types import DnoraDataType
    from dnora.grid import Grid
from dnora.cacher.caching_strategies import CachingStrategy

from dnora.process.gridded import FillNaNs
from dnora.utils.io import get_url


class DataReader(ABC):
    """Reads forcing data from some source and provide it to the object.

    The area is defined from the Grid object that is passed.
    """

    _default_folders = {}
    _default_filename = None  # If we have a source independent filename
    _default_filenames = {}  # Possible source-dependent filenames

    def post_processing(self):
        """Class to use for post processing of data"""
        return FillNaNs(0)

    @staticmethod
    def caching_strategy() -> CachingStrategy:
        """Defines what caching strategy to use"""
        return CachingStrategy.PatchInTime

    @abstractmethod
    def __call__(
        self,
        obj_type: DnoraDataType,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        folder: str,
        filename: str,
        inds: list[int],
        obj_data_vars: list[str],
        **kwargs,
    ):
        """Reads in the forcing witih grid and between start_time and end_time.

        The variables needed to be returned are:

        - coord_dict: Dictionary containing at least lon/lat or x/y. Can also contain e.g. time
            E.g. coord_dict = {'lon': np_array_of_longitude, 'lat': np_array_of_latitude}

        - data_dict:  Dictionary containing the data variables. Key-names depend on object type

            E.g. data_dict = {'u': np_array_of_u_wind_component, 'v': np_array_of_v_wind_component}

        - meta_dict: dict{key, value} will be set as attributes of the xr.Dataset. Can be empty.
            E.g. meta_dict = {'model': 'MEPS', 'institute': 'MET Norway'}

        - metaparameter_dict: dict{key, values} couples the data variables to a metaparameter (see dnora.metaparameter.parameters).
            E.g. metaparameter_dict = {'u': XWind, 'v': YWind}
        Can be empty and is automatically detected from the attributes of he dnora object class.
        Only needs to be specified in e.g. WaveSeries where the data variables are not fixed and created dynamically.

        """

        pass

    def name(self) -> str:
        return type(self).__name__

    def __str__(self) -> str:
        return self.name()

    default_data_source = DataSource.UNDEFINED

    def set_default_data_source(self, source):
        self._default_data_source = source

    def default_data_source(self) -> DataSource:
        return self._default_data_source

    def _folder(
        self,
        folder: str,
        source: DataSource,
        tile: Optional[str] = None,
        tile_name: Optional[str] = None,
        strict: bool = True,
    ) -> str:
        default_folder = self._default_folders.get(source)
        if source in [DataSource.REMOTE, DataSource.LOCAL]:
            folder = folder or default_folder
            if folder is not None:
                if tile_name is not None:
                    folder = re.sub("#TILENAME", tile_name, folder)
                if tile is not None:
                    folder = re.sub("#TILE", tile, folder)
                return folder
        else:
            if folder is None:
                raise ValueError(
                    f"Define an environmental variable DNORA_{source.name}_PATH"
                )
            if default_folder is not None:
                folder = get_url(folder, self._default_folders.get(source))
                if tile_name is not None:
                    folder = re.sub("#TILENAME", tile_name, folder)
                if tile is not None:
                    folder = re.sub("#TILE", tile, folder)

                return folder
        if strict:
            raise ValueError(f"No folder is defined for source {source}!")
        return ""

    def _filename(
        self,
        filename: str,
        source: DataSource,
        tile: Optional[str] = None,
        tile_name: Optional[str] = None,
        strict: bool = True,
    ) -> str:
        if filename is None:
            filename = self._default_filenames.get(source)
        if filename is None:
            filename = self._default_filename
        if filename is None:
            if strict:
                raise ValueError(f"No filename is defined for source {source}!")
            return ""

        if tile_name is not None:
            filename = re.sub("#TILENAME", tile_name, filename)
        if tile is not None:
            filename = re.sub("#TILE", tile, filename)

        return filename


class PointDataReader(DataReader):
    def post_processing(self):
        """Class to use for post processing of data"""
        return None

    @abstractmethod
    def get_coordinates(self, grid, start_time, source, folder, **kwargs):
        """Return a list of ALL the available coordinated in the source.

        These are needed fo the PointPicker object to choose the relevant
        point to actually read in.

        The variables needed to be returned are:

        - coord_dict: Dictionary containing at least lon/lat or x/y.
            E.g. coord_dict = {'lon': np_array_of_longitude, 'lat': np_array_of_latitude}
        """
        pass


class SpectralDataReader(PointDataReader):

    def set_convention(self, convention: Union[SpectralConvention, str]) -> None:
        if isinstance(convention, str):
            self._convention = SpectralConvention[convention.upper()]
        else:
            self._convention = convention

    def convention(self) -> SpectralConvention:
        """Return the convention of the spectra returned to the object.

        The conventions to choose from are predetermined:

        OCEAN:    Oceanic convention
                    Directional vector monotonically increasing.
                    Direction to. North = 0, East = 90.

        MET:      Meteorological convention
                    Directional vector monotonically increasing.
                    Direction from. North = 0, East = 90.

        MATH:     Mathematical convention
                    Direction from. North = 0, East = 90.
                    Direction to. North = 90, East = 0.

        MATHVEC:  Mathematical convention in vector
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        WW3:      WAVEWATCH III output convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 0, East = 90.
        """
        return self._convention


ReaderFunction = Union[DataReader, PointDataReader, SpectralDataReader]
