from __future__ import annotations
from abc import ABC, abstractmethod

# Import objects
from grid import Grid
from data_sources import DataSource
from spectral_conventions import SpectralConvention

from metaparameter import metaparameter

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..dnora_types import DnoraDataType


class DataReader(ABC):
    """Reads forcing data from some source and provide it to the object.

    The area is defined from the Grid object that is passed.
    """

    @abstractmethod
    def __call__(
        self,
        obj_type: DnoraDataType,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        folder: str,
        **kwargs,
    ):
        """Reads in the forcing witih grid and between start_time and end_time.

        The variables needed to be returned are:

        coord_dict: Dictionary containing at least lon/lat or x/y. Can also contain e.g. time
        data_dict:  Dictionary containing the data variables. Key-names depend on object type
        meta_dict: dict{key, value} will be set as attributes of the xr.Dataset. Can be empty.
        metaparameter_dict: dict{key, values} couples the data variables to a metaparametr. Can be empty.
        """

        pass

    def name(self) -> str:
        return type(self).__name__

    def default_data_source(self) -> DataSource:
        return DataSource.UNDEFINED

    def create_metaparameter_dict(self, parameter_strings: list[str]):
        metaparameter_dict = {}
        for param in parameter_strings:
            val = metaparameter.get(param)
            if val is not None:
                metaparameter_dict[param] = val
        return metaparameter_dict


class PointDataReader(DataReader):
    @abstractmethod
    def get_coordinates(self, grid, start_time, source, folder):
        """Return a list of all the available coordinated in the source.

        These are needed fo the PointPicker object to choose the relevant
        point to actually read in.

        The variables needed to be returned are:

        lon_all, lat_all, x_all, y_all
        """
        pass


class SpectralDataReader(PointDataReader):
    @abstractmethod
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
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        MATHVEC:  Mathematical convention in vector
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        WW3:      WAVEWATCH III output convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 0, East = 90.
        """
        pass

    def post_processing(self):
        return None
