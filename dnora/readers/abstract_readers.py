from abc import ABC, abstractmethod

# Import objects
from ..grid.grid import Grid
from ..data_sources import DataSource
from ..spectral_conventions import SpectralConvention
from ..dnora_object_type import DnoraDataType


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
        **kwargs
    ):
        """Reads in the forcing witih grid and between start_time and end_time.

        The variables needed to be returned are:

        coord_dict: Dictionary containing at least lon/lat or x/y. Can also contain e.g. time
        data_dict:  Dictionary containing the data variables. Key-names depend on object type
        meta_dict: dict{key, value} will be set as attributes of the xr.Dataset
        """

        pass

    def name(self) -> str:
        return type(self).__name__

    def default_data_source(self) -> DataSource:
        return DataSource.UNDEFINED


class PointDataReader(DataReader):
    @abstractmethod
    def get_coordinates(self, grid, start_time, source, folder):
        """Return a list of all the available coordinated in the source.

        These are needed fo the PointPicker object to choose the relevant
        point to actually read in.

        The variables needed to be returned are:

        coord_dict: Dictionary containing at least lon/lat or x/y. Can also contain e.g. time
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
