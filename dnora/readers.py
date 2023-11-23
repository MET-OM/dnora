from abc import ABC, abstractmethod

# Import objects
from .grd.grd_mod import Grid
from .data_sources import DataSource


class DataReader(ABC):
    """Reads forcing data from some source and provide it to the object.

    The area is defined from the Grid object that is passed.
    """

    @abstractmethod
    def __call__(
        self,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        folder: str,
        **kwargs
    ):
        """Reads in the forcing witih grid and between start_time and end_time.

        The variables needed to be returned are:

        time:   Time stamps as numpy.datetime64 array
        u:      west-to-east velocity [lon, lat, time] as numpy array
        v:      south-to-north velocity [lon, lat, time] as numpy array
        lon:    Longitude vector as numpy array (None if Cartesian)
        lat:    Latitude vector as numpy array (None if Cartesian)
        x:      Longitude vector as numpy array (None if Spherical)
        y:      Latitude vector as numpy array (None if Spherical)
        metadata: dict{key, value} will be set as attributes of the xr.Dataset
        """

        return coord_dict, data_dict, meta_dict

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

        Provide the result as for equally long nump arrays (None for missing values).
        """
        return coord_dict


class SpectralDataReader(PointDataReader):
    @abstractmethod
    def convention(self):
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
        return convention

    def post_processing(self):
        return None
