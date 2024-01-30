from enum import Enum, auto

from dnora.grid import Grid, TriGrid
from dnora.wind import Wind
from dnora.spectra import Spectra
from dnora.spectra1d import Spectra1D
from dnora.waveseries import WaveSeries
from dnora.waterlevel import WaterLevel
from dnora.current import Current
from dnora.ice import Ice

from dnora.wind.read import WindReader
from dnora.waterlevel.read import WaterLevelReader
from dnora.current.read import CurrentReader
from dnora.ice.read import IceReader

from dnora.readers import abstract_readers
from typing import Union

ReaderFunction = Union[
    abstract_readers.DataReader,
    abstract_readers.PointDataReader,
    abstract_readers.SpectralDataReader,
    WindReader,
    WaterLevelReader,
    CurrentReader,
    IceReader,
]

DnoraObject = Union[
    Grid,
    TriGrid,
    Wind,
    Spectra,
    Spectra1D,
    WaveSeries,
    WaterLevel,
    Current,
    Ice,
]


class DnoraDataType(Enum):
    GRID = Grid
    TRIGRID = TriGrid
    SPECTRA1D = (
        Spectra1D  # Spectra1D needs to be before Spectra because of flename creation
    )
    SPECTRA = Spectra
    WIND = Wind
    CURRENT = Current
    WATERLEVEL = WaterLevel
    ICE = Ice
    WAVESERIES = WaveSeries


class DnoraFileType(Enum):
    INPUTFILE = auto()


def object_type_from_string(obj_str: str) -> DnoraDataType:
    if isinstance(obj_str, DnoraDataType):
        return obj_str
    return DnoraDataType[obj_str.upper()]
