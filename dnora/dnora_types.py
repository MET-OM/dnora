from enum import Enum

from dnora.grid import Grid, TriGrid
from dnora.wind import Wind
from dnora.spectra import Spectra
from dnora.spectra1d import Spectra1D
from dnora.waveseries import WaveSeries
from dnora.waterlevel import WaterLevel
from dnora.current import Current
from dnora.ice import Ice

from dnora.wind.read import WindReader
from dnora.spectra.read import SpectraReader
from dnora.spectra1d.read import Spectra1DReader
from dnora.waveseries.read import WaveSeriesReader
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
    SpectraReader,
    Spectra1DReader,
    WaveSeriesReader,
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
    SPECTRA = Spectra
    SPECTRA1D = Spectra1D
    WIND = Wind
    CURRENT = Current
    WATERLEVEL = WaterLevel
    ICE = Ice
    WAVESERIES = WaveSeries


def object_type_from_string(obj_str: str) -> DnoraDataType:
    if isinstance(obj_str, DnoraDataType):
        return obj_str
    return DnoraDataType[obj_str.upper()]
