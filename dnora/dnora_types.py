from enum import Enum

from grid.grid import Grid, TriGrid
from wind.wind import Wind
from spectra.spectra import Spectra
from spectra1d.spectra1d import Spectra1D
from waveseries.waveseries import WaveSeries
from waterlevel.waterlevel import WaterLevel
from current.current import Current
from ice.ice import Ice

from wind.read import WindReader
from spectra.read import SpectraReader
from spectra1d.read import Spectra1DReader
from waveseries.read import WaveSeriesReader
from waterlevel.read import WaterLevelReader
from current.read import CurrentReader
from ice.read import IceReader

from readers import generic_readers
from typing import Union

ReaderFunction = Union[
    generic_readers.DataReader,
    generic_readers.PointDataReader,
    generic_readers.SpectralDataReader,
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
