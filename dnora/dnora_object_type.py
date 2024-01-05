from enum import Enum

from .grid import Grid, TriGrid
from .spectra import Spectra
from .spectra1d import Spectra1D
from .wind import Wind
from .current import Current
from .waterlevel import WaterLevel
from .ice import Ice
from .waveseries import WaveSeries


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
