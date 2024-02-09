from dnora.dnora_types import DnoraDataType
from typing import Union

from dnora.grid import Grid, TriGrid
from dnora.wind import Wind
from dnora.spectra import Spectra
from dnora.spectra1d import Spectra1D
from dnora.waveseries import WaveSeries
from dnora.waterlevel import WaterLevel
from dnora.current import Current
from dnora.ice import Ice

DnoraObject = Union[
    Grid, TriGrid, Wind, Spectra, Spectra1D, WaveSeries, WaterLevel, Current, Ice
]

dnora_objects = {
    DnoraDataType.GRID: Grid,
    DnoraDataType.TRIGRID: TriGrid,
    DnoraDataType.WIND: Wind,
    DnoraDataType.SPECTRA: Spectra,
    DnoraDataType.SPECTRA1D: Spectra1D,
    DnoraDataType.WAVESERIES: WaveSeries,
    DnoraDataType.WATERLEVEL: WaterLevel,
    DnoraDataType.CURRENT: Current,
    DnoraDataType.ICE: Ice,
}
