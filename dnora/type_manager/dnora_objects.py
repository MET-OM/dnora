from .dnora_types import DnoraDataType
from typing import Union

from dnora.grid import Grid, TriGrid
from dnora.wind import Wind
from dnora.spectra import Spectra
from dnora.spectra1d import Spectra1D
from dnora.waveseries import WaveSeries
from dnora.waterlevel import WaterLevel
from dnora.current import Current
from dnora.ice import Ice
from dnora.wavegrid import WaveGrid

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
    WaveGrid,
]

dnora_objects = {
    DnoraDataType.GRID: Grid,
    DnoraDataType.TRIGRID: TriGrid,
    DnoraDataType.WIND: Wind,
    DnoraDataType.SPECTRA: Spectra,
    DnoraDataType.SPECTRA1D: Spectra1D,
    DnoraDataType.WAVESERIES: WaveSeries,
    DnoraDataType.WAVEGRID: WaveGrid,
    DnoraDataType.WATERLEVEL: WaterLevel,
    DnoraDataType.CURRENT: Current,
    DnoraDataType.ICE: Ice,
}
