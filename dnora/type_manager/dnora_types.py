from enum import Enum, auto
from typing import Union

class DnoraDataType(Enum):
    GRID = auto()
    TRIGRID = auto()
    # Spectra1D needs to be before Spectra because of filename creation
    SPECTRA1D = auto()
    SPECTRA = auto()
    WIND = auto()
    CURRENT = auto()
    WATERLEVEL = auto()
    ICE = auto()
    WAVESERIES = auto()
    WAVEGRID = auto()
    SPECTRALGRID = auto()


class DnoraFileType(Enum):
    GRID = auto()
    TRIGRID = auto()
    SPECTRA1D = auto()
    SPECTRA = auto()
    WIND = auto()
    CURRENT = auto()
    WATERLEVEL = auto()
    ICE = auto()
    WAVESERIES = auto()
    WAVEGRID = auto()
    INPUT = auto()


def data_type_from_string(obj_str: Union[str, DnoraDataType]) -> DnoraDataType:
    if isinstance(obj_str, str):
        obj_str = DnoraDataType[obj_str.upper()]
    return obj_str


def file_type_from_string(obj_str: Union[str, DnoraFileType]) -> DnoraFileType:
    if isinstance(obj_str, str):
        obj_str = DnoraFileType[obj_str.upper()]
    return obj_str
