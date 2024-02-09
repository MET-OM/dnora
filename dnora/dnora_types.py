from enum import Enum, auto


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


class DnoraFileType(Enum):
    GRID = auto()
    SPECTRA1D = auto()
    SPECTRA = auto()
    WIND = auto()
    CURRENT = auto()
    WATERLEVEL = auto()
    ICE = auto()
    WAVESERIES = auto()
    INPUT = auto()


class DataSource(Enum):
    LOCAL = auto()
    INTERNAL = auto()
    REMOTE = auto()
    UNDEFINED = auto()


def data_type_from_string(obj_str: str | DnoraDataType) -> DnoraDataType:
    if isinstance(obj_str, DnoraDataType):
        return obj_str
    return DnoraDataType[obj_str.upper()]


def file_type_from_string(obj_str: str | DnoraFileType) -> DnoraFileType:
    if isinstance(obj_str, DnoraFileType):
        return obj_str
    return DnoraFileType[obj_str.upper()]


def data_source_from_string(data_str: str | DataSource) -> DataSource:
    if isinstance(data_str, DataSource):
        return data_str
    return DataSource[data_str.upper()]
