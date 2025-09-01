from enum import Enum, auto
from typing import Union

class DataSource(Enum):
    LOCAL = auto()
    INTERNAL = auto()
    IMMUTABLE = auto()
    REMOTE = auto()
    CREATION = auto()
    UNDEFINED = auto()


def data_source_from_string(data_str: Union[str,DataSource]) -> DataSource:
    if isinstance(data_str, DataSource):
        return data_str
    return DataSource[data_str.upper()]
