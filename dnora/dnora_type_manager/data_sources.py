from enum import Enum, auto


class DataSource(Enum):
    LOCAL = auto()
    INTERNAL = auto()
    IMMUTABLE = auto()
    REMOTE = auto()
    UNDEFINED = auto()


def data_source_from_string(data_str: str | DataSource) -> DataSource:
    if isinstance(data_str, DataSource):
        return data_str
    return DataSource[data_str.upper()]
