from enum import Enum, auto


class DataSource(Enum):
    LOCAL = auto()
    INTERNAL = auto()
    REMOTE = auto()
    UNDEFINED = auto()
