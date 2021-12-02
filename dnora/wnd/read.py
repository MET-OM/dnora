from abc import ABC,  abstractmethod
# Import objects
from ..grd.grd_mod import Grid

class ForcingReader(ABC):
    """Reads forcing data from some source and provide it to the object.

    The area is defined from the Grid object that is passed.
    """

    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, grid: Grid, start_time: str, end_time: str, expansion_factor: float):
        pass
