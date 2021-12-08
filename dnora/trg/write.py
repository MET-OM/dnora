from abc import ABC, abstractmethod

class TrGridWriter(ABC):
    """Abstract class for writing the TrGrid-object's data to files to be Used
    by the wave models.
    """

    def _preferred_format(self):
        return 'General'

    @abstractmethod
    def __call__(self, grid: Grid) -> Tuple:
        pass
