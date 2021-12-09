from abc import ABC, abstractmethod
from .trg_mod import TrGrid

class TrGridWriter(ABC):
    """Abstract class for writing the TrGrid-object's data to files to be Used
    by the wave models.
    """

    def _preferred_format(self):
        return 'General'

    @abstractmethod
    def __call__(self, grid: TrGrid) -> Tuple:
        pass

class WW3(TrGridWriter):
    """Writes the grid to WAVEWATCH III unstructured format."""
    def _preferred_format(self):
        return 'WW3msh'

    def __init__(self) -> None:
        return

    def __call__(self, grid: TrGrid, filename: str, infofilename: str, folder: str) -> Tuple:
        output_file = add_suffix(filename, 'bathy')
        output_files.append(output_file)
        output_path = add_folder_to_filename(output_file, folder)
        grid.write_status(filename=infofilename, folder=folder)

        return output_files, folder
