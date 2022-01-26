from __future__ import annotations
from abc import ABC, abstractmethod
from ..aux import add_suffix, add_folder_to_filename, clean_filename
from .. import msg
from ..defaults import list_of_placeholders
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .trg_mod import Grid
class TrGridWriter(ABC):
    """Abstract class for writing the TrGrid-object's data to files to be Used
    by the wave models.
    """

    def _preferred_format(self):
        return 'General'

    def _preferred_extension(self):
        return 'txt'

    @abstractmethod
    def __call__(self, grid: Grid) -> Tuple:
        pass

    def _im_silent(self) -> bool:
        """Return False if you want to be responsible for printing out the
        file names."""
        return True



class WW3(TrGridWriter):
    """Writes the grid to WAVEWATCH III unstructured format."""
    def _preferred_format(self):
        return 'WW3'

    def _preferred_extension(self):
        return 'msh'

    def __init__(self) -> None:
        return

    def __call__(self, grid: Grid, filename: str, infofilename: str, folder: str) -> Tuple:

        output_file = add_suffix(filename, 'bathy')
        output_path = add_folder_to_filename(output_file, folder)
        output_path = clean_filename(output_path, list_of_placeholders)

        grid.write_status(filename=infofilename, folder=folder)

        msg.to_file(output_path)

        with open(output_path,'w') as f:
            # Write header
            f.write('$MeshFormat\n')
            f.write('2 0 8\n')
            f.write('$EndMeshFormat\n')

            # Write Nodes
            f.write('$Nodes\n')
            f.write(f"{len(grid.lon()):12.0f}\n")
            #inds = list(range(1,len(lon)+1))
            for n in range(len(grid.nodes())):
                f.write(f"{grid.nodes()[n]:10.0f}{grid.lon()[n]:22.8f}{grid.lat()[n]:22.8f}{abs(grid.topo()[n]):22.8f}\n")
            f.write('$EndNodes\n')

            # Write Elements
            f.write('$Elements\n')
            f.write(f"{len(grid.tri())+len(grid.boundary_inds()):12.0f}\n")
            # Write boundary points
            ct = 0
            for n in range(len(grid.boundary_inds())):
                ct += 1
                f.write(f"{ct:10.0f}{15:10.0f}{2:10.0f}{1:10.0f}{0:10.0f}{grid.boundary_inds()[n]+1:10.0f}\n")
            for n in range(len(grid.tri())):
                ct += 1
                #161       2       3       0       1       0       1       2       3
                f.write(f"{ct:8.0f}{2:8.0f}{3:8.0f}{0:8.0f}{n+1:8.0f}{0:8.0f}{grid.tri()[n,0]+1:8.0f}{grid.tri()[n,1]+1:8.0f}{grid.tri()[n,2]+1:8.0f}\n")
            f.write('$EndElements\n')

        return output_file, folder
