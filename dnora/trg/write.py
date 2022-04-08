from __future__ import annotations
from abc import ABC, abstractmethod

import numpy as np

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

    def _im_silent(self) -> bool:
        """Return False if you want to be responsible for printing out the
        file names."""
        return True

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

    def _im_silent(self) -> bool:
        """Return False if you want to be responsible for printing out the
        file names."""
        return False

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
            arr_node = np.stack([
                grid.nodes(),
                grid.lon(),
                grid.lat(),
                grid.topo(),
            ]).transpose()
            fmt_node = '%10.0f%22.8f%22.8f%22.8f}'

            f.write('$Nodes\n')
            f.write(f"{len(grid.lon()):12.0f}\n")
            np.savetxt(f, arr_node, fmt=fmt_node)
            f.write('$EndNodes\n')

            # Write Elements
            N_bound = len(grid.boundary_inds())
            arr_bound = np.stack([
                np.arange(1, N_bound+1),
                N_bound*[15],
                N_bound*[2],
                N_bound*[1],
                N_bound*[0],
                grid.boundary_inds()+1
            ]).transpose()
            fmt_bound = '%10.0f%10.0f%10.0f%10.0f%10.0f%10.0f'

            N_ele = len(grid.tri())
            arr_ele = np.stack([
                np.arange(1, N_ele+1)+N_bound,
                N_ele*[2],
                N_ele*[3],
                N_ele*[0],
                np.arange(1, N_ele+1),
                N_ele*[0],
                grid.tri()[:, 0]+1,
                grid.tri()[:, 1]+1,
                grid.tri()[:, 2]+1,
            ]).transpose()
            fmt_ele = '%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f'

            f.write('$Elements\n')
            f.write(f"{len(grid.tri())+len(grid.boundary_inds()):12.0f}\n")
            np.savetxt(f, arr_bound, fmt=fmt_bound)
            np.savetxt(f, arr_ele, fmt=fmt_ele)
            f.write('$EndElements\n')

        return output_file, folder
