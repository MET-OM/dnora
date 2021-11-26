from __future__ import annotations # For TYPE_CHECKING
import numpy as np
from copy import copy
from .. import msg
from ..aux import check_if_folder, add_folder_to_filename, add_prefix, add_suffix
from abc import ABC, abstractmethod

from ..defaults import dflt_grd

from typing import TYPE_CHECKING, Tuple
if TYPE_CHECKING:
    from .grd_mod import Grid # Boundary object

class GridWriter(ABC):
    """Abstract class for writing the Grid-object's data to files to be Used
    by the wave models.
    """

    def _preferred_format(self):
        return 'General'

    @abstractmethod
    def __call__(self, grid: Grid) -> Tuple:
        pass

class WW3(GridWriter):
    """Writes the grid to WAVEWATCH III format."""
    def _preferred_format(self):
        return 'WW3'

    def __init__(self, matrix=False, ) -> None:
        self.matrix = matrix
        return

    def __call__(self, grid: Grid, filename: str, infofilename: str, folder: str) -> Tuple:
        msg.header(f'{type(self).__name__}: writing grid topography from {grid.name()}')

        mask_out = np.zeros(grid.topo().shape)
        mask_out[grid.land_sea_mask()] = 1
        if grid.boundary_mask().size > 0:
            msg.info(f'Setting {sum(sum(np.logical_and(grid.boundary_mask(), grid.land_sea_mask()))):d} boundary points in grid...')
            mask_out[np.logical_and(grid.boundary_mask(), grid.land_sea_mask())] = 2

        existed = check_if_folder(folder=folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {folder}")

        output_files = []
        if self.matrix:
            output_file = add_prefix(filename, 'mat')
            output_file = add_suffix(output_file, 'bathy')
            output_files.append(output_file)
            output_path = add_folder_to_filename(output_file, folder)

            msg.to_file(output_path)
            np.savetxt(output_path, grid.topo(), delimiter=',',fmt='%1.6f')

            output_file = add_prefix(filename, 'mat')
            output_file = add_suffix(output_file, 'mapsta')
            output_files.append(output_file)
            output_path = add_folder_to_filename(output_file, folder)

            msg.to_file(output_path)
            np.savetxt(output_path, mask_out, delimiter=',',fmt='%1.0f')

        else:
            output_file = add_suffix(filename, 'bathy')
            output_files.append(output_file)
            output_path = add_folder_to_filename(output_file, folder)

            msg.to_file(output_path)
            np.savetxt(output_path, grid.topo().ravel(), delimiter=',',fmt='%1.6f')

            output_file = add_suffix(filename, 'mapsta')
            output_files.append(output_file)
            output_path = add_folder_to_filename(output_file, folder)

            msg.to_file(output_path)
            np.savetxt(output_path, mask_out.ravel(), delimiter=',',fmt='%1.0f')

        grid.write_status(filename=infofilename, folder=folder)

        return output_files, folder

class SWAN(GridWriter):
    """Writes the grid to SWAN format."""
    def __init__(self, out_format='SWAN'):
        self.out_format = out_format
        return

    def _preferred_format(self):
        return self.out_format

    def __call__(self, grid: Grid, filename: str, infofilename: str, folder: str) -> None:
        msg.header(f'{type(self).__name__}: writing grid topography from {grid.name()}')

        #msg.header(f"Writing grid to SWAN format to folder: {self.folder}.")
        mask_out = np.ones(grid.topo().shape)
        mask_out[grid.land_sea_mask()] = 0
        if grid.boundary_mask().size > 0:
            msg.info(f'Setting {sum(sum(np.logical_and(grid.boundary_mask(), grid.land_sea_mask()))):d} boundary points in grid...')
            mask_out[np.logical_and(grid.boundary_mask(), grid.land_sea_mask())] = 2

        existed = check_if_folder(folder=folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {folder}")

        #output_file = grid.filename(filestring=filestring, extension='bot')
        output_path = add_folder_to_filename(filename, folder)

        msg.to_file(output_path)
        np.savetxt(output_path, grid.topo(), delimiter='\t',fmt='%1.2f')
        grid.write_status(filename=infofilename, folder=folder)

        return filename, folder

# class SWASH(SWAN):
#     """Writes the grid to SWASH format.
#
#     NB! Format is same as SWAN, and a separate class is defined only to get
#     the file name right.
#     """
#     def _preferred_format(self):
#         return 'SWASH'
#
#     def __call__(self, grid: Grid, filename: str, infofilename: str, folder: str) -> Tuple:
#         filename, folder = super().__call__(grid)
#         return filename, folder
