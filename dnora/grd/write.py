from __future__ import annotations # For TYPE_CHECKING
import numpy as np
from copy import copy
from .. import msg
from ..aux import check_if_folder, add_folder_to_filename
from abc import ABC, abstractmethod

from ..defaults import dflt_grd

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .grd_mod import Grid # Boundary object

class GridWriter(ABC):
    """Abstract class for writing the Grid-object's data to files to be Used
    by the wave models.
    """
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, grid: Grid):
        pass

class WW3(GridWriter):
    def __init__(self, folder: str = '', matrix = False, filestring: str=dflt_grd['fs']['WW3']) -> None:
        self.folder = copy(folder)
        self.matrix = matrix
        self.filestring=copy(filestring)
        return

    def __call__(self, grid: Grid) -> None:
        #msg.header(f"Writing grid to WW3 format to folder: {self.folder}.")
        msg.header(f'{type(self).__name__}: writing grid topography from {grid.name()}')

        mask_out = np.zeros(grid.topo().shape)
        mask_out[grid.land_sea_mask()] = 1
        if grid.boundary_mask().size > 0:
            msg.info(f'Setting {sum(sum(np.logical_and(grid.boundary_mask(), grid.land_sea_mask()))):d} boundary points in grid...')
            mask_out[np.logical_and(grid.boundary_mask(), grid.land_sea_mask())] = 2

        existed = check_if_folder(folder=self.folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {self.folder}")

        output_files = []
        if self.matrix:
            #output_file_bathy = create_filename_obj(filestring='mat_$Grid_bathy.txt', objects=[grid])
            output_file = grid.filename(filestring=self.filestring, prefix='mat', suffix='bathy', extension='txt')
            output_files.append(output_file)
            output_path = add_folder_to_filename(output_file, self.folder)

            msg.to_file(output_path)
            np.savetxt(output_path, grid.topo(), delimiter=',',fmt='%1.6f')

            output_file = grid.filename(filestring=self.filestring, prefix='mat', suffix='mapsta', extension='txt')
            output_files.append(output_file)
            output_path = add_folder_to_filename(output_file, self.folder)

            msg.to_file(output_path)
            np.savetxt(output_path, mask_out, delimiter=',',fmt='%1.0f')
            # output_file_mapsta = create_filename_obj(filestring='mat_$Grid_mapsta.txt', objects=[grid])
            # output_path = add_folder_to_filename(output_file_mapsta, self.folder)
            #
            # msg.to_file(output_path)
            # np.savetxt(output_path, mask_out, delimiter=',',fmt='%1.0f')

        else:
            output_file = grid.filename(filestring=self.filestring, suffix='bathy', extension='txt')
            output_files.append(output_file)
            output_path = add_folder_to_filename(output_file, self.folder)

            msg.to_file(output_path)
            np.savetxt(output_path, grid.topo().ravel(), delimiter=',',fmt='%1.6f')

            output_file = grid.filename(filestring=self.filestring, suffix='mapsta', extension='txt')
            output_files.append(output_file)
            output_path = add_folder_to_filename(output_file, self.folder)

            msg.to_file(output_path)
            np.savetxt(output_path, mask_out.ravel(), delimiter=',',fmt='%1.0f')

        grid.write_status(folder=self.folder)

        return output_files, self.folder

class SWAN(GridWriter):
    def __init__(self, folder: str = '', filestring: str=dflt_grd['fs']['SWAN']) -> None:
        self.folder = copy(folder)
        self.filestring = copy(filestring)
        return

    def __call__(self, grid: Grid):
        msg.header(f'{type(self).__name__}: writing grid topography from {grid.name()}')

        #msg.header(f"Writing grid to SWAN format to folder: {self.folder}.")
        mask_out = np.ones(grid.topo().shape)
        mask_out[grid.land_sea_mask()] = 0
        if grid.boundary_mask().size > 0:
            msg.info(f'Setting {sum(sum(np.logical_and(grid.boundary_mask(), grid.land_sea_mask()))):d} boundary points in grid...')
            mask_out[np.logical_and(grid.boundary_mask(), grid.land_sea_mask())] = 2

        existed = check_if_folder(folder=self.folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {self.folder}")

        #output_file = create_filename_obj(filestring='$Grid_SWAN.bot', objects=[grid])
        output_file = grid.filename(filestring=self.filestring, extension='bot')
        output_path = add_folder_to_filename(output_file, self.folder)

        msg.to_file(output_path)
        np.savetxt(output_path, grid.topo(), delimiter='\t',fmt='%1.2f')
        grid.write_status(folder=self.folder)

        return output_file, self.folder

class SWASH(SWAN):
    def __init__(self, folder: str = '', filestring: str=dflt_grd['fs']['SWASH']) -> None:
        self.folder = copy(folder)
        self.filestring = copy(filestring)
        return

    def __call__(self, grid: Grid) -> None:
        output_file, output_folder = super().__call__(grid)
        return output_file, output_folder
