import numpy as np
from copy import copy
from .. import msg
from ..aux import check_if_folder, create_filename_obj, add_folder_to_filename

from .grd_mod import TopoWriter # Abstract class
from .grd_mod import Grid # Grid object

class WW3(TopoWriter):
    def __init__(self, folder: str = '', matrix = False) -> None:
        self.folder = copy(folder)
        self.matrix = matrix
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

        if self.matrix:
            output_file_bathy = create_filename_obj(filestring='mat_$Grid_bathy.txt', objects=[grid])
            output_path = add_folder_to_filename(output_file_bathy, self.folder)

            msg.to_file(output_path)
            np.savetxt(output_path, -1*grid.topo(), delimiter=',',fmt='%1.6f')

            output_file_mapsta = create_filename_obj(filestring='mat_$Grid_mapsta.txt', objects=[grid])
            output_path = add_folder_to_filename(output_file_mapsta, self.folder)

            msg.to_file(output_path)
            np.savetxt(output_path, mask_out, delimiter=',',fmt='%1.0f')

        else:
            output_file_bathy = create_filename_obj(filestring='$Grid_bathy.txt', objects=[grid])
            output_path = add_folder_to_filename(output_file_bathy, self.folder)

            msg.to_file(output_path)
            np.savetxt(output_path, -1*grid.topo().ravel(), delimiter=',',fmt='%1.6f')

            output_file_mapsta =  create_filename_obj(filestring='$Grid_mapsta.txt', objects=[grid])
            output_path = add_folder_to_filename(output_file_mapsta, self.folder)

            msg.to_file(output_path)
            np.savetxt(output_path, mask_out.ravel(), delimiter=',',fmt='%1.0f')

        grid.write_status(folder=self.folder)

        # This is set as info in case an input file needs to be generated
        grid._written_as = [output_file_bathy, output_file_mapsta]
        grid._written_to = self.folder

        return

class SWAN(TopoWriter):
    def __init__(self, folder: str = '') -> None:
        self.folder = copy(folder)
        return

    def __call__(self, grid: Grid) -> None:
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

        output_file = create_filename_obj(filestring='$Grid_SWAN.bot', objects=[grid])
        output_path = add_folder_to_filename(output_file, self.folder)

        msg.to_file(output_path)
        np.savetxt(output_path, grid.topo(), delimiter='\t',fmt='%1.0f')
        grid.write_status(folder=self.folder)

        # This is set as info in case an input file needs to be generated
        grid._written_as = output_file
        grid._written_to = self.folder
        return
