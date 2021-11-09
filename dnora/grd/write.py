import numpy as np
from copy import copy
from . import msg

from .grd_mod import TopoWriter # Abstract class

class WW3(TopoWriter):
    def __init__(self, folder = ''):
        if (not folder == '') and (not folder[-1] == '/'):
            folder = folder + '/'
        self.folder = copy(folder)

    def __call__(self, grid, matrix = False):
        msg.header('Create files for regular grid')
        mask_out = np.zeros(grid.topo().shape)
        mask_out[grid.land_sea_mask()] = 1
        if grid.boundary_mask().size > 0:
            msg.info(f'Setting {sum(sum(np.logical_and(grid.boundary_mask(), grid.land_sea_mask()))):d} boundary points in grid...')
            mask_out[np.logical_and(grid.boundary_mask(), grid.land_sea_mask())] = 2

        if matrix:
            fn1 = self.folder + 'mat_'+grid.name()+'_bathy.txt'
            msg.to_file(fn1)
            np.savetxt(fn1, -1*grid.topo(), delimiter=',',fmt='%1.6f')

            fn2 = self.folder + 'mat_'+grid.name()+'_mapsta.txt'
            msg.to_file(fn2)
            np.savetxt(fn2, mask_out, delimiter=',',fmt='%1.0f')
        else:
            fn1 = self.folder + grid.name()+'_bathy.txt'
            msg.to_file(fn1)
            np.savetxt(fn1, -1*grid.topo().ravel(), delimiter=',',fmt='%1.6f')

            fn2 = self.folder + grid.name()+'_mapsta.txt'
            msg.to_file(fn2)
            np.savetxt(fn2, mask_out.ravel(), delimiter=',',fmt='%1.0f')

        grid.write_status(folder = self.folder)


class SWAN(TopoWriter):
    def __init__(self):
        pass

    def __call__(self, grid):
        msg.header('Create files for regular grid')
        mask_out = np.ones(grid.topo().shape)
        mask_out[grid.land_sea_mask()] = 0
        if grid.boundary_mask().size > 0:
            msg.info(f'Setting {sum(sum(np.logical_and(grid.boundary_mask(), grid.land_sea_mask()))):d} boundary points in grid...')
            mask_out[np.logical_and(grid.boundary_mask(), grid.land_sea_mask())] = 2

        fn1 = grid.name()+'_SWAN.bot'
        msg.to_file(fn1)
        np.savetxt(fn1, grid.topo(), delimiter='\t',fmt='%1.0f')
        grid.write_status()

