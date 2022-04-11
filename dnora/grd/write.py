from __future__ import annotations # For TYPE_CHECKING

import numpy as np
from copy import copy
from abc import ABC, abstractmethod
import utm

# Import objects
from typing import TYPE_CHECKING, Tuple
if TYPE_CHECKING:
    from .grd_mod import Grid

# Import default values and auxiliry functions
from ..defaults import dflt_grd, list_of_placeholders
from .. import msg
from ..aux import add_folder_to_filename, add_prefix, add_suffix, clean_filename

class GridWriter(ABC):
    """Abstract class for writing the Grid-object's data to files to be Used
    by the wave models.
    """

    def _preferred_extension(self):
        return 'txt'

    def _preferred_format(self):
        return 'General'

    def _im_silent(self) -> bool:
        """Return False if you want to be responsible for printing out the
        file names."""
        return True

    @abstractmethod
    def __call__(self, grid: Grid, filename: str, infofilename: str, folder: str) -> Tuple:
        pass

class BoundaryPoints(GridWriter):
    """Writes boundary points from unsutructured grid."""
    def __init__(self, include_index = False):
        self.include_index = copy(include_index)
        return

    def __call__(self, grid: Grid, filename: str, infofilename: str, folder: str) -> Tuple:
        output_file = clean_filename(filename, list_of_placeholders)
        output_path = add_folder_to_filename(output_file, folder)


        with open(output_path,'w') as f:
            if self.include_index and hasattr(grid, 'boundary_inds'):
                for n in range(len(grid.boundary_points())):
                    # Going from python 0-indexing to node 1-indexing
                    f.write(f'{grid.boundary_inds()[n]+1:.0f} {grid.boundary_points()[n,0]:13.10f} {grid.boundary_points()[n,1]:13.10f}\n')
            else:
                for n in range(len(grid.boundary_points())):
                    f.write(f'{grid.boundary_points()[n,0]:13.10f} {grid.boundary_points()[n,1]:13.10f}\n')
        print(output_path)
        return output_file, folder

class WW3(GridWriter):
    """Writes the grid to WAVEWATCH III format."""
    def _preferred_format(self):
        return 'WW3'

    def __init__(self, matrix=False) -> None:
        self.matrix = matrix
        return

    def __call__(self, grid: Grid, filename: str, infofilename: str, folder: str) -> Tuple:

        mask_out = np.zeros(grid.topo().shape)
        mask_out[grid.land_sea_mask()] = 1
        if grid.boundary_mask().size > 0:
            msg.info(f'Setting {sum(sum(np.logical_and(grid.boundary_mask(), grid.land_sea_mask()))):d} boundary points in grid...')
            mask_out[np.logical_and(grid.boundary_mask(), grid.land_sea_mask())] = 2

        output_files = []
        if self.matrix:
            output_file = add_prefix(filename, 'mat')
            output_file = add_suffix(output_file, 'bathy')
            output_files.append(output_file)
            output_path = add_folder_to_filename(output_file, folder)

            #msg.to_file(output_path)
            np.savetxt(output_path, grid.topo(), delimiter=',',fmt='%1.6f')

            output_file = add_prefix(filename, 'mat')
            output_file = add_suffix(output_file, 'mapsta')
            output_files.append(output_file)
            output_path = add_folder_to_filename(output_file, folder)

            #msg.to_file(output_path)
            np.savetxt(output_path, mask_out, delimiter=',',fmt='%1.0f')

        else:
            output_file = add_suffix(filename, 'bathy')
            output_files.append(output_file)
            output_path = add_folder_to_filename(output_file, folder)

            #msg.to_file(output_path)
            np.savetxt(output_path, grid.topo().ravel(), delimiter=',',fmt='%1.6f')

            output_file = add_suffix(filename, 'mapsta')
            output_files.append(output_file)
            output_path = add_folder_to_filename(output_file, folder)

            #msg.to_file(output_path)
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

    def _preferred_extension(self):
        return 'bot'

    def __call__(self, grid: Grid, filename: str, infofilename: str, folder: str) -> None:

        mask_out = np.ones(grid.topo().shape)
        mask_out[grid.land_sea_mask()] = 0
        if grid.boundary_mask().size > 0:
            msg.info(f'Setting {sum(sum(np.logical_and(grid.boundary_mask(), grid.land_sea_mask()))):d} boundary points in grid...')
            mask_out[np.logical_and(grid.boundary_mask(), grid.land_sea_mask())] = 2

        #output_file = grid.filename(filestring=filestring, extension='bot')
        output_path = add_folder_to_filename(filename, folder)

        #msg.to_file(output_path)
        np.savetxt(output_path, grid.topo(), delimiter='\t',fmt='%1.2f')
        grid.write_status(filename=infofilename, folder=folder)

        return filename, folder

class Xyz(GridWriter):
    """Writes the grid to Xyz-format."""
    def __init__(self, use_raw=False, utm=False):
        self._use_raw = use_raw
        self._utm = utm
        pass

    def __call__(self, grid: Grid, filename: str, infofilename: str, folder: str) -> None:

        #output_file = grid.filename(filestring=filestring, extension='bot')
        output_path = add_folder_to_filename(filename, folder)

        if self._use_raw:
            z = grid.raw_topo();
            x = grid.raw_lon()
            y = grid.raw_lat()
        else:
            z = grid.topo();
            x = grid.lon()
            y = grid.lat()



        with open(output_path, 'w') as f:
            f.write(f'{filename}\n')
            for nx, lon in enumerate(x):
                for ny, lat in enumerate(y):
                    if self._utm:
                        [lon_out, lat_out, utm_zone, utm_letter] = utm.from_latlon(lat, lon, 33, 'W')
                        fmt='.4f'
                    else:
                        lon_out=lon
                        lat_out=lat
                        fmt='.9f'
                    if ~np.isnan(z[ny,nx]):
                        f.write(f'{lon_out:{fmt}},{lat_out:{fmt}},{z[ny,nx]:.1f}\n')
        #np.savetxt(output_path, grid.topo(), delimiter='\t',fmt='%1.2f')
        #grid.write_status(filename=infofilename, folder=folder)

        return filename, folder

class REEF3D(GridWriter):
    """Writes the grid to Xyz-format in relative Cartesian grid."""
    def __init__(self, use_raw=False):
        self._use_raw = use_raw
        pass

    def _preferred_extension(self):
        return 'dat'

    def _preferred_format(self):
        return 'REEF3D'

    def __call__(self, grid: Grid, filename: str, infofilename: str, folder: str) -> None:

        output_path = add_folder_to_filename(filename, folder)

        if self._use_raw:
            z = grid.raw_topo();
            [x, y, utm_zone, utm_letter] = utm.from_latlon(grid.raw_lat(), grid.raw_lon())
            x = x-min(x) # metres from corner
            y = y-min(y)
        else:
            z = grid.topo();
            x, y = np.meshgrid(np.linspace(0,grid.nx()*grid.dx(),grid.nx()), np.linspace(0,grid.ny()*grid.dy(),grid.ny()))
            x = x.ravel()
            y = y.ravel()
            z = z.ravel()

        print('Max depth(m):',np.nanmax(z))
        z = -1*z +np.nanmax(z)# set at zero the maximum depth.
        z[np.isnan(z)] = np.nanmax(z) + 14 # + 14 for land points

        fmt='.5f'
        with open(output_path, 'w') as f:
            for i, __ in enumerate(x):
                f.write(f'{x[i]:{fmt}} {y[i]:{fmt}} {z[i]:.1f}\n')

        #np.savetxt(output_path, grid.topo(), delimiter='\t',fmt='%1.2f')
        #grid.write_status(filename=infofilename, folder=folder)

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
