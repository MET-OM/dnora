import numpy as np
import xarray as xr
from copy import copy
from .. import aux_funcs

from ..grd.boundary import MaskSetter
from ..grd.process import GridProcessor
from ..grd.mesh import Mesher, Interpolate
from .. import msg
def topography_methods(c):
    def set_mask(self, mask_setter: MaskSetter, mask_type: str=None) -> None:
        """Set a mask that represents e.g. Boundary points or spectral output
        point.

        NB! Points can overlap with land!
        """

        if mask_type is None:
            mask_type = mask_setter._mask_type()

        if mask_type is None:
            msg.advice(f"Either provide mask_type variable or provide a MaskSetter with a _mask_type method.")

        msg.header(mask_setter, f"Setting {mask_type} points...")
        print(mask_setter)

        mask = mask_setter(self.sea_mask())
        self._update_mask(mask_type, mask)

    def mesh_grid(self, mesher: Mesher=Interpolate(method = 'nearest')) -> None:
        """Meshes the raw data down to the grid definitions."""

        if not hasattr(self, 'raw'):
            msg.warning('Import topography using .import_topo() before meshing!')
            return

        msg.header(mesher, "Meshing grid bathymetry...")
        print(mesher)

        if self.is_gridded():
            xQ, yQ = np.meshgrid(self.native_x(), self.native_y())
        else:
            xQ, yQ = self.native_xy()

        if self.is_cartesian():
            x, y = self.raw().xy()
        else:
            x, y = self.raw().lonlat()

        topo = mesher(self.raw().topo().ravel(), x, y, xQ, yQ)

        self._update_datavar('topo', topo)
        self._update_sea_mask()
        print(self)

    def process_topo(self, grid_processor: GridProcessor=None) -> None:
        """Processes the raw bathymetrical data, e.g. with a filter."""
        if grid_processor is None:
            return

        msg.header(grid_processor, "Processing topography...")

        print(grid_processor)
        if self.raw().is_gridded():
            topo = grid_processor.grid(self.raw().topo(), self.raw().lon(), self.raw().lat(), self.raw().sea_mask(), self.raw().boundary_mask())
            if topo is None:
                msg.warning('Filtering of gridded topography is not implemented in this GridProcessor.')
                return
        else:
            topo = grid_processor.topo(self.raw().topo(), self.raw().lon(), self.raw().lat(), self.raw().sea_mask())
            if topo is None:
                msg.warning('Filtering of unstructured topography is not implemented in this GridProcessor.')
                return

        self.raw()._update_datavar('topo', topo)
        self.raw()._update_sea_mask()

    def process_grid(self, grid_processor: GridProcessor=None) -> None:
        """Processes the gridded bathymetrical data, e.g. with a filter."""
        if grid_processor is None:
            return

        msg.header(grid_processor, "Processing meshed grid...")
        print(grid_processor)
        if self.is_gridded():
            topo = grid_processor.grid(self.topo(), self.lon(), self.lat(), self.sea_mask(), self.boundary_mask())
            if topo is None:
                msg.warning('Filtering of gridded topography is not implemented in this GridProcessor.')
                return
        else:
            topo = grid_processor.topo(self.topo(), self.lon(), self.lat(), self.sea_mask())
            if topo is None:
                msg.warning('Filtering of unstructured topography is not implemented in this GridProcessor.')
                return

        self._update_datavar('topo', topo)
        self._update_sea_mask()

    def update_sea_mask(self):
        self._update_mask('sea', (self.ds().topo.values>0).astype(int))

    def topo(self, land: float=0., empty=False, **kwargs) -> np.ndarray:
        """Returns an array containing the meshed topography of the grid."""
        topo = self._topo(empty=empty, **kwargs)
        if topo is None or empty:
            return topo
        topo[np.logical_not(self.sea_mask(**kwargs))] = land
        return topo

    c.set_mask = set_mask
    c.mesh_grid = mesh_grid
    c.process_topo = process_topo
    c.process_grid = process_grid
    c._update_sea_mask = update_sea_mask
    c.topo = topo
    return c
