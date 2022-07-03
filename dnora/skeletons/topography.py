import numpy as np
import xarray as xr
from copy import copy
import aux_funcs

from dnora.grd.boundary import BoundarySetter
from dnora.grd.process import GridProcessor
from dnora import msg
def topography_methods(c):
    def set_boundary(self, boundary_setter: BoundarySetter) -> None:
        """Marks the points that should be treated as boundary points in the
        grid.

        The boundary points are stored in a boolean array where True values
        mark a boundary point.

        NB! No check for land points are done, so it is possible that a land
        point is marked as a boundary point. Possibly accounting for this is
        the responsibility of the GridWriter.
        """

        msg.header(boundary_setter, "Setting boundary points...")
        print(boundary_setter)

        boundary_mask = boundary_setter(self.sea_mask())
        self._update_mask('boundary', boundary_mask)

    def process_topo(self, grid_processor: GridProcessor=None) -> None:
        """Processes the raw bathymetrical data, e.g. with a filter."""
        if grid_processor is None:
            return

        msg.header(filt, "Filtering topography...")
        land_sea_mask = self.raw_topo() > 0 # Sea points set to true

        print(grid_processor)
        if self.raw.is_gridded():
            topo = grid_processor.grid(self.raw.topo(), self.raw.lon(), self.raw.lat(), self.raw.sea_mask(), self.raw.boundary_mask())
            if topo is None:
                msg.warning('Filtering of gridded topography is not implemented in this GridProcessor.')
                return
        else:
            topo = grid_processor.topo(self.raw.topo(), self.raw.lon(), self.raw.lat(), self.raw.sea_mask())
            if topo is None:
                msg.warning('Filtering of unstructured topography is not implemented in this GridProcessor.')
                return

        self.raw._update_datavar('topo', topo)

    def process_grid(self, grid_processor: GridProcessor=None) -> None:
        """Processes the gridded bathymetrical data, e.g. with a filter."""
        if grid_processor is None:
            return

        msg.header(filt, "Filtering meshed grid...")
        print(filt)
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

    def update_sea_mask(self):
        self._update_mask('sea', (self.topo()>0).astype(int))

    def topo(self, land: float=0., empty=False) -> np.ndarray:
        """Returns an array containing the meshed topography of the grid."""
        topo = self._topo(empty=empty)
        if topo is None or empty:
            return topo
        topo[np.logical_not(self.sea_mask())] = land
        return topo

    c.set_boundary = set_boundary
    c.process_topo = process_topo
    c.process_grid = process_grid
    c._update_sea_mask = update_sea_mask
    c.topo = topo
    return c
