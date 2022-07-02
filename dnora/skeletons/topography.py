import numpy as np
import xarray as xr
from copy import copy

class Topography:
    def _reset_vars(self):
        self._set(self.topo(empty=True),'topo', coords='grid')
        self._update_sea_mask(np.logical_not(np.isnan(self.topo())).astype(int))
    #
    # def _update_sea_mask(self) -> None:
    #     sea_mask =
    #     self._set(sea_mask,'sea_mask', coords='grid')

    def topo(self, land: float=0., empty: bool=False) -> np.ndarray:
        """Get the topography, with land set to a given value (default land=0.)

        Set empty=True to get an empty topo (even if it doesn't exist)."""
        if empty:
            return np.full(self.size('all'), 999.)

        topo = self._get('topo')

        if topo is None:
            return None

        topo[np.logical_not(self.sea_mask())] = land

        return topo
