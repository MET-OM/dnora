import numpy as np
import xarray as xr
from copy import copy

class Topography:
    def _reset_vars(self):
        self._set_data(self.topo(empty=True),'topo')
        self._set_data(self.boundary_mask(empty=True),'boundary_mask')
        self._update_sea_mask()

    def _update_sea_mask(self) -> None:
        sea_mask = np.logical_not(np.isnan(self.topo())).astype(int)
        self._set_data(sea_mask,'sea_mask')

    def sea_mask(self, boolean: bool=True, empty: bool=False) -> np.ndarray:
        """Returns bool array of the sea mask.

        Set boolean=False to get 0 for land and 1 for sea.
        Set empty=True to get an empty mask (even if it doesn't exist)"""

        data_type = 'bool'*boolean or 'float'

        if empty:
            return np.full(self.core_size(), 1.).astype(data_type)

        mask = self.get('sea_mask')

        if mask is None:
            return None

        return mask.astype(data_type)

    def boundary_mask(self, boolean: bool=True, empty: bool=False) -> np.ndarray:
        """Returns bool array of the boundary mask.

        Set boolean=False to get 1 for boundary points
        Set empty=True to get an empty mask (even if it doesn't exist)."""

        data_type = 'bool'*boolean or 'float'

        if empty:
            return np.full(self.core_size(), 1.).astype(data_type)

        mask = self.get('boundary_mask')

        if mask is None:
            return None

        return mask.astype(data_type)

    def topo(self, land: float=0., empty: bool=False) -> np.ndarray:
        """Get the topography, with land set to a given value (default land=0.)

        Set empty=True to get an empty topo (even if it doesn't exist)."""
        if empty:
            return np.full(self.core_size(), 999.)

        topo = self.get('topo')

        if topo is None:
            return None

        topo[np.logical_not(self.sea_mask())] = land

        return topo

    def boundary_points(self, order_by: str='lat'):
        return self.native_xy(mask=self.boundary_mask(), order_by=order_by)

    def boundary_points_xy(self, order_by: str='y'):
        return self.xy(mask=self.boundary_mask(), order_by=order_by)

    def boundary_points_lonlat(self, order_by: str='lat'):
        return self.lonlat(mask=self.boundary_mask(), order_by=order_by)

    def sea_points(self, order_by: str='lat'):
        return self.native_xy(mask=self.sea_mask(), order_by=order_by)

    def sea_points_xy(self, order_by: str='y'):
        return self.xy(mask=self.sea_mask(), order_by=order_by)

    def sea_points_lonlat(self, order_by: str='lat'):
        return self.lonlat(mask=self.sea_mask(), order_by=order_by)

    def land_points(self, order_by: str='lat'):
        return self.native_xy(mask=np.logical_not(self.sea_mask()), order_by=order_by)

    def land_points_xy(self, order_by: str='y'):
        return self.xy(mask=np.logical_not(self.sea_mask()), order_by=order_by)

    def land_points_lonlat(self, order_by: str='lat'):
        return self.lonlat(mask=np.logical_not(self.sea_mask()), order_by=order_by)
