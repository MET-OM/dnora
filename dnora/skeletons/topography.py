import numpy as np
import xarray as xr
from copy import copy

class Topography:
    def _set_topo(self, topo: np.ndarray) -> None:
        self.merge_in_ds(self.compile_to_ds(topo,'topo'))

    def _reset_vars(self):
        ds_topo = self.compile_to_ds(999*np.ones(self.size()),'topo')
        ds_bnd = self.compile_to_ds(np.zeros(self.size()),'boundary_mask')
        self.merge_in_ds([ds_topo, ds_bnd])
        self._update_sea_mask()

    def _update_sea_mask(self) -> None:
        sea_mask = np.logical_not(np.isnan(self.topo())).astype(int)
        self.merge_in_ds(self.compile_to_ds(sea_mask,'sea_mask'))

    def sea_mask(self, logical=True) -> np.ndarray:
        """Returns bool array of the sea mask.
        Set logical=False to get 0 for land and 1 for sea. """

        if hasattr(self.data, 'sea_mask'):
            mask = self.data.sea_mask.values
        else:
            mask = np.full(self.size(), 1.)

        if logical:
            mask = mask.astype(bool)
        return mask

    def boundary_mask(self, logical: bool=True) -> np.ndarray:
        """Returns bool array of the boundary mask.
        Set logical=False to get 1 for boundary points """

        if hasattr(self.data, 'boundary_mask'):
            mask = self.data.boundary_mask.values
        else:
            mask = np.full(self.size(), 1.)

        if logical:
            mask = mask.astype(bool)
        return mask


    def topo(self, land: float=-999) -> np.ndarray:
        if hasattr(self.data, 'topo'):
            topo = copy(self.data.topo.values)
            topo[np.logical_not(self.sea_mask())] = land
            return topo
        return None

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
