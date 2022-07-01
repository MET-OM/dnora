import numpy as np
import xarray as xr
from copy import copy



class PowerSpectrum:
    def _reset_vars(self):
        self._set(np.zeros(self.size()),'spec')
        self._update_mask()

    def _update_mask(self) -> None:
        self._set(np.ones(self.grid_size()).astype(bool),'mask', only_grid_coords=True)

    def mask(self, logical=True) -> np.ndarray:
        """Returns bool array of the mask.
        Set logical=False to get 0 for bad spectra and 1 for good spectra. """

        mask = self._get('mask', default_data=np.ones(self.grid_size()).astype(int))
        if logical:
            mask = mask.astype(bool)
        else:
            mask = mask.astype(int)

        return mask

    def spec(self) -> np.ndarray:
        return self._get('spec')
