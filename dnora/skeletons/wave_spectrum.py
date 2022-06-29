import numpy as np
import xarray as xr
from copy import copy



class WaveSpectrum:
    def _set_spec(self, spec: np.ndarray) -> None:
        self.merge_in_ds(self.compile_to_ds(spec,'spec'))

    def _reset_vars(self):
        ds_spec = self.compile_to_ds(np.zeros(self.size()),'spec')
        self.merge_in_ds(ds_spec)
        self._update_sea_mask()

    def _update_sea_mask(self) -> None:
        mask = np.ones(self.size()).astype(bool)
        self.merge_in_ds(self.compile_to_ds(mask,'mask'))

    def mask(self, logical=True) -> np.ndarray:
        """Returns bool array of the mask.
        Set logical=False to get 0 for bad spectra and 1 for good spectra. """

        if hasattr(self.data, 'mask'):
            mask = self.data.mask.values
        else:
            mask = np.full(self.size(), 1.)

        if logical:
            mask = mask.astype(bool)
        return mask

    def spec(self) -> np.ndarray:
        if hasattr(self.data, 'spec'):
            return self.data.spec.values
        return None

    def freq(self, angular=False):
        constant = 1.
        if angular:
            constant = 2*np.pi

        if hasattr(self, 'data') and hasattr(self.data, 'freq'):
            return self.data.freq.values*constant
        return None

class WaveSpectrum2D(WaveSpectrum):
    def dirs(self, radians=False):
        constant = 1.
        if radians:
            constant = np.pi/180

        if hasattr(self, 'data') and hasattr(self.data, 'dirs'):
            return self.data.dirs.values*constant
        return None
