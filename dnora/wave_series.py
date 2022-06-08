import xarray as xr
from typing import List
import numpy as np
import matplotlib.pyplot as plt
from typing import List
from copy import copy
import pandas as pd
from dnora.wave_parameters import WaveParameter

class WaveSeries:
    def __call__(self, parameter: str) -> np.ndarray:
        return self.data[parameter].values

    def from_spectra(self, spec: xr.Dataset, parameters: List[WaveParameter]):
        dsets = []
        for wp in parameters:
            dsets.append(wp(spec))
        self.data = xr.merge(dsets)

    def plot(self, parameter: str, x: List[int]=None) -> np.ndarray:
        self.slice_data(x=x)[parameter].plot()
        plt.show()

    def parameters(self) -> List[str]:
        return list(self.data.data_vars)

    def slice_data(self, start_time: str='', end_time: str='', x: List[int]=None):
        """Slice data in space (x) and time. Returns an xarray dataset."""

        if x is None:
            x=self.x()

        if not start_time:
            # This is not a string, but slicing works also with this input
            start_time = self.time()[0]

        if not end_time:
            # This is not a string, but slicing works also with this input
            end_time = self.time()[-1]

        sliced_data = self.data.sel(time=slice(start_time, end_time), x = x)

        return sliced_data

    def x(self):
        if hasattr(self, 'data'):
            return copy(self.data.x.values)
        else:
            return np.array([])

    def time(self):
        return pd.to_datetime(self.data.time.values)
