import xarray as xr
from typing import List
import numpy as np

class WaveSeries:
    def from_spectra(self, spec, parameters):
        dsets = []
        for wp in parameters:
            dsets.append(wp(spec))
        self.data = xr.merge(dsets)
