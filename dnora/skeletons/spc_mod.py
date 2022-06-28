from point_skeleton import PointSkeleton
import numpy as np

class Spectra(PointSkeleton):
    def __init__(self, grid, name: str="AnonymousSpectra"):
        self.grid = copy(grid)
        self._name = copy(name)
        self._convention = None
        self._history = []

    def freq(self, angular=False):
        constant = 1.
        if angular:
            constant = 2*np.pi

        if hasattr(self, 'data') and hasattr(self.data, 'freq'):
            return self.data.freq.values*constant
        else:
            return None
