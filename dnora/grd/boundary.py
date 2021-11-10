import numpy as np
from typing import List
from .. import msg
from .grd_mod import BoundarySetter # Abstract class

# BoundarySetter used by the grid object
from .grd_mod import ClearBoundary

class EdgesAsBoundary(BoundarySetter):
    def __init__(self, edges: List[str]=['N', 'S', 'E', 'W'], step: int=1):
        self.edges = edges
        if step < 1:
            raise ValueError('step cannot be smaller than 1')
        else:
            self.step = int(step)
        return

    def __call__(self, mask_size: tuple):
        boundary_mask = np.full(mask_size, False)
        # --------- North boundary ----------
        if 'N' in self.edges:
            boundary_mask[-1,::self.step] = True
        ## --------- South boundary ----------
        if 'S' in self.edges:
            boundary_mask[0,::self.step] = True
        ## --------- East boundary ----------
        if 'E' in self.edges:
            boundary_mask[::self.step,-1] = True
        ## --------- West boundary ----------
        if 'W' in self.edges:
            boundary_mask[::self.step,0] = True

        return boundary_mask

    def __str__(self):
        return(f"Setting all edges {self.edges} to boundary points using step {self.step}.")


class SetMatrix(BoundarySetter):
    def __init__(self, matrix):
        self.matrix = matrix
        return

    def __call__(self, mask_size: tuple):
        if self.matrix.shape == mask_size:
            return self.matrix
        else:
            raise Exception(f'Given mask for boundary points does not match the dimensions of the grid ({self.matrix.shape[0]}x{self.matrix.shape[1]} vs {mask_size[0]}x{mask_size[1]})')

    def __str__(self):
        return(f"Setting boundary points using the boolean matrix I was initialized with.")
