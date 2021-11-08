import numpy as np
from .. import msg

from .grd_mod import GridProcessor #Abstract class

# GridModifiers used as defaults in the Grid object methods
from .grd_mod import TrivialFilter


class SetMinDepth(GridProcessor):
    def __init__(self, min_depth, to_land = -1):
        self.to_land = to_land
        self.min_depth = min_depth
        return
    def __call__(self, data, lon, lat, land_sea_mask, boundary_mask):
        shallow_points = data > self.min_depth

        if self.to_land >= 0:
            msg.info(f"Setting points shallower than {self.min_depth} to land ({self.to_land})")
            new_data = copy(data)
            new_data[np.logical_and(shallow_points, land_sea_mask)] = self.to_land # Don't touch land points by usign self.mask
            msg.info(f"Affected {np.count_nonzero(np.logical_and(shallow_points, land_sea_mask))} points")
        else:
            msg.info(f"Setting points shallower than {self.min_depth} to {self.min_depth}")
            # Set points to the limiter
            new_data = copy(data)
            new_data[np.logical_and(shallow_points, land_sea_mask)] = self.min_depth # Don't touch land points by usign self.mask
            msg.info(f"Affected {np.count_nonzero(np.logical_and(shallow_points, land_sea_mask))} points")
        return new_data

