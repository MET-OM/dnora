import numpy as np
from coordinate_manager import CoordinateManager
from functools import partial
def add_mask(name, coords, default_value):
    def mask_decorator(c):
        def get_mask(self, empty: bool=False) -> np.ndarray:
            """Returns bool array of the mask.

            Set boolean=False to get 0 for land and 1 for sea.
            Set empty=True to get an empty mask (even if it doesn't exist)"""

            if empty:
                return np.full(self.size(coords), default_value).astype(bool)

            mask = self.ds_manager.get(f'{name}_mask')

            if mask is None:
                return None

            return mask.astype(bool)

        def get_masked_points(self, type: str='native', order_by: str='lat', strict=False):
            mask = self.get(f'{name}_mask')

            if type == 'native':
                return self.native_xy(mask=mask, order_by=order_by)
            elif type in skeleton_strings['cartesian_strings']:
                return self.xy(mask=mask, order_by=order_by, strict=strict)
            elif type in skeleton_strings['spherical_strings']:
                return self.lonlat(mask=mask, order_by=order_by, strict=strict)

        if not hasattr(c, '_coord_manager'):
            c._coord_manager =  CoordinateManager()
        c._coord_manager.add_mask(name, coords)
        exec(f'c.{name}_mask = get_mask')
        exec(f'c.{name}_points = get_masked_points')

        return c

    return mask_decorator
