import numpy as np
from coordinate_manager import CoordinateManager
from functools import partial
def add_mask(name: str, coords: str, default_value: int, opposite_name: str=None):
    def mask_decorator(c):
        def get_mask(self, empty: bool=False, **kwargs) -> np.ndarray:
            """Returns bool array of the mask.

            Set boolean=False to get 0 for land and 1 for sea.
            Set empty=True to get an empty mask (even if it doesn't exist)

            **kwargs can be used for slicing data.
            """

            if empty:
                return np.full(self.size(coords, **kwargs), default_value).astype(bool)
            mask = self.ds_manager.get(f'{name}_mask', **kwargs).values.copy()

            if mask is None:
                return None

            return mask.astype(bool)

        def get_not_mask(self, **kwargs):
            mask = get_mask(self, **kwargs)
            if mask is None:
                return None
            return np.logical_not(mask)

        def get_masked_points(self, type: str='native', order_by: str='lat', strict=False, **kwargs):
            mask = self.get(f'{name}_mask', **kwargs).copy()

            if type == 'native':
                return self.native_xy(mask=mask, order_by=order_by, **kwargs)
            elif type in self._coord_manager.cartesian_strings:
                return self.xy(mask=mask, order_by=order_by, strict=strict, **kwargs)
            elif type in self._coord_manager.spherical_strings:
                return self.lonlat(mask=mask, order_by=order_by, strict=strict, **kwargs)

        def get_not_points(self, type: str='native', order_by: str='lat', strict=False, **kwargs):
            mask = np.logical_not(self.get(f'{name}_mask', **kwargs).copy())

            if type == 'native':
                return self.native_xy(mask=mask, order_by=order_by, **kwargs)
            elif type in self._coord_manager.cartesian_strings:
                return self.xy(mask=mask, order_by=order_by, strict=strict, **kwargs)
            elif type in self._coord_manager.spherical_strings:
                return self.lonlat(mask=mask, order_by=order_by, strict=strict, **kwargs)


        if not hasattr(c, '_coord_manager'):
            c._coord_manager =  CoordinateManager()
        c._coord_manager.add_mask(name, coords)
        exec(f'c.{name}_mask = get_mask')
        exec(f'c.{name}_points = get_masked_points')
        if opposite_name is not None:
            exec(f'c.{opposite_name}_mask = get_not_mask')
            exec(f'c.{opposite_name}_points = get_not_points')
        return c

    return mask_decorator
