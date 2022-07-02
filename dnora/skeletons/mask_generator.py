import numpy as np

def add_mask(name, coords, default_value):
    def mask_decorator(c):
        def get_mask(self, boolean: bool=None, empty: bool=False) -> np.ndarray:
            """Returns bool array of the mask.

            Set boolean=False to get 0 for land and 1 for sea.
            Set empty=True to get an empty mask (even if it doesn't exist)"""

            if boolean is None:
                # If not specified take it from the default value
                boolean = isinstance(default_value, bool)

            data_type = 'bool'*boolean or 'float'

            if empty:
                return np.full(self.size(coords), default_value).astype(data_type)

            mask = self._get(f'{name}_mask')

            if mask is None:
                return None

            return mask.astype(data_type)

        def get_masked_points(self, type: str='native', order_by: str='lat', strict=False):
            mask = getattr(self, f'{name}_mask')(boolean=True)

            if type == 'native':
                return self.native_xy(mask=mask, order_by=order_by)
            elif type in self._cartesian_strings:
                return self.xy(mask=mask, order_by=order_by, strict=strict)
            elif type in self._spherical_strings:
                return self.lonlat(mask=mask, order_by=order_by, strict=strict)

        if not hasattr(c, '_mask_dict'):
            c._mask_dict = {}

        c._mask_dict[name] = (coords, default_value)
        exec(f'c.{name}_mask = get_mask')
        exec(f'c.{name}_points = get_masked_points')
        #exec(f'c._update_{name}_mask = update_mask')

        return c

    return mask_decorator
