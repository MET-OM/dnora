import numpy as np
from .coordinate_manager import CoordinateManager
from functools import partial
def add_datavar(name, coords='grid', default_value=0., stash_get=False, aftermath=False):
    """stash_get = True means that the coordinate data can be accessed
    by method ._{name}() instead of .{name}()

    This allows for alternative definitions of the get-method elsewere."""

    def datavar_decorator(c):
        def get_var(self, empty: bool=False, data_array: bool=False, **kwargs) -> np.ndarray:
            """Returns the data variable.

            Set empty=True to get an empty data variable (even if it doesn't exist).

            **kwargs can be used for slicing data.
            """
            if not self._structure_initialized():
                return None
            if empty:
                return np.full(self.size(coords, **kwargs), default_value)

            data = self.ds_manager.get(name, **kwargs)
            if data_array:
                return data.copy()
            return data.values.copy()
    
        if not hasattr(c, '_coord_manager'):
            c._coord_manager =  CoordinateManager()

        c._coord_manager.add_var(name, coords)

        if stash_get:
            exec(f'c._{name} = get_var')
        else:
            if aftermath:
                exec(f'c.{name} = partial(get_var, c)')
            else:
                exec(f'c.{name} = get_var')

        return c

    return datavar_decorator
