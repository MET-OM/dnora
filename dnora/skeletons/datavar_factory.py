import numpy as np
from coordinate_manager import CoordinateManager
from functools import partial
def add_datavar(name, coords='grid', default_value=0., stash_get=False):
    def datavar_decorator(c):
        def get_var(self, empty: bool=False) -> np.ndarray:
            """Returns the data variable.

            Set empty=True to get an empty data variable (even if it doesn't exist)"""

            if empty:
                return np.full(self.size(coords), default_value)

            return self.ds_manager.get(name)

        if not hasattr(c, '_coord_manager'):
            c._coord_manager =  CoordinateManager()

        c._coord_manager.add_var(name, coords)

        if stash_get:
            exec(f'c._{name} = get_var')
        else:
            exec(f'c.{name} = get_var')

        return c

    return datavar_decorator
