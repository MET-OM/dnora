import numpy as np

def add_datavar(name, coords, default_value):
    def datavar_decorator(c):
        def get_var(self, empty: bool=False) -> np.ndarray:
            """Returns the data variable.

            Set empty=True to get an empty data variable (even if it doesn't exist)"""

            if empty:
                return np.full(self.size(coords), default_value)

            return self._get(name)

        if not hasattr(c, '_datavar_dict'):
            c._datavar_dict = {}

        c._datavar_dict[name] = (coords, default_value)
        exec(f'c.{name} = get_var')
        #exec(f'c.{name}_points = get_masked_points')
        #exec(f'c._update_{name} = update_var')

        return c

    return datavar_decorator
