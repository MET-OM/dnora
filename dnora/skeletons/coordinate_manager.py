from typing import Callable

class CoordinateManager:
    """Keeps track of coordinates and data variables that are added to classes
    by the decorators.

    Also contains list of strings of fixed coordinates used by Skeleton."""

    cartesian_coords = ['x', 'y']
    spherical_coords = ['lon', 'lat']
    cartesian_strings = ['x','y','xy']
    spherical_strings = ['lon','lat','lonlat']
    spatial_coords = ['y','x','lat','lon','inds']

    def __init__(self):
        self.coords = {}
        self.coords['grid'] = []
        self.coords['gridpoint'] = []
        self.coords['initial'] = []

        self.vars = {}
        self.vars['added'] = {}
        self.vars['initial'] = {}

        self.masks = {}
        self.masks['added'] = {}


    def add_var(self, name: str, coords: str, get_empty: Callable) -> None:
        self.vars['added'][name] = (coords, get_empty)

    def add_mask(self, name: str, coords: str, get_empty: Callable) -> None:
        self.masks['added'][name] = (coords, get_empty)

    def add_coord(self, name: str, grid_coord: bool) -> None:
        if grid_coord:
            self.coords['grid'].append(name)
        else:
            self.coords['gridpoint'].append(name)

    def add_initial_vars(self, initial_vars) -> None:
        if not isinstance(initial_vars, dict):
            raise ValueError('initial_vars needs to be a dict of tuples!')
        self.vars['initial'] = initial_vars

    def add_initial_coords(self, initial_coords) -> None:
        if not isinstance(initial_coords, list):
            raise ValueError('initial_coords needs to be a list of strings!')
        self.coords['initial'] = initial_coords

    def initial_vars(self):
        return self.vars['initial']

    def initial_coords(self):
        return self.coords['initial']

    def added_vars(self):
        return self.vars['added']

    def added_masks(self):
        return self.masks['added']

    def added_coords(self, type: str='all') -> list[str]:
        """Returns list of coordinates that have been added to the fixed
        Skeleton coords.

        'all': All added coordinates
        'grid': coordinates for the grid (e.g. z, time)
        'gridpoint': coordinates for a grid point (e.g. frequency, direcion or time)
        """
        if type not in ['all', 'grid', 'gridpoint']:
            print("Variable type needs to be 'all', 'grid' or 'gridpoint'.")
            return None

        if type == 'all':
            return self.added_coords('grid') + self.added_coords('gridpoint')
        return self.coords[type]
