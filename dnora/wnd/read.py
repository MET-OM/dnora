from abc import ABC,  abstractmethod
from copy import copy
import numpy as np
import xarray as xr
from subprocess import call
from .. import msg
from ..aux import create_time_stamps, u_v_from_dir, expand_area, lon_in_km

import os, glob

#from .wnd_mod import Forcing # Forcing object

from ..grd.grd_mod import Grid # Grid object

class ForcingReader(ABC):
    """Reads forcing data from some source and provide it to the object.

    The area is defined from the Grid object that is passed.
    """

    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, grid: Grid, start_time: str, end_time: str, expansion_factor: float):
        pass


