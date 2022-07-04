from point_skeleton import PointSkeleton
from gridded_skeleton import GriddedSkeleton
import numpy as np
from power_spectrum import PowerSpectrum
from coordinate_factory import add_time, add_frequency, add_direction
from datetime import datetime, timedelta
from mask_factory import add_mask
from datavar_factory import add_datavar
from enum import Enum, auto

class BoundaryConvention(Enum):
    """Conventions of BoundarySpectra (2D)"""

    """
    Oceanic convention
    Directional vector monotonically increasing.
    Direction to. North = 0, East = 90."""
    OCEAN = auto()

    """
    Meteorological convention
    Directional vector monotonically increasing.
    Direction from. North = 0, East = 90."""
    MET = auto()


    """
    Mathematical convention
    Directional vector of type: [90 80 ... 10 0 350 ... 100]
    Direction to. North = 90, East = 0."""
    MATH = auto()

    """
    Mathematical convention in vector
    Directional vector of type: [90 80 ... 10 0 350 ... 100]
    Direction to. North = 90, East = 0."""
    MATHVEC = auto()

    """
    WAVEWATCH III output convention
    Directional vector of type: [90 80 ... 10 0 350 ... 100]
    Direction to. North = 0, East = 90."""
    WW3 = auto()

class SpectralConvention(Enum):
    """Conventions of Spectra (1D)"""

    """
    Oceanic convention
    Direction to. North = 0, East = 90."""
    OCEAN = auto()

    """
    Meteorological convention
    Direction from. North = 0, East = 90."""
    MET = auto()

    """
    Mathematical convention
    Direction to. North = 90, East = 0."""
    MATH = auto()


@add_mask(name='bad', coords='all', default_value=0)
@add_datavar(name='spec', default_value=0.)
@add_frequency(grid_coord=False)
@add_time(grid_coord=True)
class Spectra(PointSkeleton):
    def __init__(self, grid=None, x=None, y=None, lon=None, lat=None, name: str="AnonymousSpectra"):
        self._name = name
        if grid is not None:
            x, y = grid.xy(strict=True)
            lon, lat = grid.lonlat(strict=True)
        t = np.arange(datetime(2020,1,1), datetime(2020,1,2), timedelta(hours=1)).astype(datetime)
        self._init_structure(x, y, lon, lat, time=t, freq=np.array([0.1,0.2,0.3]))


@add_mask(name='bad', coords='all', default_value=0)
@add_datavar(name='spec', default_value=0.)
@add_direction(grid_coord=False)
@add_frequency(grid_coord=False)
@add_time(grid_coord=False)
class Boundary(PointSkeleton):
    def __init__(self, grid=None, x=None, y=None, lon=None, lat=None, name: str="AnonymousSpectra"):
        self._name = name
        if grid is not None:
            x, y = grid.xy(strict=True)
            lon, lat = grid.lonlat(strict=True)
        t = np.arange(datetime(2020,1,1), datetime(2021,2,2), timedelta(hours=1)).astype(datetime)
        self._init_structure(x, y, lon, lat, time=t, freq=np.array([0.1,0.2,0.3]), dirs=np.array([0, 45, 90, 135, 180, 225, 270, 315]))
