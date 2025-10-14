from enum import Enum
from typing import Union


class SpectralConvention(Enum):
    """Conventions of BoundarySpectra (2D)"""

    """
    Oceanic convention
    Directional vector monotonically increasing.
    Direction to. North = 0, East = 90."""
    OCEAN = "ocean"

    """
    Meteorological convention
    Directional vector monotonically increasing.
    Direction from. North = 0, East = 90."""
    MET = "met"

    """
    Mathematical convention
    Directional vector monotonically increasing.
    Direction to. North = 90, East = 0."""
    MATH = "math"

    """
    Mathematical convention in vector
    Directional vector of type: [90 80 ... 10 0 350 ... 100]
    Direction to. North = 90, East = 0."""
    MATHVEC = "mathvec"

    """
    WAVEWATCH III output convention
    Directional vector of type: [90 80 ... 10 0 350 ... 100]
    Direction to. North = 0, East = 90."""
    WW3 = "ww3"

    UNDEFINED = "undefined"


def spectral_convention_from_string(
    convention_str: Union[str, SpectralConvention],
) -> SpectralConvention:
    if convention_str is None:
        return None
    if isinstance(convention_str, SpectralConvention):
        return convention_str
    return SpectralConvention[convention_str.upper()]
