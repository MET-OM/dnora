from enum import Enum, auto

class SpectralConvention(Enum):
    """Conventions of BoundarySpectra (2D)"""

    """
    Oceanic convention
    Directional vector monotonically increasing.
    Direction to. North = 0, East = 90."""
    OCEAN = 'ocean'

    """
    Meteorological convention
    Directional vector monotonically increasing.
    Direction from. North = 0, East = 90."""
    MET = 'met'

    """
    Mathematical convention
    Directional vector of type: [90 80 ... 10 0 350 ... 100]
    Direction to. North = 90, East = 0."""
    MATH = 'math'

    """
    Mathematical convention in vector
    Directional vector of type: [90 80 ... 10 0 350 ... 100]
    Direction to. North = 90, East = 0."""
    MATHVEC = 'mathvec'

    """
    WAVEWATCH III output convention
    Directional vector of type: [90 80 ... 10 0 350 ... 100]
    Direction to. North = 0, East = 90."""
    WW3 = 'ww3'


def convention_from_string(convention: str) -> SpectralConvention:
    """Returns the SpectralConvention corresponding to a string.

    Possible (case insensitive) strings are:
    ['ocean', 'met', 'math', 'mathvec', 'ww3']

    This is implemented so that the user can set convetnions without having
    to import the Enum.

    DNORA internally should use only Enums.
    """
    if convention is None:
        return None

    if isinstance(convention, SpectralConvention):
        return convention

    convention = convention.lower()
    if convention == SpectralConvention.OCEAN.value:
        return SpectralConvention.OCEAN
    elif convention == SpectralConvention.MET.value:
        return SpectralConvention.MET
    elif convention == SpectralConvention.MATH.value:
        return SpectralConvention.MATH
    elif convention == SpectralConvention.MATHVEC.value:
        return SpectralConvention.MATHVEC
    elif convention == SpectralConvention.WW3.value:
        return SpectralConvention.WW3
    else:
        return None


def convert_2d_to_1d(convention_2d: SpectralConvention):
    """Determine the correct 1d convention from a 2d convention.

    If a 1d convention is given, the identical convetnion is returned."""
    if convention_2d in [SpectralConvention.OCEAN, SpectralConvention.WW3]:
        convention_1d = SpectralConvention.OCEAN
    elif convention_2d in [SpectralConvention.MATH, SpectralConvention.MATHVEC]:
        convention_1d = SpectralConvention.MATH
    else:
        convention_1d = SpectralConvention.MET

    return convention_1d
