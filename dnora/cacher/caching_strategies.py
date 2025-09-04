from enum import Enum, auto
from .patching_functions import single_patch, patch_in_time


class CachingStrategy(Enum):
    """Strategies for the cacher to read in data that was not found in cahce

    SinglePatch:    Find the minimum extent in time and space that covers the missing data
                    Great for readers that have a large overhead when calles (e.g. ERA5).

    PatchInTime:    Find patches in time (possibly several) that cover the times when we have ANY missing data in space.
                    Works optimally if cached data is extended in time only, but results in re-reading all data
                    if even one spatial tile is missing.

    PatchInSpace:   Finds patches in time (possibly several) that covers the spatial tiles where we have ANY missing data in time.
                    Works optimally if cached data is extended in space only, but results in re-reading all data
                    if even one hour is missing for all tiles.

    TileByTile:     Reads in every space-time tile separately. Guaranteed to not re-read any existing data, but
                    can result in a lot of individual calls to the reader.

    DontCacheMe:    Disables all caching functionality for particular reader. Useful for readers that create data (e.g. ConstantGriddedData)
                    or readers that can have variable sources (e.g. Netcdf or SpectraToWaveSeries)
    """

    SinglePatch = auto()
    PatchInTime = auto()
    PatchInSpace = auto()
    TileByTile = auto()
    DontCacheMe = auto()


caching_strategies = {
    CachingStrategy.SinglePatch: single_patch,
    CachingStrategy.PatchInTime: patch_in_time,
}
