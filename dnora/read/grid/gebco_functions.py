
from dnora.cacher.tiles import create_tiles


def get_covering_tiles(area, expansion_factor):
    """Finds the tiles needed to cover the given area"""
    lons, lats, __ = create_tiles(area, '2020-01-01', '2020-01-01', expansion_factor)
    return lons, lats
    
def get_tile_names(lons, lats):
    """Create tile names from longitude and latitude tuples"""
    tiles = []
    for lo, la  in zip(lons, lats):
        tiles.append(f"{lo[0]:03.0f}E_{la[0]:03.0f}N")

    return tiles


