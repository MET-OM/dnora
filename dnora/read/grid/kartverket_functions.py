from .tiling_functions import read_limits

def find_tiles(lon, lat):
    tiles, lons, lats = read_limits("kartverket50m_tile_coords.txt")
    needed_tiles = []

    for t, lo, la in zip(tiles, lons, lats):
        lon_match = lo[0]<lon[1] and lo[1]>lon[0]
        lat_match = la[0]<lat[1] and la[1]>lat[0]
        if lon_match and lat_match:
            needed_tiles.append(t)

    return needed_tiles