
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


def get_url(year: int):
    """Returns correct opendap url for given year"""
    if year < 2019:
        raise ValueError(f"GEBCO datasets only available for years 2019 onwards!")
    if year in [2021,2022,2023,2024]:
        return f"https://dap.ceda.ac.uk/thredds/dodsC/bodc/gebco/global/gebco_{year}/ice_surface_elevation/netcdf/GEBCO_{year}_CF.nc"
    else:
        return f"https://dap.ceda.ac.uk/thredds/dodsC/bodc/gebco/global/gebco_{year}/ice_surface_elevation/netcdf/GEBCO_{year}.nc"
    
    