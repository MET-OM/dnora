import xarray as xr
from copy import copy

from .grd_mod import TopoReader # Abstract class

# Readers used as defaults in the Grid object methods
from .grd_mod import EmptyTopo

class EMODNET2018(TopoReader):
    """Reads data from EMODNET"""
    def __init__(self, expansion_factor = 1.2, tile = 'C5', folder = '/lustre/storeB/project/fou/om/WW3/bathy/emodnet_115m_x_115m'):
        self.source=f'{folder}/{tile}_2018.dtm'
        self.expansion_factor = expansion_factor
        return

    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float):
        # If we limit ourselves to exactly the grid, we will get nans at the edges in the interpolation. Add 10% tolerance around all edges.
        tolerance_lon = (lon_max-lon_min)*(self.expansion_factor - 1)*0.5
        tolerance_lat = (lat_max-lat_min)*(self.expansion_factor - 1)*0.5

        ds = xr.open_dataset(self.source).sel(COLUMNS=slice(lon_min-tolerance_lon, lon_max+tolerance_lon), LINES=slice(lat_min-tolerance_lat, lat_max+tolerance_lat))
        topo = ds.DEPTH.values
        topo_lon = ds.COLUMNS.values
        topo_lat = ds.LINES.values
        return topo, topo_lon, topo_lat

    def __str__(self):
        return(f"Reading EMODNET topography from {self.source}.")


class ForceFeed(TopoReader):
    """Simply passes on the data it was fed upon initialization"""
    def __init__(self, topo, topo_lon, topo_lat):
        self.topo = copy(topo)
        self.topo_lon = copy(topo_lon)
        self.topo_lat = copy(topo_lat)
        return

    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float):
        # Just use the values it was forcefed on initialization
        topo = copy(self.topo)
        topo_lon = copy(self.topo_lon)
        topo_lat = copy(self.topo_lat)
        return topo, topo_lon, topo_lat

    def __str__(self):
        return("Passing on the topography I was initialized with.")
