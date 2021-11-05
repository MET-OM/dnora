from abc import ABC, abstractmethod
import xarray as xr
from copy import copy
import numpy as np

class TopoFetcher(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, lon_min, lon_max, lat_min, lat_max):
        pass


class EMODNET2018(TopoFetcher):
    """Reads data from EMODNET"""
    def __init__(self, expansion_factor = 1.2, tile = 'C5', folder = '/lustre/storeB/project/fou/om/WW3/bathy/emodnet_115m_x_115m'):
        self.source=f'{folder}/{tile}_2018.dtm'
        self.expansion_factor = expansion_factor
        return

    def __call__(self, lon_min, lon_max, lat_min, lat_max):
        # If we limit ourselves to exactly the grid, we will get nans at the edges in the interpolation. Add 10% tolerance around all edges.
        tolerance_lon = (lon_max-lon_min)*(self.expansion_factor - 1)*0.5
        tolerance_lat = (lat_max-lat_min)*(self.expansion_factor - 1)*0.5

        ds = xr.open_dataset(self.source).sel(COLUMNS=slice(lon_min-tolerance_lon, lon_max+tolerance_lon), LINES=slice(lat_min-tolerance_lat, lat_max+tolerance_lat))
        topo = ds.DEPTH.values
        topo_lon = ds.COLUMNS.values
        topo_lat = ds.LINES.values
        return topo, topo_lon, topo_lat


class EmptyTopo(TopoFetcher):
    """Creates an empty topography. Called when setting initial spacing."""
    def __init__(self, grid):
        self.grid = copy(grid)
        pass

    def __call__(self, lon_min, lon_max, lat_min, lat_max):
        # Creates a trivial topography with all water points
        topo = np.ones((self.grid.data.ny,self.grid.data.nx))*-9999
        topo_lon = copy(self.grid.lon())
        topo_lat = copy(self.grid.lat())
        return topo, topo_lon, topo_lat


class ForceFeed(TopoFetcher):
    """Simply passes on the data it was fed upon initialization"""
    def __init__(self, topo, topo_lon, topo_lat):
        self.topo = copy(topo)
        self.topo_lon = copy(topo_lon)
        self.topo_lat = copy(topo_lat)
        return

    def __call__(self, lon_min, lon_max, lat_min, lat_max):
        # Just use the values it was forcefed on initialization
        topo = copy(self.topo)
        topo_lon = copy(self.topo_lon)
        topo_lat = copy(self.topo_lat)
        return topo, topo_lon, topo_lat
