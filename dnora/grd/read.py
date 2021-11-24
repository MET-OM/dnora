import xarray as xr
from copy import copy
from abc import ABC, abstractmethod
from ..aux import expand_area
import numpy as np
from typing import Tuple

class TopoReader(ABC):
    """Abstract class for reading the bathymetry.

    The dimensions and orientation of the bathymetry array that is returned to
    the object should be:

    rows = latitude and colums = longitude (i.e.) shape = (nr_lat, nr_lon).

    North = [-1,:]
    South = [0,:]
    East = [:,-1]
    West = [:,0]
    """

    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float):
        """Reads the bathymetrical information from a source and returns the data.

        This method is called from within the Grid-object
        """
        pass

    @abstractmethod
    def __str__(self):
        """Describes what topography is read and from where.

        This is called by the Grid-objeect to provide output to the user.
        """
        pass


class EmptyTopo(TopoReader):
    """Creates an empty topography. Called when setting initial spacing."""
    def __init__(self, grid):
        self.grid = grid
        pass

    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float):
        """Creates a trivial topography with all water points."""
        topo = np.ones((self.grid.data.ny,self.grid.data.nx))*9999
        topo_lon = copy(self.grid.lon())
        topo_lat = copy(self.grid.lat())
        return topo, topo_lon, topo_lat

    def __str__(self):
        return("Creating an empty topography with depth values 9999.")

class EMODNET2018(TopoReader):
    """Reads data from EMODNET bathymetry.

    The data needs to be downloaded from https://portal.emodnet-bathymetry.eu/,
    since no API to the database exists.
    """

    def __init__(self, expansion_factor: float=1.2, tile: str='C5', folder: str='/lustre/storeB/project/fou/om/WW3/bathy/emodnet_115m_x_115m') -> Tuple:
        self.source=f'{folder}/{tile}_2018.dtm'
        self.expansion_factor = expansion_factor
        return

    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float):
        # Area is expanded a bit to not get in trouble in the meshing stage
        # when we interpoolate or filter
        lon0, lon1, lat0, lat1 = expand_area(lon_min, lon_max, lat_min, lat_max, self.expansion_factor)

        ds = xr.open_dataset(self.source).sel(COLUMNS=slice(lon0, lon1), LINES=slice(lat0, lat1))

        topo = ds.DEPTH.values

        # Set depth to positive values and land to -999
        land_mask = topo > 0
        topo = -1*topo
        topo[land_mask] = -999

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

    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float) -> Tuple:
        # Just use the values it was forcefed on initialization
        topo = copy(self.topo)
        topo_lon = copy(self.topo_lon)
        topo_lat = copy(self.topo_lat)

        return topo, topo_lon, topo_lat

    def __str__(self):
        return("Passing on the topography I was initialized with.")
