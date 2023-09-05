from __future__ import annotations

import xarray as xr
import pandas as pd
import utm
from copy import copy
from abc import ABC, abstractmethod
import numpy as np
from typing import Tuple, TYPE_CHECKING
if TYPE_CHECKING:
    from .grd_mod import Grid
# Import aux_funcsiliry functions
from ..aux_funcs import expand_area
import warnings
from .. import msg

class TopoReader(ABC):
    """Abstract class for reading the bathymetry. """

    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float):
        """Reads the bathymetrical information from a source and returns the data.

        !!!! DEPTH VALUES ARE POSITIVE AND EVERYTHING ELSE (including 0)
        IS INTERPRETED AS LAND. LON, LAT IN WGS84 !!!



        Two format options are available:

        1) Matrix, where topo_lon, topo_lat are vectors and topo is the
        corresponding matrix.

        The dimensions and orientation of the bathymetry array that is returned to
        the object should be:

        rows = latitude and colums = longitude (i.e.) shape = (nr_lat, nr_lon).

        North = [-1,:]
        South = [0,:]
        East = [:,-1]
        West = [:,0]

        2) Xyz, where topo, topo_lon, topo_lat are all vectors with the same length.

        ----

        This method is called from within the Grid-object
        """
        return topo, topo_lon, topo_lat

    @abstractmethod
    def __str__(self):
        """Describes what topography is read and from where.

        This is called by the Grid-object to provide output to the user.
        """
        pass


class EmptyTopo(TopoReader):
    """Creates an empty topography. Called when setting initial spacing."""
    def __init__(self, grid: Grid=None, nx=0, ny=0):
        if grid is not None:
            self.size = (grid.ny(),grid.nx())
        else:
            self.size = (ny, nx)
        pass

    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float):
        """Creates a trivial topography with all water points."""
        topo = np.ones(self.size)*9999
        topo_lon = np.linspace(lon_min, lon_max, self.size[1])
        topo_lat = np.linspace(lat_min, lat_max, self.size[0])
        return topo, topo_lon, topo_lat

    def __str__(self):
        return("Creating an empty topography with depth values 9999.")

class EMODNET2018(TopoReader):
    """Reads data from EMODNET bathymetry.

    The data needs to be downloaded from https://portal.emodnet-bathymetry.eu/,
    since no API to the database exists.
    """

    def __init__(self, expansion_factor: float=1.2, tile: str='C5', folder: str='/lustre/storeB/project/fou/om/WW3/bathy/emodnet_115m_x_115m'):
        self.source=f'{folder}/{tile}_2018.dtm'
        self.expansion_factor = expansion_factor
        return

    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float) -> Tuple:
        # Area is expanded a bit to not get in trouble in the meshing stage
        # when we interpoolate or filter
        lon0, lon1, lat0, lat1 = expand_area(lon_min, lon_max, lat_min, lat_max, self.expansion_factor)

        ds = xr.open_dataset(self.source).sel(COLUMNS=slice(lon0, lon1), LINES=slice(lat0, lat1))

        topo = ds.DEPTH.values

        # Negative valies and NaN's are land
        topo = -1*topo

        topo_lon = ds.COLUMNS.values
        topo_lat = ds.LINES.values

        return topo, topo_lon, topo_lat

    def __str__(self):
        return(f"Reading EMODNET topography from {self.source}.")

class KartverketNo50m(TopoReader):
    """Reads data from Kartverket bathymetry.

        High resolution bathymetry dataset for the whole Norwegian Coast.
        Can be found at:
        https://kartkatalog.geonorge.no/metadata/dybdedata-terrengmodeller-50-meters-grid-landsdekkende/bbd687d0-d34f-4d95-9e60-27e330e0f76e

        For reading several files at once, supply the 'tile' argument with a glob pattern, e.g. 'B*'.

        Contributed by: https://github.com/emiliebyer
        """

    def __init__(self, expansion_factor: float=1.2, utmzone: int=33,
                 tile: str='B1008', folder: str='/lustre/storeB/project/fou/om/WW3/bathy/kartverket_50m_x_50m') -> Tuple:
        self.source=f'{folder}/{tile}_grid50_utm33.xyz'
        self.expansion_factor = expansion_factor
        self.utmzone = utmzone

        return

    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float):
        # Area is expanded a bit to not get in trouble in the meshing stage
        # when we interpoolate or filter
        lon0, lon1, lat0, lat1 = expand_area(lon_min, lon_max, lat_min, lat_max, self.expansion_factor)

        print(f'Expansion factor: {self.expansion_factor}')

        import dask.dataframe as dd
        df = dd.read_csv(self.source, sep= ' ', header=None)
        df.columns = ['x','y','z']
        x = np.array(df['x'].astype(float))
        y = np.array(df['y'].astype(float))
        z = np.array(df['z'].astype(float))

        # Converting from utm to latitude and longitude
        lat, lon = utm.to_latlon(x, y, self.utmzone, northern=True, strict=False)

        # Applying given max and min values for lat and lon
        mask_lat = np.logical_and(lat0 <= lat, lat <= lat1)
        mask_lon = np.logical_and(lon0 <= lon, lon <= lon1)
        mask = np.logical_and(mask_lat, mask_lon)

        #topo_lon = lon[mask]
        x = x[mask]
        y = y[mask]
        topo = z[mask]

        # Adding NaN values for land area
        x_grid = np.unique(x) #np.arange(min(x), max(x), 50)
        y_grid = np.unique(y) #np.arange(min(y), max(y), 50)
        grid = np.empty((len(y_grid),len(x_grid)))
        grid[:] = np.nan

        # Adding depth values for all defined points, leaving the rest as NaN
        for i in range(len(topo)):
            xpos = np.where(x_grid==x[i])[0][0]
            ypos = np.where(y_grid==y[i])[0][0]
            grid[ypos][xpos] = topo[i]

        # Converting grid back to list of points usind Pandas
        df = pd.DataFrame(grid)
        df = df.unstack().reset_index()
        df['x'] = [x_grid[i] for i in df['level_0']]
        df['y'] = [y_grid[i] for i in df['level_1']]

        topo_lat, topo_lon = utm.to_latlon(df['x'], df['y'], 33, northern=True, strict = False)
        topo = df[0]

        # Masking again to get rid of additional values
        mask_lat = np.logical_and(lat0 < topo_lat, topo_lat < lat1)
        mask_lon = np.logical_and(lon0 < topo_lon, topo_lon < lon1)
        mask = np.logical_and(mask_lat, mask_lon)

        topo_lat = topo_lat[mask]
        topo_lon = topo_lon[mask]
        topo = topo[mask]


        return topo, topo_lon, topo_lat

    def __str__(self):
        return(f"Reading Kartverket topography from {self.source}.")

class GEBCO2021(TopoReader):
    """ Reads the GEBCO_2021 gridded bathymetric data set.
        A global terrain model for ocean and land,
        providing elevation data, in meters, on a 15 arc-second interval grid.
        Reference: GEBCO Compilation Group (2021) GEBCO 2021 Grid (doi:10.5285/c6612cbe-50b3-0cff- e053-6c86abc09f8f)
        Data (in netCDF format) can be downloaded here: https://www.gebco.net/data_and_products/gridded_bathymetry_data/

        Contributed by: https://github.com/emiliebyer
    """
    def __init__(self, expansion_factor: float=1.2, tile: str='n61.0_s59.0_w4.0_e6.0', folder: str='/lustre/storeB/project/fou/om/WW3/bathy/gebco2021'):
        self.source=f'{folder}/gebco_2021_{tile}.nc'
        self.expansion_factor = expansion_factor
        return

    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float):
        # Area is expanded a bit to not get in trouble in the meshing stage
        # when we interpoolate or filter
        lon0, lon1, lat0, lat1 = expand_area(lon_min, lon_max, lat_min, lat_max, self.expansion_factor)

        ds = xr.open_dataset(self.source).sel(lon=slice(lon0, lon1), lat=slice(lat0, lat1))

        elevation = ds.elevation.values.astype(float)

        # Negative valies and NaN's are land
        topo = -1*elevation

        topo_lon = ds.lon.values.astype(float)
        topo_lat = ds.lat.values.astype(float)

        return topo, topo_lon, topo_lat

    def __str__(self):
        return(f"Reading GEBCO2021 topography from {self.source}.")

class GEBCO2022(TopoReader):
    """ Reads the GEBCO_2021 gridded bathymetric data set.
        A global terrain model for ocean and land,
        providing elevation data, in meters, on a 15 arc-second interval grid.
        Reference: GEBCO Compilation Group (2021) GEBCO 2021 Grid (doi:10.5285/c6612cbe-50b3-0cff- e053-6c86abc09f8f)
        Data (in netCDF format) can be downloaded here: https://www.gebco.net/data_and_products/gridded_bathymetry_data/

        Contributed by: https://github.com/emiliebyer
    """
    def __init__(self, expansion_factor: float=1.2, tile: str='n61.0_s59.0_w4.0_e6.0', folder: str='/lustre/storeB/project/fou/om/WW3/bathy/gebco2021'):
        self.source=f'{folder}/gebco_2022_{tile}.nc'
        self.expansion_factor = expansion_factor
        return

    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float):
        # Area is expanded a bit to not get in trouble in the meshing stage
        # when we interpoolate or filter
        lon0, lon1, lat0, lat1 = expand_area(lon_min, lon_max, lat_min, lat_max, self.expansion_factor)

        ds = xr.open_dataset(self.source).sel(lon=slice(lon0, lon1), lat=slice(lat0, lat1))

        elevation = ds.elevation.values.astype(float)

        # Negative valies and NaN's are land
        topo = -1*elevation

        topo_lon = ds.lon.values.astype(float)
        topo_lat = ds.lat.values.astype(float)

        return topo, topo_lon, topo_lat

    def __str__(self):
        return(f"Reading GEBCO2022 topography from {self.source}.")


class Merge(TopoReader):
    """Merges raw topography from several grids"""

    def __init__(self, list_of_grids=None):
        self.list_of_grids = copy(list_of_grids)
        return
    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float) -> Tuple:

        topo = np.array([])
        topo_lon = np.array([])
        topo_lat = np.array([])
        for grid in self.list_of_grids:
            topo = np.append(topo,grid.raw_topo())
            topo_lon = np.append(topo_lon,grid.raw_lon())
            topo_lat = np.append(topo_lat,grid.raw_lat())

        return topo, topo_lon, topo_lat

    def __str__(self):
        return("Merging data from several grids.")

class EMODNET2020(TopoReader):
    """Reads bathymetry from multiple EMODNET tiles in netcdf format.

    For reading several files at once, supply the 'tile' argument with a glob pattern, e.g. 'C*'.

    Contributed by: https://github.com/poplarShift
    """

    def __init__(self, tile: str='C5', expansion_factor: float=1.2, folder: str='/lustre/storeB/project/fou/om/WW3/bathy/emodnet2020') -> Tuple:
        self.source = f'{folder}/{tile}_202*.nc' # for 2020 and 2022 version
        self.expansion_factor = expansion_factor
        return

    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float):
        # Area is expanded a bit to not get in trouble in the meshing stage
        # when we interpoolate or filter
        lon0, lon1, lat0, lat1 = expand_area(lon_min, lon_max, lat_min, lat_max, self.expansion_factor)

        def _crop(ds):
            """
            EMODNET tiles overlap by two cells on each boundary.
            """
            return ds.isel(lon=slice(2, -2), lat=slice(2, -2))

        import dask
        with dask.config.set(**{'array.slicing.split_large_chunks': True}):
            with xr.open_mfdataset(self.source, preprocess=_crop) as ds:
                ds = ds.sel(lon=slice(lon0, lon1), lat=slice(lat0, lat1))
                topo = -1 * ds.elevation.values
                topo_lon = ds.lon.values
                topo_lat = ds.lat.values
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

class MshFile(TopoReader):
    """Reads topography data from msh-file"""

    def __init__(self, filename: str, expansion_factor: float=1.2):
        self.filename = copy(filename)
        self.expansion_factor = expansion_factor
        return

    def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float) -> Tuple:
        import meshio

        mesh = meshio.read(self.filename)

        topo_lon = mesh.points[:,0]
        topo_lat = mesh.points[:,1]
        topo = mesh.points[:,2]

        lon0, lon1, lat0, lat1 = expand_area(lon_min, lon_max, lat_min, lat_max, self.expansion_factor)
        mask1=np.logical_and(topo_lon>=lon0, topo_lon<=lon1)
        mask2=np.logical_and(topo_lat>=lat0, topo_lat<=lat1)
        mask = np.logical_and(mask1, mask2)
        topo_lon = topo_lon[mask]
        topo_lat = topo_lat[mask]
        topo=topo[mask]

        return topo, topo_lon, topo_lat

    def __str__(self):
        return(f"Reading topography from {self.filename}.")
