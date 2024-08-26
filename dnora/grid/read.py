from __future__ import annotations

import xarray as xr
from abc import ABC, abstractmethod
import numpy as np
from typing import TYPE_CHECKING
import json
from urllib.request import urlopen
import re

if TYPE_CHECKING:
    from dnora.grid import Grid, TriGrid

from dnora.aux_funcs import expand_area, get_url
from dnora import msg
from typing import Union
import os
from dnora.dnora_type_manager.data_sources import DataSource
import dask
from pathlib import Path
import meshio
from dnora.readers.abstract_readers import DataReader
from dnora.dnora_type_manager.dnora_types import DnoraDataType

# from dnora.defaults import read_environment_variable
from .emodnet_functions import find_tile, get_covering_tiles, download_tile

# class TopoReader(ABC):
#     """Abstract class for reading the bathymetry."""

#     @abstractmethod
#     def __call__(self, grid: Union[Grid, TriGrid], source: DataSource, **kwargs):
#         """Reads the bathymetrical information from a source and returns the data.

#         !!!! DEPTH VALUES ARE POSITIVE AND EVERYTHING ELSE (including 0)
#         IS INTERPRETED AS LAND.

#         LON, LAT IN WGS84

#         X, Y IN UTM !!!


#         Two format options are available:

#         1) Matrix, where lon, lat or x, y are vectors and topo is the
#         corresponding matrix.

#         The dimensions and orientation of the bathymetry array that is returned to
#         the object should be:

#         rows = latitude and colums = longitude (i.e.) shape = (nr_lat, nr_lon).

#         North = [-1,:]
#         South = [0,:]
#         East = [:,-1]
#         West = [:,0]

#         2) Xyz, where topo, lon, lat, x, y are all vectors with the same length.

#         topo_x, topo_y, zone_number, zone_letter None in spherical data
#         topo_lon, topo_lat None in cartesian data

#         ----

#         This method is called from within the Grid-object

#         return (
#             topo,
#             coord_dict,
#             zone_number,
#             zone_letter,
#             metadata,
#         )
#         """

#     @abstractmethod
#     def __str__(self):
#         """Describes what topography is read and from where.

#         This is called by the Grid-object to provide output to the user.
#         """
#         pass

#     def default_data_source(self) -> DataSource:
#         return DataSource.LOCAL

#     def name(self) -> str:
#         return type(self).__name__


# class ConstantTopo(DataReader):
#     """Creates an empty topography. Called when setting initial spacing."""

#     def __init__(self, depth: float = 999.0):
#         self.depth = float(depth)

#     def __call__(
#         self, grid: Union[Grid, TriGrid], source: DataSource, folder: str, **kwargs
#     ):
#         """Creates a trivial topography with all water points."""
#         topo = np.full(grid.size(), self.depth)
#         zone_number, zone_letter = grid.utm()
#         coord_dict = {
#             "lon": grid.lon(strict=True),
#             "lat": grid.lat(strict=True),
#             "x": grid.x(strict=True),
#             "y": grid.y(strict=True),
#         }
#         data_dict = {
#             "topo": topo,
#             "zone_number": zone_number,
#             "zone_letter": zone_letter,
#         }
#         meta_dict = {"source": type(self).__name__}
#         metaparameter_dict = {}

#         return coord_dict, data_dict, meta_dict, metaparameter_dict

#     def __str__(self):
#         return f"Creating an constant topography with depth values {self.depth}."


class EMODNET(DataReader):
    """Reads bathymetry from multiple EMODNET tiles in netcdf format.

    Contributed by: https://github.com/poplarShift
    """

    def default_data_source(self) -> DataSource:
        return DataSource.LOCAL

    @staticmethod
    def _get_files(folder, tiles, year):
        fn = []
        for tile in tiles:
            fn.append(f"{tile}_{year}.nc")

        return get_url(folder, fn, get_list=True)

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid: Union[Grid, TriGrid],
        start_time,
        end_time,
        source: DataSource,
        expansion_factor: float = 1.2,
        folder: str = None,
        year: int = 2022,
        **kwargs,
    ) -> tuple:

        folder = get_url(folder, f"EMODNET/{year:.0f}")
        # Area is expanded a bit to not get in trouble in the meshing stage
        # when we interpoolate or filter

        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)
        msg.plain(
            f"Downloading bathymetry for: {lon[0]:10.7f}-{lon[1]:10.7f}, {lat[0]:10.7f}-{lat[1]:10.7f}."
        )

        # Identefy tiles:
        tile_nw = find_tile(lon[0], lat[1])
        tile_se = find_tile(lon[1], lat[0])

        if not tile_nw or not tile_se:
            msg.warning(f"Area not coverd by EMODNET!!!")
            return
        self.tiles = get_covering_tiles(tile_se=tile_se, tile_nw=tile_nw)

        msg.plain(f"Using tiles: {self.tiles}")
        self.files = self._get_files(folder, self.tiles, year)
        # Check if tiles exist locally
        tiles_to_download = []
        for file in self.files:
            if not os.path.isfile(file):
                tiles_to_download.append(Path(file).name[0:2])

        if tiles_to_download:
            if not os.path.exists(folder):
                if not os.path.isdir(folder):
                    msg.plain(f"Creating folder {folder}")
                os.makedirs(folder)
            for tile in tiles_to_download:
                msg.from_file(tile)
                download_tile(tile=tile, year=year, folder=folder)

        def _crop(ds):
            """
            EMODNET tiles overlap by two cells on each boundary.
            """
            return ds.isel(lon=slice(2, -2), lat=slice(2, -2))

        with dask.config.set(**{"array.slicing.split_large_chunks": True}):
            with xr.open_mfdataset(self.files, preprocess=_crop) as ds:
                msg.from_file(self.files)
                ds = ds.sel(lon=slice(lon[0], lon[1]), lat=slice(lat[0], lat[1]))
                topo = -1 * ds.elevation.values
                coord_dict = {"lon": ds.lon.values, "lat": ds.lat.values}

                data_dict = {"topo": topo}
                meta_dict = {"source": f"EMODNET{year:.0f}"}

                return coord_dict, data_dict, meta_dict

    def __str__(self):
        return f"Reading EMODNET topography from {self.files}."


class KartverketNo50m(DataReader):
    """Reads data from Kartverket bathymetry.

    High resolution bathymetry dataset for the whole Norwegian Coast.
    Can be found at:
    https://kartkatalog.geonorge.no/metadata/dybdedata-terrengmodeller-50-meters-grid-landsdekkende/bbd687d0-d34f-4d95-9e60-27e330e0f76e

    For reading several files at once, supply the 'tile' argument with a glob pattern, e.g. 'B*'.

    Contributed by: https://github.com/emiliebyer
    """

    def default_data_source(self) -> DataSource:
        return DataSource.LOCAL

    def _folder(self, folder: str):
        return get_url(folder, "KartverketNo50m")

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid: Union[Grid, TriGrid],
        source: DataSource,
        folder: str,
        expansion_factor: float = 1.2,
        zone_number: int = 33,
        tile: str = "B1008",
        **kwargs,
    ) -> tuple:
        # Area is expanded a bit to not get in trouble in the meshing stage
        # when we interpoolate or filter

        folder = self._folder(folder)

        self.source = get_url(folder, f"{tile}_grid50_utm{zone_number}.xyz")
        x, y = expand_area(grid.edges("x"), grid.edges("y"), expansion_factor)

        print(f"Expansion factor: {expansion_factor}")

        import dask.dataframe as dd

        df = dd.read_csv(self.source, sep=" ", header=None)
        df.columns = ["x", "y", "z"]
        topo_x = np.array(df["x"].astype(float))
        topo_y = np.array(df["y"].astype(float))
        z = np.array(df["z"].astype(float))

        mask_x = np.logical_and(x[0] <= topo_x, topo_x <= x[1])
        mask_y = np.logical_and(y[0] <= topo_y, topo_y <= y[1])
        mask = np.logical_and(mask_x, mask_y)

        topo_x = topo_x[mask]
        topo_y = topo_y[mask]
        topo = z[mask]

        coord_dict = {"x": topo_x, "y": topo_y}
        data_dict = {"topo": topo, "zone_number": zone_number, "zone_letter": "W"}
        meta_dict = {"source": "Kartverket50m"}

        return coord_dict, data_dict, meta_dict

    def __str__(self):
        return f"Reading Kartverket topography from {self.source}."


class GEBCO(DataReader):

    def default_data_source(self) -> DataSource:
        return DataSource.LOCAL

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid: Union[Grid, TriGrid],
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        expansion_factor: float = 1.2,
        tile: str = "n61.0_s59.0_w4.0_e6.0",
        year: int = 2023,
        **kwargs,
    ) -> tuple:

        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)
        # url = f'https://api.odb.ntu.edu.tw/gebco?mode=zonly&sample=1&jsonsrc={"type":"Polygon","coordinates":[[[{lon[0]},[{lat[1]}]],[{lon[1]},{lat[1]}],[{lon[1]},{lon[0]}],[{lon[0]},{lat[0]}],[{lon[0]},{lat[1]}]]]}'
        url = 'https://api.odb.ntu.edu.tw/gebco?mode=zonly&sample=1&jsonsrc={"type":"Polygon","coordinates":[[[LON0,LAT1],[LON1,LAT1],[LON1,LAT0],[LON0,LAT0],[LON0,LAT1]]]}'
        url = re.sub("LON0", f"{lon[0]:.2f}", url)
        url = re.sub("LON1", f"{lon[1]:.2f}", url)
        url = re.sub("LAT0", f"{lat[0]:.2f}", url)
        url = re.sub("LAT1", f"{lat[1]:.2f}", url)
        msg.from_file(url)
        response = urlopen(url)

        data = json.loads(response.read())

        # Negative valies and NaN's are land
        topo = -1 * np.array(data["z"])

        topo_lon = np.array(data["longitude"])
        topo_lat = np.array(data["latitude"])

        coord_dict = {"lon": topo_lon, "lat": topo_lat}
        data_dict = {"topo": topo}
        meta_dict = {"source": f"GEBCO2023", "through": "https://api.odb.ntu.edu.tw/"}

        return coord_dict, data_dict, meta_dict


class GEBCOold(DataReader):
    """Reads the GEBCO_2021 gridded bathymetric data set.
    A global terrain model for ocean and land,
    providing elevation data, in meters, on a 15 arc-second interval grid.
    Reference: GEBCO Compilation Group (2021) GEBCO 2021 Grid (doi:10.5285/c6612cbe-50b3-0cff- e053-6c86abc09f8f)
    Data (in netCDF format) can be downloaded here: https://www.gebco.net/data_and_products/gridded_bathymetry_data/

    Contributed by: https://github.com/emiliebyer
    """

    def _folder(self, folder: str, year: int):
        return get_url(folder, f"GEBCO/{year}")

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid: Union[Grid, TriGrid],
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        expansion_factor: float = 1.2,
        tile: str = "n61.0_s59.0_w4.0_e6.0",
        year: int = 2021,
        **kwargs,
    ) -> tuple:

        folder = self._folder(folder, year)
        self.source = get_url(folder, f"gebco_{year}_{tile}.nc")
        # Area is expanded a bit to not get in trouble in the meshing stage
        # when we interpoolate or filter

        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

        ds = xr.open_dataset(self.source).sel(
            lon=slice(lon[0], lon[1]), lat=slice(lat[0], lat[1])
        )

        elevation = ds.elevation.values.astype(float)

        # Negative valies and NaN's are land
        topo = -1 * elevation

        topo_lon = ds.lon.values.astype(float)
        topo_lat = ds.lat.values.astype(float)

        coord_dict = {"lon": topo_lon, "lat": topo_lat}
        data_dict = {"topo": topo}
        meta_dict = {"source": f"GEBCO{year}"}
        metaparameter_dict = {}

        return coord_dict, data_dict, meta_dict, metaparameter_dict

    def __str__(self):
        return f"Reading GEBCO topography from {self.source}."


# class Merge(TopoReader):
#     """Merges raw topography from several grids"""
#
#     def __init__(self, list_of_grids=None):
#         self.list_of_grids = copy(list_of_grids)
#         return
#     def __call__(self, lon_min: float, lon_max: float, lat_min: float, lat_max: float) -> Tuple:
#
#         topo = np.array([])
#         topo_lon = np.array([])
#         topo_lat = np.array([])
#         for grid in self.list_of_grids:
#             topo = np.append(topo,grid.raw_topo())
#             topo_lon = np.append(topo_lon,grid.raw_lon())
#             topo_lat = np.append(topo_lat,grid.raw_lat())
#
#         return topo, topo_lon, topo_lat
#
#     def __str__(self):
#         return("Merging data from several grids.")


# class ForceFeed(TopoReader):
#     """Simply passes on the data it was fed upon initialization"""

#     def __init__(self, topo):
#         self.topo = topo
#         return

#     def __call__(
#         self, grid: Union[Grid, TriGrid], source: DataSource, **kwargs
#     ) -> tuple:
#         # Just use the values it was forcefed on initialization
#         topo = self.topo

#         return
#             topo,
#             grid.lon(strict=True),
#             grid.lat(strict=True),
#             grid.x(strict=True),
#             grid.y(strict=True),
#             None,
#             None,
#             {"source": "ForceFeed"},
#         )

#     def __str__(self):
#         return "Passing on the topography I was initialized with."


class MshFile(DataReader):
    """Reads topography data from msh-file"""

    def default_data_source(self) -> DataSource:
        return DataSource.LOCAL

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid: Union[Grid, TriGrid],
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        filename: str = None,
        expansion_factor: float = 1.2,
        zone_number: int = None,
        zone_letter: str = None,
    ) -> tuple:
        self.filename = filename
        mesh = meshio.read(filename)

        topo_x = mesh.points[:, 0]
        topo_y = mesh.points[:, 1]
        topo = mesh.points[:, 2]

        if zone_number is None or zone_letter is None:
            xedges, yedges = expand_area(
                grid.edges("lon"), grid.edges("lat"), expansion_factor
            )
        else:
            xedges, yedges = expand_area(
                grid.edges("x"), grid.edges("y"), expansion_factor
            )

        mask1 = np.logical_and(topo_x >= xedges[0], topo_x <= xedges[1])
        mask2 = np.logical_and(topo_y >= yedges[0], topo_y <= yedges[1])
        mask = np.logical_and(mask1, mask2)
        topo_x = topo_x[mask]
        topo_y = topo_y[mask]
        topo = topo[mask]
        if zone_number is None or zone_letter is None:
            msg.plain(
                f"No utm-zone, e.g. zone_number=33, zone_letter='W', provided. Assuming data in {self.filename} is in lon-lat!"
            )
            coord_dict = {"lon": topo_x, "lat": topo_y}
        else:
            coord_dict = {"x": topo_x, "y": topo_y}

        data_dict = {
            "topo": topo,
            "zone_number": zone_number,
            "zone_letter": zone_letter,
        }
        meta_dict = {"source": self.filename}

        return coord_dict, data_dict, meta_dict

    def __str__(self):
        return f"Reading topography from {self.filename}."
