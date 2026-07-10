from __future__ import annotations

import xarray as xr
import numpy as np
from typing import TYPE_CHECKING
import json
from urllib.request import urlopen
import re

if TYPE_CHECKING:
    from dnora.grid import Grid, TriGrid
from dnora.read.abstract_readers import DataReader
from dnora.utils.io import get_url
from dnora import utils
from dnora import msg
from typing import Union
import os
from dnora.type_manager.data_sources import DataSource
import dask
from pathlib import Path
import meshio

from dnora.type_manager.dnora_types import DnoraDataType
from .emodnet_functions import find_tile, get_covering_tiles, download_tile
from . import kartverket_functions as kartverket50m_functions
from . import gebco_functions
import dask.dataframe as dd


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
        lon, lat = utils.grid.expand_area(
            grid.edges("lon"), grid.edges("lat"), expansion_factor
        )
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

    Original version contributed by: https://github.com/emiliebyer
    """

    def default_data_source(self) -> DataSource:
        return DataSource.LOCAL

    def _folder(self, folder: str):
        return get_url(folder, "KartverketNo50m")

    @staticmethod
    def _get_files(folder, tiles):
        fn = []
        for tile in tiles:
            fn.append(f"{tile}_grid50_utm33.xyz")

        return get_url(folder, fn, get_list=True)

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid: Union[Grid, TriGrid],
        source: DataSource,
        folder: str,
        expansion_factor: float = 1.2,
        **kwargs,
    ) -> tuple:
        # Area is expanded a bit to not get in trouble in the meshing stage
        # when we interpoolate or filter


        folder = self._folder(folder)

        #self.source = get_url(folder, f"{tile}_grid50_utm{zone_number}.xyz")
        # grid.utm.set((zone_number, "W"))
        # Determine tiles
        print(f"Expansion factor: {expansion_factor}")
        lon, lat = utils.grid.expand_area(
            grid.edges("lon"),
            grid.edges("lat"),
            expansion_factor=expansion_factor,
        )
        x, y = utils.grid.expand_area(
            grid.edges("x", utm=(33,'W')),
            grid.edges("y", utm=(33,'W')),
            expansion_factor=expansion_factor,
            cartesian=True
        )
        tiles = kartverket50m_functions.find_tiles(lon, lat)
        
        msg.plain(f"Using tiles: {tiles}")
        files = self._get_files(folder, tiles)
        # Check if tiles exist locally
        tiles_to_download = []
        for file in files:
            if not os.path.isfile(file):
                tiles_to_download.append(Path(file).name[0:5])
        if tiles_to_download:
            msg.plain(f"Downloading tiles: {tiles_to_download}")
            if not os.path.exists(folder):
                if not os.path.isdir(folder):
                    msg.plain(f"Creating folder {folder}")
                os.makedirs(folder)
            email = kartverket50m_functions.get_geonorge_email()
            for tile in tiles_to_download:
                msg.blank()
                msg.process(f"Downloading tile {tile}")
                kartverket50m_functions.download_tile(tile=tile, folder=folder, email=email)

        
        df_list = []

        msg.blank()
        for file in files:
            msg.from_file(file)

            # Read file in chunks using Dask DataFrame
            df = dd.read_csv(file, sep=" ", header=None, dtype="float32")

            # Assign meaningful column names
            df.columns = ["x", "y", "z"]

            # Apply filtering based on x and y ranges
            df_filtered = df[
                (df["x"] >= x[0]) & (df["x"] <= x[1]) & 
                (df["y"] >= y[0]) & (df["y"] <= y[1])
            ]

            # Append the filtered DataFrame to the list
            df_list.append(df_filtered)

        msg.process("Merging tiles...")

        # Merge all filtered DataFrames into one
        if df_list:
            # Use Dask to concatenate all filtered DataFrames
            merged_df = dd.concat(df_list)

            # Compute the result (convert Dask DataFrame to Pandas DataFrame)
            
            merged_df = merged_df.compute()

            # Split into x, y, z arrays if needed
            
            x_all = merged_df["x"].to_numpy()
            y_all = merged_df["y"].to_numpy()
            z_all = merged_df["z"].to_numpy()

        else:
            # Empty arrays in case no data exists
            x_all = np.empty((0,))
            y_all = np.empty((0,))
            z_all = np.empty((0,))       
        
            
        coord_dict = {"x": x_all, "y": y_all}
        data_dict = {"topo": z_all, "zone_number": 33, "zone_letter": "W"}
        meta_dict = {"source": "Kartverket50m", 'url': 'https://kartkatalog.geonorge.no/metadata/dybdedata-terrengmodeller-50-meters-grid-landsdekkende/bbd687d0-d34f-4d95-9e60-27e330e0f76e'}

        return coord_dict, data_dict, meta_dict

    def __str__(self):
        return f"Reading Kartverket50m topography"


class GEBCO(DataReader):
    """Reads the GEBCO bathymetry
    """

    def default_data_source(self) -> DataSource:
        return DataSource.LOCAL

    def _folder(self, folder: str, year: int):
        return get_url(folder, f"GEBCO/{year}")

    @staticmethod
    def _get_files(folder, tiles, year: int):
        fn = []
        for tile in tiles:
            fn.append(f"{tile}_gebco_{year}.nc")

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
        year: int = 2026,
        **kwargs,
    ) -> tuple:

        folder = self._folder(folder, year)
        
        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        lons, lats = gebco_functions.get_covering_tiles(grid, expansion_factor)
        tiles = gebco_functions.get_tile_names(lons, lats)
        
        
        files = self._get_files(folder, tiles, year)
        # Check if tiles exist locally
        tiles_to_download = []
        for n, file in enumerate(files):
            if not os.path.isfile(file):
                tiles_to_download.append((Path(file), lons[n], lats[n], tiles[n]))
        
        if tiles_to_download:
            url = f"https://dap.ceda.ac.uk/thredds/dodsC/bodc/gebco/global/gebco_{year}/ice_surface_elevation/netcdf/GEBCO_{year}.nc"
            ds = xr.open_dataset(url)
            msg.blank()
            msg.process(f"Downloading {len(tiles_to_download)} tiles...")
            if not os.path.exists(folder):
                if not os.path.isdir(folder):
                    msg.plain(f"Creating folder {folder}")
                os.makedirs(folder)
            msg.from_file(url)
            for tile in tiles_to_download:
                msg.to_file(tile[0])
                ds.sel(lon=slice(*tile[1]), lat=slice(*tile[2])).to_netcdf(tile[0])

        msg.process('Reading local files...')
        lon, lat = utils.grid.expand_area(
            grid.edges("lon"),
            grid.edges("lat"),
            expansion_factor=expansion_factor,
        )
        with dask.config.set(**{"array.slicing.split_large_chunks": True}):
            with xr.open_mfdataset(files) as ds:
                msg.from_multifile(files)
                ds = ds.sel(lon=slice(lon[0], lon[1]), lat=slice(lat[0], lat[1]))

                coord_dict = {key: ds.get(key) for key in ["lon", "lat", "x", "y"]}
                data_dict = {"topo": -ds.get("elevation").values}
                meta_dict = ds.attrs

                return coord_dict, data_dict, meta_dict
            
    def __str__(self):
        return f"Reading GEBCO bathymetry."


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
        **kwargs,
    ) -> tuple:
        self.filename = filename
        mesh = meshio.read(filename)

        topo_x = mesh.points[:, 0]
        topo_y = mesh.points[:, 1]
        topo = mesh.points[:, 2]

        if zone_number is None or zone_letter is None:
            xedges, yedges = utils.grid.expand_area(
                grid.edges("lon"), grid.edges("lat"), expansion_factor
            )
        else:
            xedges, yedges = utils.grid.expand_area(
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
