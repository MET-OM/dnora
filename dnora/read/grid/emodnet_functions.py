import numpy as np
from urllib.request import urlretrieve
import zipfile
import os.path

import progressbar
from .tiling_functions import read_limits

pbar = None


def show_progress(block_num, block_size, total_size):
    global pbar
    if pbar is None:
        pbar = progressbar.ProgressBar(maxval=total_size)
        pbar.start()

    downloaded = block_num * block_size
    if downloaded < total_size:
        pbar.update(downloaded)
    else:
        pbar.finish()
        pbar = None




def find_tile(lon, lat):
    tiles, lons, lats = read_limits("emodnet_dtm_coords.txt")
    west_of = lon <= lons[:, 1]
    east_of = lon >= lons[:, 0]
    north_of = lat >= lats[:, 0]
    south_of = lat <= lats[:, 1]
    lon_ok = np.logical_and(east_of, west_of)
    lat_ok = np.logical_and(north_of, south_of)
    if not np.any(np.logical_and(lon_ok, lat_ok)):
        return []

    ind = np.where(np.logical_and(lon_ok, lat_ok))[0][0]

    # print(
    #     f"Found {lon}, {lat} inside tile {tiles[ind]}: {lons[ind,0]}-{lons[ind,1]}, {lats[ind,0]}-{lats[ind,1]}"
    # )
    tile = tiles[ind]

    return tile


def get_covering_tiles(tile_nw, tile_se):
    list_of_letters = ["B", "C", "D", "E", "F", "G", "H"]
    bad_tiles = {"D8", "G4", "H4", "H5", "H6", "H7"}
    numbers = np.arange(int(tile_nw[1]), int(tile_se[1]) + 1)
    ind0 = list_of_letters.index(tile_nw[0])
    ind1 = list_of_letters.index(tile_se[0])
    letters = list_of_letters[ind0 : ind1 + 1]
    tiles = []
    for l in letters:
        for n in numbers:
            tiles.append(f"{l}{n}")
    return list(set(tiles) - bad_tiles)


def download_tile(tile: str, year: int, folder: str):
    versions = {2022: "v11", 2020: "v10"}
    url = f"https://downloads.emodnet-bathymetry.eu/{versions[year]}/{tile}_{year:.0f}.nc.zip"
    filename = f"{folder}/{tile}_{year:.0f}.nc.zip"

    print(f"Downloading {url}...")
    urlretrieve(url, filename, show_progress)

    if os.path.isfile(filename):
        print(f"Extreacting {filename}...")
        with zipfile.ZipFile(filename, "r") as zip_ref:
            zip_ref.extractall(folder)
        os.remove(filename)
