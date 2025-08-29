import numpy as np
import pandas as pd

from pathlib import Path

import re
import os


#





def check_if_file(filename: str, halt=False) -> bool:
    """Checks if a file exists and halts with an error if halt=True"""
    exists = os.path.isfile(filename)
    if not exists and halt:
        raise FileNotFoundError(f"DNORA cannot find file {filename}")
    return exists


def check_if_folder(folder: str, create: bool = True) -> bool:
    """Creates a folder if it does not exist, and returns True if it
    already existed."""

    if folder == "":
        existed = True
    else:
        existed = os.path.isdir(folder)

    if not existed:
        os.mkdir(folder)

    return existed


# -----------------------------------------------------------------------------
# MISC STAND ALONE FUNCTIONS
# # -----------------------------------------------------------------------------
# def read_ww3_info(
#     filename,
# ) -> tuple[float, float, float, float, float, float, int, int]:
#     """Read grid specification from the GridName_info.txt file"""
#     with open(filename, "r") as f:
#         lines = f.readlines()

#     for n in range(len(lines)):
#         line = lines[n].split()

#         if len(line):
#             if line[0] == "lon:":
#                 lon_min = float(line[1])
#                 lon_max = float(line[3][0:-1])
#                 lat_min = float(line[5])
#                 lat_max = float(line[7])
#             elif line[0] == "dlon,":
#                 dlon = float(line[3][0:-1])
#                 dlat = float(line[4])
#             elif line[0] == "nx,":
#                 nx = int(line[3])
#                 ny = int(line[5])
#     return lon_min, lon_max, lat_min, lat_max, dlon, dlat, nx, ny


def u_v_from_speed_dir(ws, wdir) -> tuple[float, float]:
    """Converts wind speed and direction (from) to u and v components."""

    # see http://tornado.sfsu.edu/geosciences/classes/m430/Wind/WindDirection.html
    u = -ws * (np.sin(np.deg2rad(wdir)))
    v = -ws * (np.cos(np.deg2rad(wdir)))

    return u, v


# def speed_dir_from_u_v(u, v) -> tuple[float, float]:
#     """Convert component to speed and direction (from)"""
#     ws = (u**2 + v**2) ** 0.5
#     wdir = np.mod((90 - np.rad2deg(np.arctan2(v, u))) + 180, 360)
#     return ws, wdir


def get_url(
    folder: str,
    filename: str,
    time_stamp: pd.DatetimeIndex | str = None,
    get_list: bool = False,
) -> str:
    if not isinstance(filename, list):
        filename = [filename]
    if time_stamp is not None:
        time_stamp = pd.to_datetime(time_stamp)
        filename = [time_stamp.strftime(fn) for fn in filename]
        folder = time_stamp.strftime(folder)
    url = []
    for fn in filename:
        url_temp = Path(folder).joinpath(fn)
        url_temp = re.sub(f"https:/", "https://", str(url_temp), 1)
        url_temp = re.sub(f"http:/", "http://", str(url_temp), 1)
        url_temp = re.sub(f"ftp:/", "ftp://", str(url_temp), 1)
        if time_stamp is not None:
            for floor_hour in range(1, 24):
                hfloor = int(np.floor(time_stamp.hour / floor_hour) * floor_hour)
                url_temp = re.sub(f"\[{floor_hour}\]", f"{hfloor:02.0f}", str(url_temp))

        url.append(url_temp)
    if len(url) == 1 and not get_list:
        return os.path.expanduser(url[0])
    else:
        return [os.path.expanduser(u) for u in url]


def set_metaparameters_in_object(obj, metaparameter_dict, data_dict):
    for key, value in data_dict.items():
        metaparameter = metaparameter_dict.get(
            key
        )  # Check if metaparameter provided by reader

        if metaparameter is None:
            # DNORA object usually has specified the metaparameters
            if hasattr(obj, "meta_dict"):
                metaparameter = obj.meta_dict.get(key)

        if metaparameter is not None:
            obj.set_metadata(metaparameter.meta_dict(), name=key)

    return obj
