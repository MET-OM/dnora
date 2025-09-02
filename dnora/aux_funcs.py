import numpy as np
import pandas as pd

from pathlib import Path

import re
import os
from typing import Union



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


def get_url(
    folder: str,
    filename: str,
    time_stamp: Union[pd.DatetimeIndex, str] = None,
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
        # If we are on a Windows machine, then / will be replaced by \
        # If we have an url, then we still want /
        if 'http' in str(url_temp) or 'ftp' in str(url_temp):
            url_temp = url_temp.as_posix()
        else:
            url_temp =str(url_temp)
        
        url_temp = re.sub("https:/", "https://", url_temp, 1)
        url_temp = re.sub("http:/", "http://", url_temp, 1)
        url_temp = re.sub("ftp:/", "ftp://", url_temp, 1)
        
        if time_stamp is not None:
            for floor_hour in range(1, 24):
                hfloor = int(np.floor(time_stamp.hour / floor_hour) * floor_hour)
                url_temp = re.sub(f"\[{floor_hour}\]", f"{hfloor:02.0f}", url_temp)

        url.append(url_temp)
    if len(url) == 1 and not get_list:
        return os.path.expanduser(url[0])
    else:
        return [os.path.expanduser(u) for u in url]
