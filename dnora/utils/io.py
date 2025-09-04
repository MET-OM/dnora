import os
import re
from pathlib import Path
import pandas as pd
import numpy as np
from typing import Union
def set_nested_value(dictionary, keys, value):
    for key in keys[:-1]:  # Iterate through all keys except the last one
        dictionary = dictionary[key]
    dictionary[keys[-1]] = value  # Set the value for the last key


def read_ww3_nml(filename: str) -> dict:
    """Reads a WAVEWATCH III namelist file (FORTRAN style namelist)"""
    with open(filename, "r") as file:
        nml_dict = {}
        line = "!"
        while line:
            line = file.readline()
            if len(line) < 1:
                continue
            if line[0] == "&":
                main_key = line.replace("\n", "").replace("&", "")
                nml_dict[main_key] = {}
                line = file.readline()
                while line[0] != "/":
                    keys, value = line.replace(" ", "").split("=")
                    keys = keys.split("%")
                    dictionary = nml_dict[main_key]
                    for key in keys[:-1]:
                        if dictionary.get(key) is None:
                            dictionary[key] = {}
                        dictionary = dictionary[key]

                    value = value.replace("\n", "")
                    set_nested_value(nml_dict[main_key], keys, value)
                    line = file.readline()

    return nml_dict


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
