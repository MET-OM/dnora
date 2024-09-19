import xarray as xr
import pandas as pd
import numpy as np
from .file_structure import FileStructure
import os, glob
from dnora.aux_funcs import get_url
from dnora import msg
from dnora.type_manager.dnora_types import DnoraDataType

from typing import Callable
import geo_parameters as gp


def read_first_ds(
    folder: str,
    filename: str,
    start_time: str,
    file_structure: FileStructure = None,
) -> xr.Dataset:
    if file_structure is not None:
        _, _, file_times = file_structure.create_time_stamps(start_time, start_time)
        file_time = file_times[0]
    else:
        file_time = start_time
    url = get_url(folder, filename, file_time)
    ds = xr.open_dataset(url).isel(time=[0])

    return ds


def read_one_ds(
    start_time: str,
    end_time: str,
    file_times: list[str],
    urls: list[str],
    hours_per_file: int,
    n: int,
    expected_shape: tuple[int],
    ds_creator_function: callable = None,
) -> tuple[xr.Dataset, str]:
    """This functions reads one Dataset and crops it.
    If the expected file is not found, it goes to the previous one nad reads data with a longer lead time.

    This function can be used for e.g. forecast data where we have overlapping data in time in different files

    ds_creator_fuction takes arguments (start_time, end_time, url) and returns an xr.Dataset
    The ds_creator_function might use normal xarray or e.g. fimex and is therefore injected as a callable
    """

    ct = 0
    keep_trying = True
    try_next_file = False
    while keep_trying:

        try:
            url = urls[n - ct]
            ds = ds_creator_function(start_time, end_time, url)
            if not file_is_consistent(ds, expected_shape):
                msg.plain(f"SKIPPING, file inconsistent: {url}")
                try_next_file = True
            else:
                try_next_file = False
                keep_trying = False

        except (OSError, RuntimeError):  # xr gives OSError, fimex gives RuntimeError
            msg.plain(f"SKIPPING, file not found: {url}")
            try_next_file = True

        if try_next_file:
            if data_left_to_try_with(hours_per_file, n, ct, file_times, end_time):
                ct += 1
            else:
                msg.plain(f"SKIPPING, no data left to try with!")
                ds = None
                keep_trying = False

    return ds, url


def get_constant_url(folder, filename, file_times) -> list[str]:
    """Applies the same folder and filename to all file_times to get url.
    folder and file_name can contain %Y etc. that will be replaced"""
    return [get_url(folder, filename, file_time) for file_time in file_times]


def read_ds_list(
    start_times: pd.DatetimeIndex,
    end_times: pd.DatetimeIndex,
    file_times: pd.DatetimeIndex,
    folder: str,
    filename: str,
    ds_creator_function: Callable,
    url_function: Callable = get_constant_url,
    hours_per_file: int = None,
) -> list[xr.Dataset]:
    """Reads a list of xr.Datasets using the time stamps and folder/filename given.

    If one file is missing, the function tries to patch if the files have overlap (if hours_per_file is given).
    the ds_creator function takes arguments (start_time, end_time, url) and returns an xr.Dataset.
    """

    urls = url_function(folder, filename, file_times)
    ds_list = []
    expected_shape = None
    for n in range(len(file_times)):
        msg.plain(f"Reading data for: {start_times[n]}-{end_times[n]}")

        ds, url = read_one_ds(
            start_times[n],
            end_times[n],
            file_times,
            urls,
            hours_per_file,
            n,
            expected_shape,
            ds_creator_function,
        )

        if ds is not None:
            msg.from_file(url)
            if not ds_list:
                keys = list(ds.sizes.keys())
                expected_shape = tuple([ds[var].size for var in keys if var != "time"])
            ds_list.append(ds)

    return ds_list


def create_dicts(ds, var_mapping: dict):
    coord_dict = {
        "time": ds[var_mapping["time"]].values.squeeze(),
        "lon": ds[var_mapping["lon"]].values,
        "lat": ds[var_mapping["lat"]].values,
        "freq": ds[var_mapping["freq"]].values.squeeze(),
        "dirs": ds[var_mapping["dirs"]].values.squeeze(),
    }
    data_dict = {"spec": ds[var_mapping["spec"]].data}

    meta_dict = ds.attrs
    return coord_dict, data_dict, meta_dict


def create_coord_dict(wanted_coords: list, ds: xr.Dataset, alias_mapping: dict) -> dict:
    """Creates a dict of wanted coordinates from xr.Dataset and identifies which variables in the ds has been used for coordinates"""
    coord_dict = {}
    ds_coord_strings = []
    for coord in wanted_coords:
        if not hasattr(ds, coord):
            ds_coord = alias_mapping[coord]
        else:
            ds_coord = coord
        ds_coord_strings.append(ds_coord)
        coord_dict[coord] = ds[ds_coord].values.squeeze()
    return coord_dict, ds_coord_strings


def create_data_dict(wanted_vars: list, ds: xr.Dataset, alias_mapping: dict) -> dict:
    data_dict = {}
    for var in wanted_vars:
        param = alias_mapping.get(var)
        if param is None:
            standard_name = ds[var].attrs.get("standard_name")
            param = gp.get(standard_name)

        if param is None:
            raise ValueError(f"Could not find dnora name for variable {var}!")

        data_dict[param] = ds[var].data
    return data_dict


def data_left_to_try_with(hours_per_file, n, ct, file_times, end_time) -> bool:
    """Checks if we can go back one file and still have data covering the entire period.

    E.g. Each files conains 72 hours but we want to read in 6 hour chunks.
    If one ifle is missing we can use hours 7-12 in the previous file etc."""
    # Structure is not defined. E.g. if we have montly files and know we don't have any overlap
    if hours_per_file is None:
        return False

    if n - ct <= 0:
        return False

    if pd.Timestamp(end_time) - pd.Timestamp(file_times[n - ct - 1]) > pd.Timedelta(
        hours_per_file, "hours"
    ):
        return False

    return True


def file_is_consistent(
    ds: xr.Dataset,
    expected_shape: tuple[int],
) -> bool:
    """Checks if dataset is consistent with the expected points"""

    # always trust the first file that is read
    if expected_shape is None:
        return True

    keys = list(ds.sizes.keys())
    given_shape = tuple([ds[var].size for var in keys if var != "time"])
    if given_shape == expected_shape:
        return True
    else:
        return False


def setup_temp_dir(data_type: DnoraDataType, reader_name: str) -> None:
    """Sets up a temporery directory for fimex files and cleans out possible old files"""
    temp_folder = f"dnora_{data_type.name.lower()}_temp"
    if not os.path.isdir(temp_folder):
        os.mkdir(temp_folder)
        print("Creating folder %s..." % temp_folder)

    msg.plain("Removing old files from temporary folder...")
    for f in glob.glob(f"dnora_{data_type.name.lower()}_temp/{reader_name}*.nc"):
        os.remove(f)
