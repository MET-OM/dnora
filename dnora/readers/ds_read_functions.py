import xarray as xr
import pandas as pd
import numpy as np
from .file_structure import FileStructure

from dnora.aux_funcs import get_url
from dnora import msg


def read_first_ds(
    file_structure: FileStructure, start_time: str, folder: str, filename: str
) -> xr.Dataset:
    start_times, end_times, file_times = file_structure.create_time_stamps(
        start_time, start_time
    )
    url = get_url(folder, filename, file_times[0])
    ds = xr.open_dataset(url).isel(time=[0])

    return ds


def read_one_ds(
    start_time: str,
    end_time: str,
    file_times: list[str],
    urls: list[str],
    hours_per_file: int,
    n: int,
    expected_lon: np.ndarray,
    expected_lat: np.ndarray,
    lon_str: str,
    lat_str: str,
    ds_creator_function: callable,
) -> tuple[xr.Dataset, str]:
    """This functions reads one Dataset and crops it.
    If the expected file is not found, it goes to the previous one nad reads data with a longer lead time.

    This function can be used for e.g. forecast data where we have overlapping data in time in different files

    ds_creator_fuction takes arguments (start_time, end_time, url) and returns an xr.Dataset
    The ds_creator_function might use normal xarray or e.g. fimex and is therefore injected as a callable
    """

    file_time = file_times[n]
    ct = 1
    keep_trying = True
    try_next_file = False
    while keep_trying:

        try:
            url = urls[n]
            ds = ds_creator_function(start_time, end_time, url)
            if file_is_consistent(
                ds, expected_lon, expected_lat, url, lon_str, lat_str
            ):
                try_next_file = False
                keep_trying = False
            else:
                try_next_file = True
        except OSError:
            try_next_file = True

        if try_next_file:
            if data_left_to_try_with(hours_per_file, n, ct, file_times, end_times):
                file_time = file_times[n - ct]
                ct += 1
            else:
                ds = None
                keep_trying = False

    return ds, url


def read_ds_list(
    start_times: pd.DatetimeIndex,
    end_times: pd.DatetimeIndex,
    file_times: pd.DatetimeIndex,
    folder: str,
    filename: str,
    ds_creator_function: callable,
    hours_per_file: int = None,
    lon_str: str = "longitude",
    lat_str: str = "latitude",
) -> list[xr.Dataset]:

    urls = [get_url(folder, filename, file_time) for file_time in file_times]
    ds_list = []
    expected_lon, expected_lat = None, None
    for n in range(len(file_times)):
        msg.plain(f"Reading data for: {start_times[n]}-{end_times[n]}")

        ds, url = read_one_ds(
            start_times[n],
            end_times[n],
            file_times,
            urls,
            hours_per_file,
            n,
            expected_lon,
            expected_lat,
            lon_str,
            lat_str,
            ds_creator_function,
        )

        if ds is not None:
            msg.from_file(url)
            if not ds_list:
                expected_lon = ds[lon_str].values
                expected_lat = ds[lat_str].values
            ds_list.append(ds)

    return ds_list


def create_dicts(ds, var_mapping: dict):
    coord_dict = {
        "time": ds[var_mapping["time"]].values.squeeze(),
        "lon": ds[var_mapping["lon"]].values.squeeze(),
        "lat": ds[var_mapping["lat"]].values.squeeze(),
        "freq": ds[var_mapping["freq"]].values.squeeze(),
        "dirs": ds[var_mapping["dirs"]].values.squeeze(),
    }
    data_dict = {"spec": ds[var_mapping["spec"]].data}

    meta_dict = ds.attrs
    return coord_dict, data_dict, meta_dict


def data_left_to_try_with(hours_per_file, n, ct, file_times, end_times) -> bool:
    """Checks if we can go back one file and still have data covering the entire period.

    E.g. Each files conains 72 hours but we want to read in 6 hour chunks.
    If one ifle is missing we can use hours 7-12 in the previous file etc."""
    # Structure is not defined. E.g. if we have montly files and know we don't have any overlap
    if hours_per_file is None:
        return False

    if n - ct < 0:
        return False

    if pd.Timestamp(end_times[n]) - pd.Timestamp(file_times[n - ct]) > pd.Timedelta(
        hours_per_file, "hours"
    ):
        return False

    return True


def file_is_consistent(
    ds: xr.Dataset,
    expected_lon: np.ndarray,
    expected_lat: np.ndarray,
    url: str,
    lon_str: str,
    lat_str: str,
) -> bool:
    """Checks if dataset is consistent with the expected points"""
    if ds is None:
        msg.plain(f"SKIPPING, file not found: {url}")
        return False

    # always trust the first file that is read
    if expected_lon is None:
        return True

    if (ds[lon_str] == expected_lon).all() and (ds[lat_str] == expected_lat).all():
        return True
    else:
        msg.plain(f"SKIPPING, file inconsistent: {url}")
        return False
