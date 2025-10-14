import xarray as xr
import pandas as pd
import numpy as np
from .file_structure import FileStructure
import os, glob
from dnora.utils.io import get_url
from dnora import msg
from dnora.type_manager.dnora_types import DnoraDataType

from typing import Callable
import geo_parameters as gp
import ftplib

import os


def get_rc_credentials(rc_file: str) -> tuple[str]:

    # Define the path to your configuration file
    config_path = os.path.expanduser(f"~/{rc_file}")

    # Initialize a dictionary to store the credentials
    config = {}

    try:
        with open(config_path, "r") as file:
            for line in file:
                # Strip whitespace and split the line into key and value
                key, value = line.strip().split(":", 1)
                config[key.strip()] = value.strip()

        # Retrieve the username and password
        ftp_username = config["ftp_username"]
        ftp_password = config["ftp_password"]

    except FileNotFoundError:
        print(f"Configuration file not found at {config_path}")

    return ftp_username, ftp_password


def ftp_read(
    start_time: pd.Timestamp,
    end_time: pd.Timestamp,
    url: str,
    ## Partial variables from ProductReader
    lon: tuple[float],
    lat: tuple[float],
    name: str,
    data_type: DnoraDataType,
    data_vars: list[str],
    ## Partial variables in ProductConfiguration
    **kwargs,
):
    temp_dir = setup_temp_dir(data_type, name)

    username, password = get_rc_credentials(".nchmf_ecmwfrc")
    ftp = ftplib.FTP("fog.met.no")
    ftp.login(username, password)
    ftp.cwd("/")

    local_file = f"{temp_dir}/{name}_{url}"
    with open(local_file, "wb") as lf:
        ftp.retrbinary(f"RETR {url}", lf.write)
    ds = xr.open_dataset(local_file)[data_vars]
    return ds.sel(time=slice(start_time, end_time))


def find_time_var_in_ds(ds):
    """Tries to identify the time variable in a dataset"""
    if "time" in list(ds.coords):
        return "time"
    time_vars = [v for v in list(ds.coords) if "time" in v]
    if not time_vars:
        raise KeyError(f"Cant find a time variable in {list(ds.coords)}")
    else:
        return time_vars[0]


def basic_xarray_read(
    start_time: pd.Timestamp,
    end_time: pd.Timestamp,
    url: str,
    time_var: str = None,
    dt_for_time_stamps_in_hours: float = 0,
    lon: tuple[float] = None,
    lat: tuple[float] = None,
    inds: list[int] = None,
    inds_var: str = "inds",
    data_vars: list[str] = None,
    chunks=None,
    time_reference_str: str = "",
    **kwargs,
):
    """Reads a single file into an Xarray Dataset:

    1) The file is vut between start_time and end_time (variable 'time')
        - If there is no variable 'time', then first variable containing string 'time' is used
        - Automatic detection can be overridden by providing keyword time_var = 'xxxx'
    2) If dt_for_time_stamps_in_hours > 0, then no attempt to decode times in the netcdf is made
        - We will instead create time stamps with given time step
    3) If lon and lat is given, then the dataset is cut to those limits
    4) If inds is given, then dataset is cut to those indeces
        - If indexd variable is not 'inds', it can be specified with e.g. inds='station'
    5) If data_vars is given, then dataset is cut to those data variables
    """
    with xr.open_dataset(
        url,
        decode_times=(not dt_for_time_stamps_in_hours and not time_reference_str),
        chunks=chunks,
    ) as f:
        time_var = time_var or find_time_var_in_ds(f)
        if dt_for_time_stamps_in_hours:
            times = pd.date_range(
                start_time, end_time, freq=f"{dt_for_time_stamps_in_hours}h"
            )
            f[time_var] = times
        elif time_reference_str:
            f = f.assign_attrs({"units": time_reference_str})
            f[time_var] = f[time_var].assign_attrs({"units": time_reference_str})
            f = xr.decode_cf(f)

        # We want exact minutes to slice accurately
        f[time_var] = f[time_var].dt.round("min")
        ds = f.sel(**{time_var: slice(start_time, end_time)})

        # Make sure longitude is between -180 and 180 (E.g. GFS has 0 to 360)
        lon_str = gp.grid.Lon.find_me_in_ds(ds, return_first=True)
        if not lon_str:
            lon_str = "lon" if "lon" in ds.coords else "longitude"
        lat_str = gp.grid.Lat.find_me_in_ds(ds, return_first=True)
        if not lat_str:
            lat_str = "lat" if "lat" in ds.coords else "latitude"
        if np.where(ds[lon_str].data > 180)[0].size > 0:
            ds = ds.assign_coords(
                **{
                    lon_str: np.where(
                        ds[lon_str].data > 180, ds[lon_str].data - 360, ds[lon_str].data
                    )
                },
            )
            ds = ds.sortby(lon_str)

        if lon is not None and lat is not None:
            if len([c for c in ds.get(lat_str).shape if c > 1]) == 1:
                ds = ds.sel(**{lat_str: slice(*lat)})
            if len([c for c in ds.get(lon_str).shape if c > 1]) == 1:
                ds = ds.sel(**{lon_str: slice(*lon)})

        if inds is not None and hasattr(ds, inds_var):
            ds = ds.isel(**{inds_var: inds})

        if data_vars:
            ds = ds[data_vars]

    return ds


def read_list_of_spatial_ds(folder: str, filename: str):
    """Reads in a list of Datasets when different points are scattered over different files.
    Assumes all files cover the entire time period"""

    url = get_url(folder, filename)
    files = glob.glob(url)
    if not files:
        raise FileNotFoundError(f"Cannot find any files {url}!")
    ds_list = []
    for file in files:
        ds = xr.open_dataset(file)
        ds_list.append(ds)

    return ds_list


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

        except (
            OSError,
            RuntimeError,
            ftplib.error_perm,
        ) as e:  # xr gives OSError, fimex gives RuntimeError
            msg.plain(f'SKIPPING! Got error "{e}" while reading {url}')
            try_next_file = True

        if try_next_file:
            if data_left_to_try_with(hours_per_file, n, ct, file_times, end_time):
                ct += 1
            else:
                msg.plain(f"SKIPPING, no data left to try with!")
                ds = None
                keep_trying = False

    return ds, url


def read_ds_list(
    start_times: pd.DatetimeIndex,
    end_times: pd.DatetimeIndex,
    file_times: pd.DatetimeIndex,
    folder: str,
    filename: str,
    ds_creator_function: Callable = basic_xarray_read,
    url_function: Callable = None,
    hours_per_file: int = None,
    lead_time: int = None,
) -> list[xr.Dataset]:
    """Reads a list of xr.Datasets using the time stamps and folder/filename given.

    If one file is missing, the function tries to patch if the files have overlap (if hours_per_file is given).
    the ds_creator function takes arguments (start_time, end_time, url) and returns an xr.Dataset.
    """
    urls = url_function(
        folder, filename, file_times, start_times=start_times, lead_time=lead_time
    )
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
                expected_shape = tuple(
                    [ds[var].size for var in keys if "time" not in var]
                )
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

    # Take ct-1 since we want to check what the file contains if we increate ct with one in the next loop
    end_of_previous_file = pd.Timestamp(file_times[n - ct - 1]) + pd.Timedelta(
        hours_per_file - 1, "hours"
    )
    last_time_to_be_read = pd.Timestamp(end_time)

    if last_time_to_be_read > end_of_previous_file:
        # if pd.Timestamp(end_time) - pd.Timestamp(file_times[n - ct - 1]) > pd.Timedelta(
        #    hours_per_file, "hours"
        # ):
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
    given_shape = tuple([ds[var].size for var in keys if "time" not in var])
    if given_shape == expected_shape:
        return True
    else:
        return False


def setup_temp_dir(
    data_type: DnoraDataType, reader_name: str, clean_old_files: bool = True
) -> str:
    """Sets up a temporery directory for fimex files and cleans out possible old files"""
    temp_folder = f"dnora_{data_type.name.lower()}_temp"
    if not os.path.isdir(temp_folder):
        os.mkdir(temp_folder)
        print("Creating folder %s..." % temp_folder)
    if clean_old_files:
        msg.plain("Removing old files from temporary folder...")
        for f in glob.glob(f"dnora_{data_type.name.lower()}_temp/{reader_name}*.*"):
            os.remove(f)

    return temp_folder
