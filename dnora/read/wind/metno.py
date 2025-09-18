from copy import copy
import numpy as np
import xarray as xr
from subprocess import call
import os, glob
import pandas as pd
from functools import partial
from dnora.read.file_structure import FileStructure
import re
import geo_parameters as gp

# Import objects
from dnora.grid import Grid

# Import aux_funcsiliry functions
from dnora import msg
from dnora.utils.io import get_url
from dnora import utils

from dnora.type_manager.data_sources import DataSource
from dnora.type_manager.dnora_types import DnoraDataType
from dnora.read.abstract_readers import DataReader
from dnora.read.product_readers import ProductReader
from dnora.read.product_configuration import ProductConfiguration
from dnora.read.fimex_functions import ds_fimex_read
from dnora.read.ds_read_functions import read_ds_list, setup_temp_dir
import calendar
from dnora.read.fimex_functions import ds_fimex_read
from dnora.read.depreciation_decorator import deprecated_class_call


@deprecated_class_call("MET Norway's", "metno", "wind")
class NORA3(ProductReader):
    """Reads wind data (from monthly files 'arome3km_1hr_YYYYMM.nc') of the NORA3 hindcast directly from MET Norways servers.

    The NORA3 HARMONIE-AROME high-resolution (ca 3 km) hindcast for the
    North Sea, the Norwegian Sea, and the Barents Sea.

    Haakenstad, H., Breivik, Ø., Furevik, B. R., Reistad, M., Bohlinger, P., &
    Aarnes, O. J. (2021). NORA3: A Nonhydrostatic High-Resolution Hindcast of
    the North Sea, the Norwegian Sea, and the Barents Sea,
    Journal of Applied Meteorology and Climatology, 60(10), 1443-1464,
    DOI: 10.1175/JAMC-D-21-0029.1
    """

    product_configuration = ProductConfiguration(
        filename="arome3km_1hr_%Y%m.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/nora3_subset_atmos/atm_hourly_v2",
            DataSource.INTERNAL: "NORA3/atmosphere/atm_hourly",
        },
        ds_creator_function=partial(
            ds_fimex_read,
            resolution_in_km=3,
        ),
        data_vars=["wind_speed", "wind_direction"],
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure(stride="month")


@deprecated_class_call("MET Norway's", "metno", "wind")
class MyWave3km(ProductReader):
    """Reads wind data from the MyWave 3km hindcast directly from MET Norways
    servers. You should probably use NORA3 because:

    The wind data is from NORA3 (see the NORA3 reader for a data reference), is taken
    from the wave model output in this reader. This means that model land points have no data.
    """

    product_configuration = ProductConfiguration(
        filename="%Y%m%d_MyWam3km_hindcast.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/%Y/%m",
        },
        ds_creator_function=partial(
            ds_fimex_read,
            resolution_in_km=3,
        ),
        data_vars=["ff", "dd"],
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure(
        stride=24,
        hours_per_file=24,
        last_file="",
        lead_time=0,
    )


def get_meps_urls(folder: str, filename: str, file_times: list[str], **kwargs):
    """This is passed to the read_ds_list. We need it because the folder and filename that makes up th URL changes is time"""
    urls = []

    for file_time in file_times:
        if file_time >= pd.Timestamp("2020-02-04T12:00"):
            remote_filename = filename
        else:
            remote_filename = re.sub("det", "subset", filename)
        # Looks like oldarchive is gone
        # if file_time >= pd.Timestamp("2020-01-01T00:00"):
        #     remote_folder = folder
        # else:
        #     remote_folder = re.sub("meps25epsarchive", "mepsoldarchive", folder)
        remote_folder = folder
        urls.append(get_url(remote_folder, remote_filename, file_time))
    return urls


def meps_extra_fimex_commands(start_time, end_time, url) -> list[str]:
    """Determines the possible extra fimex commands needed to process the MEPS netcdf.
    Repend on the url ('det'/'subset')

    start_time and end_time accepted because of standard convention"""
    if "subset" in url:
        return [
            "--extract.reduceDimension.name=ensemble_member",
            "--extract.reduceDimension.start=1",
            "--extract.reduceDimension.end=1",
        ]
    else:
        return []


@deprecated_class_call("MET Norway's", "metno", "wind")
class MEPS(ProductReader):
    """Reads wind data from MET Norways MEPS forecast.

    The data is from a 2.5 km AROME model.
    """

    product_configuration = ProductConfiguration(
        filename=f"meps_det_2_5km_%Y%m%dT%HZ.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/meps25epsarchive/%Y/%m/%d",
        },
        ds_creator_function=partial(
            ds_fimex_read,
            resolution_in_km=2.5,
            extra_commands=meps_extra_fimex_commands,
        ),
        data_vars=["x_wind_10m", "y_wind_10m"],
        default_data_source=DataSource.REMOTE,
        url_function=get_meps_urls,
    )

    file_structure = FileStructure(
        stride=6,
        hours_per_file=67,
        last_file="",
        lead_time=0,
    )


def get_nora3fp_urls(
    folder: str,
    filename: str,
    file_times: list[str],
    start_times: list[str],
    lead_time: int,
) -> list[str]:
    """This is passed to the read_ds_list. We need it because the folder and filename that makes up th URL changes is time"""
    urls = []
    for file_time, start_time in zip(file_times, start_times):
        h0 = int(file_time.hour) % 6
        subfolder = file_time.strftime("%Y/%m/%d/") + (
            file_time - np.timedelta64(h0, "h")
        ).strftime("%H")

        remote_folder = re.sub("SUBFOLDER", subfolder, folder)

        first_ind = lead_time
        ind = int((start_time.hour - first_ind) % 6) + first_ind

        time_stamp = (
            file_time.strftime("%Y%m%d")
            + (file_time - np.timedelta64(h0, "h")).strftime("%H")
            + f"_{ind:03d}"
        )

        remote_filename = re.sub("TIMESTAMP", time_stamp, filename)

        urls.append(get_url(remote_folder, remote_filename))
    return urls


@deprecated_class_call("MET Norway's", "metno", "wind")
class NORA3_fp(ProductReader):
    """Reads wind data of the NORA3 hindcast directly from MET Norways servers.

    The NORA3 HARMONIE-AROME high-resolution (ca 3 km) hindcast for the
    North Sea, the Norwegian Sea, and the Barents Sea.

    Haakenstad, H., Breivik, Ø., Furevik, B. R., Reistad, M., Bohlinger, P., &
    Aarnes, O. J. (2021). NORA3: A Nonhydrostatic High-Resolution Hindcast of
    the North Sea, the Norwegian Sea, and the Barents Sea,
    Journal of Applied Meteorology and Climatology, 60(10), 1443-1464,
    DOI: 10.1175/JAMC-D-21-0029.1
    """

    product_configuration = ProductConfiguration(
        filename=f"fcTIMESTAMP_fp.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/nora3/SUBFOLDER",
        },
        ds_creator_function=partial(
            ds_fimex_read,
            resolution_in_km=3,
        ),
        data_vars=["wind_speed", "wind_direction"],
        default_data_source=DataSource.REMOTE,
        url_function=get_nora3fp_urls,
    )

    file_structure = FileStructure(
        stride=1, hours_per_file=1, last_file="", lead_time=4
    )


@deprecated_class_call("MET Norway's", "metno", "wind")
class CLIMAREST(ProductReader):
    product_configuration = ProductConfiguration(
        filename="wind_HCLIM43_MPIESM12LR_3hr_%Y_%m.nc",
        ds_creator_function=partial(
            ds_fimex_read,
            resolution_in_km=3.0,
        ),
        data_vars=["uas", "vas"],
        default_data_source=DataSource.LOCAL,
        ds_aliases={"uas": gp.wind.XWind, "vas": gp.wind.YWind},
    )

    file_structure = FileStructure(
        stride="month",
    )
