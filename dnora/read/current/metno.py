from copy import copy
import os
from dnora.grid import Grid

# Import abstract classes
from dnora.type_manager.data_sources import DataSource
from dnora.read.file_structure import FileStructure

# Import aux_funcsiliry functions
from dnora.utils.io import get_url
from dnora import utils

from functools import partial
from dnora.read.fimex_functions import ds_fimex_read
from dnora.read.ds_read_functions import basic_xarray_read
from dnora.read.product_readers import ProductReader
from dnora.read.product_configuration import ProductConfiguration
import pandas as pd
import re
from dnora.process.gridded import FillNaNs
import geo_parameters as gp
from dnora.read.depreciation_decorator import deprecated_class_call


def get_norkyst800_urls(folder: str, filename: str, file_times: list[str], **kwargs):
    """This is passed to the read_ds_list. We need it because the folder and filename that makes up th URL changes is time"""
    urls = []

    for file_time in file_times:
        if file_time >= pd.Timestamp("2018-01-01T00:00"):
            remote_folder = folder
        else:
            remote_folder = re.sub(
                "fou-hi/norkyst800m-1h", "sea/norkyst800mv0_1h/", folder
            )
        urls.append(get_url(remote_folder, filename, file_time))
    return urls


@deprecated_class_call("MET Norway's", "metno", "current")
class NorKyst800(ProductReader):
    """Reads ocean_current data of the NorKyst800 archieve directly from MET Norways servers.

    NorKyst-800 (Norwegian Coast 800m) is a numerical, high-resolution, ocean modelling
    system covering the Norwegian Coast.

    Albretsen, J., Sperrevik, A.K., Staalstrøm, A., Sandvik, A.D., Vikebø, F., Asplin, L., 2011.
    NorKyst-800 Rapport nr. 1: Brukermanual og tekniske beskrivelser. NorKyst-800 Report
    No. 1: User Manual and technical descriptions.
    """

    product_configuration = ProductConfiguration(
        filename="NorKyst-800m_ZDEPTHS_his.an.%Y%m%d00.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h",
        },
        ds_creator_function=partial(
            ds_fimex_read,
            resolution_in_km=0.8,
        ),
        data_vars=["u", "v"],
        default_data_source=DataSource.REMOTE,
        url_function=get_norkyst800_urls,
        ds_aliases={"u": gp.ocean.XCurrent, "v": gp.ocean.YCurrent},
        ds_pre_processor=lambda ds: (ds.isel(depth=0), {}),
    )

    file_structure = FileStructure(
        stride=24,
        hours_per_file=24,
    )

    def post_processing(self):
        return FillNaNs(0)


@deprecated_class_call("MET Norway's", "metno", "current")
class NorFjords160(ProductReader):
    """ """

    product_configuration = ProductConfiguration(
        filename="norfjords_160m_his_%Y%m%d01_surface_interp.nc",
        default_folders={
            DataSource.INTERNAL: "SWAN/Bjornafjorden2/ROMS/",
        },
        default_data_source=DataSource.INTERNAL,
        ds_creator_function=basic_xarray_read,
        time_var="ocean_time",
        ds_aliases={"u": gp.ocean.XCurrent, "v": gp.ocean.YCurrent},
        ds_pre_processor=lambda ds: (
            ds,
            {"lon": ds.lon.values[:, 0], "lat": ds.lat.values[0, :]},
        ),
    )

    file_structure = FileStructure(stride=24, hours_per_file=24, offset=1)

    def post_processing(self):
        return FillNaNs(0)
