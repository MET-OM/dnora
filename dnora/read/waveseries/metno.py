from dnora.type_manager.data_sources import DataSource
from dnora.read.waveseries import WW3Unstruct
import geo_parameters as gp
from dnora.read.abstract_readers import PointDataReader
from dnora.read.ds_read_functions import read_ds_list, read_first_ds
import xarray as xr
import numpy as np
from dnora import utils
from dnora.read.product_readers import SpectralProductReader
from dnora.read.product_configuration import ProductConfiguration
from functools import partial
from dnora.read.ds_read_functions import basic_xarray_read
from dnora.read.file_structure import FileStructure
import geo_parameters as gp
from dnora.read.aliases import WW3_DS_ALIASES
from dnora.read.depreciation_decorator import deprecated_class_call


@deprecated_class_call("MET Norway's", "metno", "waveseries")
class NORAC(SpectralProductReader):
    product_configuration = ProductConfiguration(
        filename="ww3.%Y%m.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/norac_wave/field",
            DataSource.INTERNAL: "sfiblues/wave_hindcast/hindcast_v2/field",
        },
        ds_creator_function=partial(basic_xarray_read, inds_var="node"),
        default_data_source=DataSource.REMOTE,
        ds_aliases=WW3_DS_ALIASES,
    )

    file_structure = FileStructure(
        stride="month",
    )


@deprecated_class_call("MET Norway's", "metno", "waveseries")
class E39(SpectralProductReader):
    product_configuration = ProductConfiguration(
        filename="%Y%m_E39_#TILENAME_wave.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/obs/buoy-svv-e39/%Y/%m",
        },
        tile_names={
            "A": "A_Sulafjorden",
            "B": "B_Sulafjorden",
            "B1": "B1_Sulafjorden",
            "C": "C_Sulafjorden",
            "C1": "C1_Sulafjorden",
            "D": "D_Breisundet",
            "F": "F_Vartdalsfjorden",
            "G": "G_Halsafjorden",
        },
        ds_creator_function=partial(basic_xarray_read),
        default_data_source=DataSource.REMOTE,
        core_aliases={gp.wave.Hs: "Hm0", gp.wave.Dirm: "mdir"},
    )

    file_structure = FileStructure(
        stride="month",
    )
