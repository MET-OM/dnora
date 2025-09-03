from dnora.type_manager.data_sources import DataSource
from dnora.read.waveseries import SWANnc
import geo_parameters as gp
from dnora.read.abstract_readers import PointDataReader
from dnora.read.ds_read_functions import read_ds_list, read_first_ds
from .waveseries_readers import ds_xarray_read
import xarray as xr
import numpy as np
from dnora import utils
from dnora.read.depreciation_decorator import deprecated_class_call


@deprecated_class_call("NCHMF", "nchmf", "spectra")
class SWAN4km(SWANnc):
    stride = 24  # int (for hourly), or 'month'
    hours_per_file = 84  # int (if not monthly files)

    _default_filename = "SWAN%Y%m%d%H.nc"

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE
