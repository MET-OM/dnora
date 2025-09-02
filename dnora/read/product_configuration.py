from dataclasses import dataclass, field
from dnora.type_manager.spectral_conventions import SpectralConvention
from dnora.type_manager.data_sources import DataSource
from dnora.read.ds_read_functions import basic_xarray_read
from dnora.type_manager.dnora_types import DnoraDataType
from typing import Callable
from dnora.utils.io import get_url
from typing import Optional

def get_constant_url(folder, filename, file_times, **kwargs) -> list[str]:
    """Applies the same folder and filename to all file_times to get url.
    folder and file_name can contain %Y etc. that will be replaced"""
    return [get_url(folder, filename, file_time) for file_time in file_times]


"""ds_creator function will get called once with partial using the following arguments:

    lon=lon, #tuple
    lat=lat, #tuple
    data_type=obj_type, #DnoraDataType
    name=self.name(),
    program=program, # 'fimex'/'pyfimex', ignore if not needed

    It will get called with the following arguments when reading data:

    start_time, 
    end_time, 
    url # string to filename
    
    Needs to return an xr.Dataset"""


@dataclass
class ProductConfiguration:
    filename: str = "model_output_%Y%m.nc"
    default_filenames: dict[DataSource, str] = field(default_factory=dict)
    default_folders: dict[DataSource, str] = field(default_factory=dict)
    tile: Optional[str] = None
    tile_names: dict[str, str] = field(default_factory=dict)
    default_data_source: DataSource = DataSource.UNDEFINED
    convention: SpectralConvention = SpectralConvention.UNDEFINED
    ds_creator_function: callable = basic_xarray_read
    data_vars: list[str] = field(default_factory=list)
    ds_aliases: dict[str, str] = field(default_factory=dict)
    core_aliases: dict[str, str] = field(default_factory=dict)
    # only_vars: field(default_factory=dict)
    # ignore_vars: field(default_factory=dict)
    # dynamic: field(default_factory=dict)
    time_var: str = None
    url_function: Callable = field(default=get_constant_url)
    ds_pre_processor: Callable = field(default=(lambda ds: (ds, {})))

    # def get_core_aliases(self, obj_type):
    #     return self.core_aliases.get(obj_type)

    def get_default_filename(self, source):
        return self.default_filenames.get(source)
