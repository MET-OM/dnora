from dataclasses import dataclass, field
from dnora.type_manager.spectral_conventions import SpectralConvention
from dnora.type_manager.data_sources import DataSource
from dnora.read.ds_read_functions import basic_xarray_read
from dnora.type_manager.dnora_types import DnoraDataType


@dataclass
class ProductConfiguration:
    filename: str = "model_output_%Y%m.nc"
    default_filenames: dict[DataSource, str] = field(default_factory=dict)
    default_folders: dict[DataSource, str] = field(default_factory=dict)
    tile: str | None = None
    tile_names: dict[str, str] = field(default_factory=dict)
    default_data_source: DataSource = DataSource.UNDEFINED
    convention: SpectralConvention = SpectralConvention.UNDEFINED
    ds_creator_function: callable = basic_xarray_read
    ds_aliases: dict[str, str] = field(default_factory=dict)
    core_aliases: dict[DnoraDataType, dict] = field(default_factory=dict)
    # only_vars: field(default_factory=dict)
    # ignore_vars: field(default_factory=dict)
    # dynamic: field(default_factory=dict)
    time_var: str = None

    def get_core_aliases(self, obj_type):
        return self.core_aliases.get(obj_type)

    def get_default_filename(self, source):
        return self.default_filenames.get(source)
