# Import abstract classes and needed instances of them
from dnora.process.spectra import RemoveEmpty
from dnora.type_manager.spectral_conventions import SpectralConvention

# Import aux_funcsiliry functions
from dnora.type_manager.data_sources import DataSource
from dnora.read.spectra import SWAN as SWANBase
from dnora.read.depreciation_decorator import deprecated_class_call


@deprecated_class_call("NCHMF", "nchmf", "spectra")
class SWAN(SWANBase):
    """covers the East Sea"""

    stride: int = 24
    hours_per_file: int = 73
    offset: int = 12

    def post_processing(self):
        return RemoveEmpty()
