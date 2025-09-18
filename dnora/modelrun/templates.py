# from ..run import model_executers
from dnora.modelrun.modelrun import ModelRun
from dnora.read.generic import ConstantData
from dnora.type_manager.dnora_types import DnoraDataType
import dnora
import warnings
import dnora.read.wind, dnora.read.spectra, dnora.read.spectra1d, dnora.read.waveseries, dnora.read.ice, dnora.read.waterlevel, dnora.read.current, dnora.read.grid

# DEPRECIATE = ["nora3", "era5", "noaa", "nchmf"]
DEPRECIATE = []


def deprecated_class(data_source: str, package: str):
    """
    Decorator to issue a deprecation warning whenever an instance of the class is used.
    """

    def decorator(cls):
        # Create a wrapper class that behaves like the original class
        classname = cls.__name__

        class Wrapper(cls):
            def __init__(self, *args, **kwargs):
                # Issue the warning when the class is instantiated
                message = (
                    f"{data_source} products have been moved to package dnora-{package} ('pip install dnora-{package}'). The 'dnora.modelrun.{classname}' class is deprecated and will be removed in a future release. "
                    f"Please use 'dnora_{package}.modelrun.{classname}' instead."
                )

                if package in DEPRECIATE:
                    warnings.warn(
                        message,
                        DeprecationWarning,
                        stacklevel=2,  # Points to the user's code
                    )
                super().__init__(*args, **kwargs)

        return Wrapper

    return decorator


class Constant(ModelRun):
    _reader_dict = {
        DnoraDataType.WIND: ConstantData(),
        DnoraDataType.SPECTRA: ConstantData(),
        DnoraDataType.SPECTRA1D: ConstantData(),
        DnoraDataType.ICE: ConstantData(),
        DnoraDataType.WATERLEVEL: ConstantData(),
        DnoraDataType.CURRENT: ConstantData(),
        DnoraDataType.WAVESERIES: ConstantData(),
    }


# @deprecated_class("MET Norway's", "metno")
class NORA3(ModelRun):
    _reader_dict = {
        DnoraDataType.SPECTRA: dnora.read.spectra.metno.NORA3(),
        DnoraDataType.WIND: dnora.read.wind.metno.NORA3(),
        DnoraDataType.ICE: dnora.read.ice.metno.NORA3(),
    }


# @deprecated_class("MET Norway's", "metno")
class WAM4km(ModelRun):
    _reader_dict = {
        DnoraDataType.SPECTRA: dnora.read.spectra.metno.WAM4km(),
        DnoraDataType.WIND: dnora.read.wind.metno.MEPS(),
    }


# @deprecated_class("MET Norway's", "metno")
class WW3_4km(ModelRun):
    _reader_dict = {
        DnoraDataType.SPECTRA: dnora.read.spectra.metno.WW3_4km(),
        DnoraDataType.WIND: dnora.read.wind.metno.MEPS(),
    }


# @deprecated_class("MET Norway's", "metno")
class CLIMAREST(ModelRun):
    _reader_dict = {
        DnoraDataType.SPECTRA: dnora.read.spectra.metno.CLIMAREST(),
        DnoraDataType.WIND: dnora.read.wind.metno.CLIMAREST(),
        DnoraDataType.ICE: dnora.read.ice.metno.CLIMAREST(),
    }


# @deprecated_class("ERA5", "era5")
class ERA5(ModelRun):
    _reader_dict = {
        DnoraDataType.SPECTRA: dnora.read.spectra.ec.ERA5(),
        DnoraDataType.WIND: dnora.read.wind.ec.ERA5(),
        DnoraDataType.WATERLEVEL: dnora.read.waterlevel.ec.GTSM_ERA5(),
    }


# @deprecated_class("NOAA", "noaa")
class PacIOOS(ModelRun):
    _reader_dict = {
        DnoraDataType.WIND: dnora.read.wind.noaa.PacIOOS(),
    }


# @deprecated_class("NOAA", "noaa")
class NCEP(ModelRun):
    _reader_dict = {
        DnoraDataType.WIND: dnora.read.wind.noaa.NCEP(),
    }


# @deprecated_class("NCHMF", "nchmf")
class NCHMF(ModelRun):
    _reader_dict = {
        # DnoraDataType.SPECTRA: dnora.read.spectra.nchmf.SWAN(),
        DnoraDataType.WIND: dnora.read.wind.nchmf.ECMWF(),
        DnoraDataType.WATERLEVEL: dnora.read.waterlevel.cmems.Global(),
    }
